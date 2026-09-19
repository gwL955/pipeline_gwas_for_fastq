#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""★ 资源探测与规划：CPU / 内存实测（含 cgroup / WSL2 配额识别），
workers × 每样本线程 × sort/GATK 内存推导。覆盖 16 线程/20G ~ 100 线程/900G 及以上。

推导链（auto 档）：
  探测 cpu/mem → 预留（2 核 + max(2G, 5% 内存)）→ 按核数分档定每样本线程 T 与
  workers(核) → 按内存收紧 workers（单样本峰值 = 索引 17G + T×sort缓冲 + GATK + 杂项）
  → 反推 SORT_MEM（128M-2G 钳位）与 GATK_MEM（1g-8g 钳位）
  → 按步骤类型分化 workers（DEC-35）：比对类沿用 workers（bwa 索引峰值口径）；
  GATK 单线程类 workers_gatk（每路峰值 1.3×gatk 堆、不含 bwa 索引）；
  IO 类 workers_io（fastp/fastqc/md5/合并，每路 ~2G）；HC hmm 线程反推保证
  workers_gatk×hmm ≤ usable_cores（防超订阅）
  → 快速失败：workers=1 时峰值仍 > 可用内存则报错退出（防跑到一半 OOM）。
"""

import os
import math
import json

import config

# ── 步骤类型分化的推导常量（DEC-35，非告警阈值故不入 config TH 表） ──────
GATK_MEM_OVERHEAD = 1.3   # GATK 每 worker 实际峰值 ≈ 堆×1.3（JVM 元空间/GC/IO 开销）
IO_WORKER_PEAK_GB = 2.0   # IO 类每路峰值：fastqc JVM ~1G + fastp 缓冲（FASTP_BUF_GB）


# ── 探测 ────────────────────────────────────────────────────────────────
def _cgroup_paths():
    """本进程 cgroup v2 相对路径链（自身 → 祖先 → 根），用于逐层取最紧限额"""
    paths = []
    try:
        with open("/proc/self/cgroup") as f:
            for line in f:
                parts = line.strip().split(":", 2)
                if len(parts) == 3 and parts[2]:
                    rel = parts[2].lstrip("/")
                    cur = rel
                    while cur:
                        if cur not in paths:
                            paths.append(cur)
                        cur = cur.rsplit("/", 1)[0] if "/" in cur else ""
    except OSError:
        pass
    return paths


def detect_cpu():
    """可用核 = min(os.cpu_count, sched_getaffinity(taskset 场景), cgroup cpu.max 配额)
    （WSL2/容器场景必须识别，不得只看宿主机核数）"""
    cpu = os.cpu_count() or 1
    try:
        cpu = min(cpu, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    base = "/sys/fs/cgroup"
    for rel in _cgroup_paths() + [""]:
        # cgroup v2
        f = os.path.join(base, rel, "cpu.max")
        if os.path.isfile(f):
            try:
                with open(f) as fh:
                    parts = fh.read().split()
                if len(parts) == 2 and parts[0] != "max":
                    quota = int(parts[0]) / int(parts[1])
                    cpu = min(cpu, quota)
            except (OSError, ValueError):
                pass
        # cgroup v1
        fq = os.path.join(base, "cpu", rel, "cpu.cfs_quota_us")
        fp = os.path.join(base, "cpu", rel, "cpu.cfs_period_us")
        if os.path.isfile(fq) and os.path.isfile(fp):
            try:
                q, p = int(open(fq).read()), int(open(fp).read())
                if q > 0:
                    cpu = min(cpu, q / p)
            except (OSError, ValueError):
                pass
    return max(1, int(math.floor(cpu)))


def detect_mem_gb():
    """可用内存 GB = min(MemAvailable, cgroup memory 限额)（含 v1/v2）"""
    avail = None
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    avail = int(line.split()[1]) / 1e6   # kB → GB(十进制)
                    break
    except OSError:
        pass
    if avail is None:
        avail = 16.0
    base = "/sys/fs/cgroup"
    for rel in _cgroup_paths() + [""]:
        f = os.path.join(base, rel, "memory.max")            # v2
        if os.path.isfile(f):
            try:
                v = open(f).read().strip()
                if v != "max":
                    avail = min(avail, int(v) / 1e9)
            except (OSError, ValueError):
                pass
        f = os.path.join(base, "memory", rel, "memory.limit_in_bytes")  # v1
        if os.path.isfile(f):
            try:
                v = int(open(f).read())
                if 0 < v < (1 << 60):
                    avail = min(avail, v / 1e9)
            except (OSError, ValueError):
                pass
    return round(avail, 1)


# ── 规划 ────────────────────────────────────────────────────────────────
class ResourcePlan:
    def __init__(self, **kw):
        self.__dict__.update(kw)

    def to_dict(self):
        d = dict(self.__dict__)
        for k in ("cpu_detected", "mem_detected"):
            if k in d and not isinstance(d[k], (int, float)):
                d[k] = str(d[k])
        return d

    def table(self):
        """计划透明：探测值→预留→各类 workers/T/SORT_MEM/GATK_MEM 完整推导表"""
        lines = [
            "─── 资源计划（resource plan）───",
            f"  profile           : {self.profile}",
            f"  探测 CPU/内存     : {self.cpu_detected} 线程 / {self.mem_detected} GB",
            f"  系统预留          : {self.reserve_cores} 核 + {self.reserve_mem_gb} GB",
            f"  可用(扣预留)      : {self.usable_cores} 线程 / {self.usable_mem_gb} GB",
            f"  每样本线程 T      : {self.threads}",
            f"  workers 比对类    : 核 {self.workers_cpu} ; 内存 {self.workers_mem}"
            f" → 取较小 = {self.workers}（Step2 bwa，索引 17G 峰值口径）",
            f"  workers GATK 类   : 核 {self.workers_gatk_cpu} ; 内存 {self.workers_gatk_mem}"
            f" → 取较小 = {self.workers_gatk}（Step3/4/5-HC/6，"
            f"每路 ~{self.gatk_peak_gb} GB = {GATK_MEM_OVERHEAD}×{self.gatk_mem}，不含 bwa 索引）",
            f"  workers IO 类     : {self.workers_io}（Step0/1 md5·合并·fastp·fastqc，"
            f"每路 ~{IO_WORKER_PEAK_GB} GB，≤ GATK 类）",
            f"  sort -m(每线程)   : {self.sort_mem}",
            f"  GATK -Xmx         : {self.gatk_mem}",
            f"  bwa/fastp/fastqc  : -t {self.threads} / --thread {self.fastp_threads} / -t {self.fastqc_threads}",
            f"  HC pair-hmm 线程  : {self.hc_hmm_threads}"
            f"（GATK 类 {self.workers_gatk}×{self.hc_hmm_threads} ≤ 可用 {self.usable_cores} 核，防超订阅）"
            f" ; mosdepth -t {self.mosdepth_threads}",
            f"  cohort 步骤 -Xmx  : {self.cohort_mem}（串行，可用内存 60% 钳位）",
            f"  单样本峰值内存    : {self.peak_per_sample_gb} GB（bwa 索引 {config.BWA_INDEX_MEM_GB}G"
            f" + {self.threads}×{self.sort_mem} + {self.gatk_mem} + 杂项 {config.FASTP_BUF_GB}G）",
            f"  预计总占用        : 比对 {self.workers} × {self.peak_per_sample_gb} = "
            f"{round(self.workers * self.peak_per_sample_gb, 1)} GB ; GATK 类 {self.workers_gatk} × "
            f"{self.gatk_peak_gb} = {round(self.workers_gatk * self.gatk_peak_gb, 1)} GB"
            f" / 可用 {self.mem_detected} GB",
        ]
        return "\n".join(lines)


def _threads_by_cores(cores):
    if cores >= 64:
        return 24
    if cores >= 32:
        return 12
    if cores >= 16:
        return 8
    return 4


def plan(profile="auto", workers_override=None, threads_override=None,
         max_memory_gb=None):
    """生成资源计划。优先级：--workers > --threads/--max-memory > --resource-profile > 自动探测"""
    cpu_det, mem_det = detect_cpu(), detect_mem_gb()

    if profile == "low":
        cpu_det, mem_det = config.PROFILE_LOW["cpu"], config.PROFILE_LOW["mem_gb"]
    elif profile == "high":
        cpu_det, mem_det = config.PROFILE_HIGH["cpu"], config.PROFILE_HIGH["mem_gb"]
    # auto：实测硬件，超过 100 线程/900G 照常使用（不设上限），只扣系统预留
    if max_memory_gb:
        mem_det = min(mem_det, float(max_memory_gb))

    reserve_cores = min(config.RESERVE_CORES, max(0, cpu_det - 1))
    reserve_mem = max(2.0, round(0.05 * mem_det, 1))
    usable_cores = max(1, cpu_det - reserve_cores)
    usable_mem = round(max(1.0, mem_det - reserve_mem), 1)

    T = _threads_by_cores(usable_cores)
    if threads_override:
        T = max(1, int(threads_override))

    # 内存收紧：先按最小参数估峰值，求 workers(内存)
    peak_min = config.BWA_INDEX_MEM_GB + T * (config.SORT_MEM_MIN / 1024) \
        + config.GATK_MEM_MIN_GB + config.FASTP_BUF_GB
    workers_cpu = max(1, usable_cores // T)
    workers_mem = max(1, int(usable_mem // peak_min))
    workers = min(workers_cpu, workers_mem)

    # 反推内存参数：按每 worker 预算分配（bwa 索引固定 → GATK 取剩余的 1/4 →
    # sort 只吃剩余的 1/4——实测 bwa-mem2 运行时开销+页缓存 ≈ 索引之外 2-4G，
    # 留足冗余防 cgroup/裸机 OOM；大预算被 2G 钳位不影响高配）
    budget = round(usable_mem / workers, 1)
    gatk_gb = int(round((budget - config.BWA_INDEX_MEM_GB - config.FASTP_BUF_GB) / 4))
    gatk_gb = max(config.GATK_MEM_MIN_GB, min(config.GATK_MEM_MAX_GB, gatk_gb))
    sort_mb = int(round((budget - config.BWA_INDEX_MEM_GB - config.FASTP_BUF_GB - gatk_gb)
                        / T * 1024 * 0.25))
    sort_mb = max(config.SORT_MEM_MIN, min(config.SORT_MEM_MAX, sort_mb))

    peak = round(config.BWA_INDEX_MEM_GB + T * sort_mb / 1024 + gatk_gb + config.FASTP_BUF_GB, 1)

    # ── 按步骤类型分化 workers（DEC-35）────────────────────────────────
    # 实测瓶颈：GATK 单线程工具（MarkDuplicates/BQSR/HC/HsMetrics）以比对类
    # workers 并行时整机 CPU ~10%（比对类 workers 被每路 17G bwa 索引峰值钉死）；
    # GATK 阶段索引不在工作集，每路峰值 ≈ 1.3×gatk 堆。
    # 比对类（Step2）：沿用 workers（bwa 饱和设计，96% CPU，勿动）；
    # GATK 类（Step3/4/5-HC/6）：CPU 每路预算 max(2, hmm 档) 核（hmm 线程 +
    #   JVM GC/IO 余量）——50 核机 48//4=12 路 × hmm4 = 48 核满载不超订；
    # IO 类（Step0 md5/合并、Step1 fastp/fastqc）：每路 ~2G，至多与 GATK 类同路
    #   （收益递减且压缩线程挤占）。
    hmm_tier = min(4, max(1, T // 2))          # 期望档：与现行口径一致（T=12→4，T=4→2）
    gatk_peak = round(GATK_MEM_OVERHEAD * gatk_gb, 1)
    workers_gatk_cpu = max(1, usable_cores // max(2, hmm_tier))
    workers_gatk_mem = max(1, int(usable_mem // gatk_peak))
    workers_gatk = min(workers_gatk_cpu, workers_gatk_mem)
    workers_io = max(1, min(int(usable_mem // IO_WORKER_PEAK_GB), workers_gatk))

    if workers_override:
        workers = workers_gatk = workers_io = max(1, int(workers_override))
    # HC hmm 反推（防超订阅）：workers_gatk×hmm ≤ usable_cores；档位不超过期望档
    hc_hmm = min(hmm_tier, max(1, usable_cores // workers_gatk))

    # 快速失败：单样本峰值装不下（如机器不足 20G 装不下 bwa-mem2 人类索引）
    if peak > mem_det:
        raise SystemExit(
            f"[FATAL] 资源规划失败：单样本峰值内存 {peak} GB > 可用内存 {mem_det} GB"
            f"（bwa-mem2 索引常驻 {config.BWA_INDEX_MEM_GB} GB）。"
            f"请换更大内存机器或使用不含比对的步骤（--step <2）。"
        )
    # 内存方向再钳一次 workers（用户覆盖则只告警不强制）
    if workers * peak > mem_det:
        if workers_override:
            warn = (f"[WARN] --workers={workers} 超出内存预算：{workers}×{peak} GB > "
                    f"{mem_det} GB，存在 OOM 风险")
            print(warn, flush=True)
        else:
            workers = max(1, int(mem_det // peak))
    # GATK 类同口径兜底（仅用户覆盖可能越界；自动路径已由 min 收紧）
    if workers_gatk * gatk_peak > mem_det and workers_override:
        print(f"[WARN] --workers={workers_gatk} 的 GATK 类超出内存预算："
              f"{workers_gatk}×{gatk_peak} GB > {mem_det} GB，存在 OOM 风险", flush=True)

    cohort_gb = int(round(0.6 * mem_det))
    cohort_gb = max(config.COHORT_MEM_MIN_GB, min(config.COHORT_MEM_MAX_GB, cohort_gb))

    return ResourcePlan(
        profile=profile,
        cpu_detected=cpu_det, mem_detected=mem_det,
        reserve_cores=reserve_cores, reserve_mem_gb=reserve_mem,
        usable_cores=usable_cores, usable_mem_gb=usable_mem,
        threads=T, workers=workers,
        workers_cpu=workers_cpu, workers_mem=workers_mem,
        workers_gatk=workers_gatk, workers_io=workers_io,
        workers_gatk_cpu=workers_gatk_cpu, workers_gatk_mem=workers_gatk_mem,
        gatk_peak_gb=gatk_peak,
        sort_mem=f"{sort_mb}M",   # samtools 只接受整数+单位（2.0G 会被误解析为 2 字节）
        gatk_mem=f"{gatk_gb}g",
        fastp_threads=min(16, T), fastqc_threads=min(8, max(2, T // 3)),
        hc_hmm_threads=hc_hmm, mosdepth_threads=min(8, max(2, T // 3)),
        cohort_mem=f"{cohort_gb}g",
        peak_per_sample_gb=peak,
    )


def plan_json(p):
    return json.dumps(p.to_dict(), ensure_ascii=False, indent=2)
