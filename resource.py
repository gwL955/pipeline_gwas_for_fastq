#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""★ 资源探测与规划：CPU / 内存实测（含 cgroup / WSL2 配额识别），
workers × 每样本线程 × sort/GATK 内存推导。覆盖 16 线程/20G ~ 100 线程/900G 及以上。

推导链（auto 档）：
  探测 cpu/mem → 预留（2 核 + max(2G, 5% 内存)）→ 按核数分档定每样本线程 T 与
  workers(核) → 按内存收紧 workers（单样本峰值 = 索引 17G + T×sort缓冲 + GATK + 杂项）
  → 反推 SORT_MEM（128M-2G 钳位）与 GATK_MEM（1g-8g 钳位）
  → 快速失败：workers=1 时峰值仍 > 可用内存则报错退出（防跑到一半 OOM）。
"""

import os
import math
import json

import config


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
        """计划透明：探测值→预留→workers/T/SORT_MEM/GATK_MEM 完整推导表"""
        lines = [
            "─── 资源计划（resource plan）───",
            f"  profile           : {self.profile}",
            f"  探测 CPU/内存     : {self.cpu_detected} 线程 / {self.mem_detected} GB",
            f"  系统预留          : {self.reserve_cores} 核 + {self.reserve_mem_gb} GB",
            f"  可用(扣预留)      : {self.usable_cores} 线程 / {self.usable_mem_gb} GB",
            f"  每样本线程 T      : {self.threads}",
            f"  workers(核推导)   : {self.workers_cpu} ; workers(内存推导): {self.workers_mem}"
            f" → 取较小 = {self.workers}",
            f"  sort -m(每线程)   : {self.sort_mem}",
            f"  GATK -Xmx         : {self.gatk_mem}",
            f"  bwa/fastp/fastqc  : -t {self.threads} / --thread {self.fastp_threads} / -t {self.fastqc_threads}",
            f"  HC pair-hmm 线程  : {self.hc_hmm_threads} ; mosdepth -t {self.mosdepth_threads}",
            f"  cohort 步骤 -Xmx  : {self.cohort_mem}（串行，可用内存 60% 钳位）",
            f"  单样本峰值内存    : {self.peak_per_sample_gb} GB（bwa 索引 {config.BWA_INDEX_MEM_GB}G"
            f" + {self.threads}×{self.sort_mem} + {self.gatk_mem} + 杂项 {config.FASTP_BUF_GB}G）",
            f"  预计总占用        : {self.workers} × {self.peak_per_sample_gb} = "
            f"{round(self.workers * self.peak_per_sample_gb, 1)} GB / 可用 {self.mem_detected} GB",
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

    if workers_override:
        workers = max(1, int(workers_override))

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

    cohort_gb = int(round(0.6 * mem_det))
    cohort_gb = max(config.COHORT_MEM_MIN_GB, min(config.COHORT_MEM_MAX_GB, cohort_gb))

    return ResourcePlan(
        profile=profile,
        cpu_detected=cpu_det, mem_detected=mem_det,
        reserve_cores=reserve_cores, reserve_mem_gb=reserve_mem,
        usable_cores=usable_cores, usable_mem_gb=usable_mem,
        threads=T, workers=workers,
        workers_cpu=workers_cpu, workers_mem=workers_mem,
        sort_mem=f"{sort_mb}M",   # samtools 只接受整数+单位（2.0G 会被误解析为 2 字节）
        gatk_mem=f"{gatk_gb}g",
        fastp_threads=min(16, T), fastqc_threads=min(8, max(2, T // 3)),
        hc_hmm_threads=min(4, max(1, T // 2)), mosdepth_threads=min(8, max(2, T // 3)),
        cohort_mem=f"{cohort_gb}g",
        peak_per_sample_gb=peak,
    )


def plan_json(p):
    return json.dumps(p.to_dict(), ensure_ascii=False, indent=2)
