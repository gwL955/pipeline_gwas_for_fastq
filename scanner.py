#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""输入布局扫描（三种布局，每式一个独立识别器函数）+ md5 校验 +
Lane 合并 + samples.tsv 生成。批次是最小分析单元，本模块按单个批次工作。"""

import os
import re
import hashlib
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import config
from runner import nonempty

ILLUMINA_RE = re.compile(
    r"^(?P<sample>.+)_S(?P<snum>\d+)_L(?P<lane>\d{3})_R(?P<read>[12])_001\.fastq\.gz$")
OUTSOURCED_RE = re.compile(r"_R(?P<read>[12])\.fastq\.gz$")
# 外送平铺整名锚定；与 ILLUMINA_RE 互斥（结尾 _R#.fastq.gz 与 _001.fastq.gz 不可能兼得）
FLAT_OUTSOURCED_RE = re.compile(r"^(?P<sample>.+)_R(?P<read>[12])\.fastq\.gz$")

# 布局中文名（冲突提示/samples.tsv note 列共用）
LAYOUT_LABELS = {"illumina": "平铺", "outsourced": "外送",
                 "outsourced_flat": "外送平铺"}

# 布局识别后直接剔除的非样本条目（DEC-31，样本名不区分大小写）：Illumina 下机
# 自带 Undetermined（BCLConvert 未匹配 index 的 reads，常为千万级）不是实验
# 样本，不得进入分析——曾一路进联合分型/基因型矩阵/Output 污染整批（RUN-45）
IGNORED_SAMPLES = {"undetermined"}


class SampleInfo:
    def __init__(self, sample, layout):
        self.sample = sample
        self.layout = layout          # "illumina" / "outsourced" / "outsourced_flat"
        self.r1 = []                  # 宿主绝对路径列表（按 Lane 排序）
        self.r2 = []
        self.snum = None              # Illumina S 号（追溯用）
        self.lanes = []               # Illumina Lane 号列表
        self.invalid_reason = None

    @property
    def valid(self):
        return self.invalid_reason is None

    def note(self):
        if self.layout == "illumina":
            return f"{self.snum};Lane={'-'.join(self.lanes)}"
        return f"{LAYOUT_LABELS.get(self.layout, self.layout)};{len(self.r1)}文件"


# ── 布局识别器（DEC-23：每式一个独立函数并登记于 LAYOUT_SCANNERS，
#    新增输入格式只追加函数，不动 scan_batch 主体）───────────────────────
# 样本名来源：Illumina 平铺=文件名去 _S#_L###_R[12]_001 尾；外送子目录=子目录名；
# 外送平铺=文件名去 _R[12] 尾（RUN-36 新增——外送交付常直接平铺于批次目录）。


def scan_illumina_flat(flat, subdirs):
    """布局 1：Illumina 平铺（BCLConvert 命名，多 Lane）"""
    samples = {}
    for fp in flat:
        m = ILLUMINA_RE.match(os.path.basename(fp))
        if not m:
            continue
        sm = m.group("sample")
        si = samples.setdefault(sm, SampleInfo(sm, "illumina"))
        if si.snum is None:
            si.snum = f"S{m.group('snum')}"
        if m.group("lane") not in si.lanes:
            si.lanes.append(m.group("lane"))
        (si.r1 if m.group("read") == "1" else si.r2).append(fp)
    return samples


def scan_outsourced_subdir(flat, subdirs):
    """布局 2：外送子目录 <样本名>/<样本名>_R1.fastq.gz（样本名=子目录名）"""
    samples = {}
    for dp in subdirs:
        r1s, r2s = [], []
        for name in sorted(os.listdir(dp)):
            m = OUTSOURCED_RE.search(name)
            if m:
                (r1s if m.group("read") == "1" else r2s).append(os.path.join(dp, name))
        if not r1s and not r2s:
            continue
        si = SampleInfo(os.path.basename(dp), "outsourced")
        si.r1, si.r2 = r1s, r2s
        samples[si.sample] = si
    return samples


def scan_outsourced_flat(flat, subdirs):
    """布局 3：外送平铺 <样本名>_R1.fastq.gz（文件直接放批次目录）"""
    samples = {}
    for fp in flat:
        m = FLAT_OUTSOURCED_RE.match(os.path.basename(fp))
        if not m:
            continue
        sm = m.group("sample")
        si = samples.setdefault(sm, SampleInfo(sm, "outsourced_flat"))
        (si.r1 if m.group("read") == "1" else si.r2).append(fp)
    return samples


LAYOUT_SCANNERS = (scan_illumina_flat, scan_outsourced_subdir, scan_outsourced_flat)


class ScanResult(tuple):
    """(valid, invalid) 二元组 + 附带 ignored 清单（DEC-31）。
    既有 `valid, invalid = scan_batch(...)` 解包不受影响；ignored 经
    `.ignored` 属性访问（Undetermined 等非样本条目，不算 invalid、不告警）。"""

    def __new__(cls, valid, invalid, ignored):
        self = super().__new__(cls, (valid, invalid))
        self.ignored = ignored
        return self


def scan_batch(batch_dir):
    """扫描批次目录 → (有效样本 dict, 无效样本 dict{sm: reason})，
    附 `.ignored` 名单（DEC-31：IGNORED_SAMPLES 命中的非样本条目，如 Illumina
    下机自带的 Undetermined——剔除后不进分析，独立清单供 run_summary 追溯）。
    布局按 LAYOUT_SCANNERS 依序识别：先认者优先，同一样本名被后到的
    布局再认出 → 后者记冲突无效；不匹配任何布局的文件不构成样本（忽略）。
    输入校验：R1/R2 文件数不一致、单端、0 字节 → 标记无效并跳过。"""
    samples, invalid = {}, {}
    flat, subdirs = [], []
    for name in sorted(os.listdir(batch_dir)):
        p = os.path.join(batch_dir, name)
        if os.path.isfile(p):
            flat.append(p)
        elif os.path.isdir(p):
            subdirs.append(p)

    for scan in LAYOUT_SCANNERS:
        for sm, si in scan(flat, subdirs).items():
            if sm in samples:
                invalid[sm] = f"样本名与{LAYOUT_LABELS.get(samples[sm].layout, samples[sm].layout)}布局冲突"
            elif sm not in invalid:
                samples[sm] = si

    ignored = sorted(sm for sm in samples if sm.lower() in IGNORED_SAMPLES)
    for sm in ignored:
        samples.pop(sm)

    # 校验
    for sm, si in sorted(samples.items()):
        if si.layout == "illumina":
            si.r1.sort(), si.r2.sort()
            si.lanes.sort()
        if len(si.r1) != len(si.r2):
            si.invalid_reason = f"R1({len(si.r1)})/R2({len(si.r2)}) 文件数不一致"
        elif not si.r1:
            si.invalid_reason = "无 fastq 输入"
        else:
            zero = [os.path.basename(f) for f in si.r1 + si.r2
                    if not nonempty(f)]
            if zero:
                si.invalid_reason = f"0 字节文件: {zero}"
        if not si.valid:
            invalid[sm] = si.invalid_reason

    valid = {sm: si for sm, si in samples.items() if si.valid}
    return ScanResult(valid, invalid, ignored)


# ── md5 校验（并行；有清单且失败 → 调用方 P0 中断；无清单跳过，DEC-36）─────
_md5_lock = threading.Lock()


def _md5_one(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def find_md5_manifest(batch_dir):
    """定位批次 md5 清单（DEC-36）：文件名 md5 开头、txt 结尾（均不区分大小写）、
    体积 < config.MD5_MANIFEST_MAX_BYTES（TH-37=500KB，防同名大文本数据文件误认）。
    多候选取字典序首个（其余由调用方 WARN 点名忽略）。→ (清单路径|None, 候选名列表)"""
    try:
        names = sorted(os.listdir(batch_dir))
    except OSError:
        return None, []
    cands = []
    for n in names:
        p = os.path.join(batch_dir, n)
        low = n.lower()
        if not (low.startswith("md5") and low.endswith("txt") and os.path.isfile(p)):
            continue
        try:
            if os.path.getsize(p) >= config.MD5_MANIFEST_MAX_BYTES:
                continue
        except OSError:
            continue
        cands.append(n)
    return (os.path.join(batch_dir, cands[0]) if cands else None), cands


def verify_md5(batch_dir, logger, workers=4):
    """定位 md5 清单（DEC-36：md5*开头/txt 结尾/<TH-37）并校验（并行）。
    返回 (失败样本集, 校验文件数, 状态)：状态 = OK / FAIL / SKIPPED——
    无清单 → SKIPPED 跳过校验（不视为错误）；有清单且失败 → FAIL（调用方 P0 阻断）"""
    md5_file, cands = find_md5_manifest(batch_dir)
    if not md5_file:
        logger.info("无 md5 清单（识别口径：md5 开头/txt 结尾/<500KB），跳过校验")
        return set(), 0, "SKIPPED"
    if len(cands) > 1:
        logger.warn(f"多个 md5 清单候选，取 {cands[0]}（忽略: {', '.join(cands[1:])}）")
    entries = []
    with open(md5_file, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split(None, 1)
            if len(parts) == 2:
                entries.append((parts[0], os.path.normpath(
                    os.path.join(batch_dir, parts[1].strip().lstrip("*")))))
    logger.info(f"{os.path.basename(md5_file)} 存在，校验 {len(entries)} 个文件（{workers} 并行）")
    failed_samples, n_ok, n_bad = set(), 0, 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_md5_one, p): (md5, p) for md5, p in entries}
        for fut in as_completed(futs):
            md5, p = futs[fut]
            try:
                actual = fut.result()
            except OSError as e:
                actual, e_ = None, e
            else:
                e_ = None
            with _md5_lock:
                if actual == md5:
                    n_ok += 1
                else:
                    n_bad += 1
                    # 样本名按布局推导（RUN-33 修复：平铺布局曾误取批次目录名，
                    # 致 md5_failed 与 valid 无交集——失败样本从未被剔除；
                    # RUN-36 补外送平铺式）
                    if os.path.normpath(os.path.dirname(p)) == os.path.normpath(batch_dir):
                        m = ILLUMINA_RE.match(os.path.basename(p)) \
                            or FLAT_OUTSOURCED_RE.match(os.path.basename(p))
                        failed_samples.add(m.group("sample") if m else os.path.basename(p))
                    else:
                        failed_samples.add(os.path.basename(os.path.dirname(p)))
                    logger.error(f"md5 校验失败: {p} 期望 {md5} 实际 {actual} ({e_})")
    state = "FAIL" if failed_samples else "OK"
    logger.result(f"md5 校验完成[{state}]: OK={n_ok} FAIL={n_bad}")
    return failed_samples, len(entries), state


# ── Lane 合并 + samples.tsv ─────────────────────────────────────────────
def merge_sample(si, merged_dir, runner, logger):
    """多 Lane 用 cat 顺序合并为 fastq_merged/<样本>_R1/_R2.fastq.gz"""
    out_r1 = os.path.join(merged_dir, f"{si.sample}_R1.fastq.gz")
    out_r2 = os.path.join(merged_dir, f"{si.sample}_R2.fastq.gz")
    rc = runner.run(
        f"cat {' '.join(repr(f) for f in si.r1)} > {repr(out_r1)} && "
        f"cat {' '.join(repr(f) for f in si.r2)} > {repr(out_r2)}",
        logger=logger, outputs=[out_r1, out_r2])
    if runner.dry_run:
        return True
    return rc == 0 and nonempty(out_r1) and nonempty(out_r2)


def merged_paths(merged_dir, sm):
    return (os.path.join(merged_dir, f"{sm}_R1.fastq.gz"),
            os.path.join(merged_dir, f"{sm}_R2.fastq.gz"))


def merge_all(samples, merged_dir, runner, log_for, workers):
    """并行合并；返回 (成功 dict[sm→(r1,r2)], 失败 dict[sm→reason])"""
    if not runner.dry_run:   # dry-run 零落盘：不建目录（命令仅打印）
        os.makedirs(merged_dir, exist_ok=True)   # cat 输出目录必须先存在
    ok, failed = {}, {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(merge_sample, si, merged_dir, runner,
                            log_for(sm)): sm for sm, si in samples.items()}
        for fut in as_completed(futs):
            sm = futs[fut]
            try:
                ok[sm] = merged_paths(merged_dir, sm) if fut.result() else None
            except Exception as e:   # noqa: BLE001
                ok[sm] = None
                logger = log_for(sm)
                if logger:
                    logger.error(f"合并异常: {e}")
    for sm, v in list(ok.items()):
        if v is None:
            failed[sm] = "Lane 合并失败"
            del ok[sm]
    return ok, failed


def write_samples_tsv(path, merged, samples):
    """三列主结构（样本名/R1/R2）+ 追溯列（S 号/Lane）"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("sample\tR1\tR2\tnote\n")
        for sm in sorted(merged):
            r1, r2 = merged[sm]
            f.write(f"{sm}\t{r1}\t{r2}\t{samples[sm].note()}\n")
