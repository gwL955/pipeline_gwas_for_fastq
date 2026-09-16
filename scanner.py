#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""输入布局扫描（Illumina 平铺多 Lane / 外送子目录两种布局）+ md5 校验 +
Lane 合并 + samples.tsv 生成。批次是最小分析单元，本模块按单个批次工作。"""

import os
import re
import hashlib
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

from runner import nonempty

ILLUMINA_RE = re.compile(
    r"^(?P<sample>.+)_S(?P<snum>\d+)_L(?P<lane>\d{3})_R(?P<read>[12])_001\.fastq\.gz$")
OUTSOURCED_RE = re.compile(r"_R(?P<read>[12])\.fastq\.gz$")


class SampleInfo:
    def __init__(self, sample, layout):
        self.sample = sample
        self.layout = layout          # "illumina" / "outsourced"
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
        return f"外送;{len(self.r1)}文件"


def scan_batch(batch_dir):
    """扫描批次目录 → (有效样本 dict, 无效样本 dict{sm: reason})。
    输入校验：R1/R2 文件数不一致、单端、0 字节 → 标记无效并跳过。"""
    samples, invalid = {}, {}
    flat, subdirs = [], []
    for name in sorted(os.listdir(batch_dir)):
        p = os.path.join(batch_dir, name)
        if os.path.isfile(p):
            flat.append(p)
        elif os.path.isdir(p):
            subdirs.append(p)

    # 布局 1：Illumina 平铺（BCLConvert 命名，多 Lane）
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

    # 布局 2：外送子目录 <样本名>/<样本名>_R1.fastq.gz
    for dp in subdirs:
        r1s, r2s = [], []
        for name in sorted(os.listdir(dp)):
            m = OUTSOURCED_RE.search(name)
            if m:
                (r1s if m.group("read") == "1" else r2s).append(os.path.join(dp, name))
        if not r1s and not r2s:
            continue
        sm = os.path.basename(dp)
        if sm in samples:
            invalid[sm] = "样本名与平铺布局冲突"
            continue
        si = SampleInfo(sm, "outsourced")
        si.r1, si.r2 = r1s, r2s
        samples[sm] = si

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
    return valid, invalid


# ── md5 校验（并行；失败样本终止分析并进入通知） ────────────────────────
_md5_lock = threading.Lock()


def _md5_one(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_md5(batch_dir, logger, workers=4):
    """批次内存在 md5sum.txt 时校验（并行）。返回 (失败样本集, 校验文件数)"""
    md5_file = os.path.join(batch_dir, "md5sum.txt")
    if not os.path.isfile(md5_file):
        return set(), 0
    entries = []
    with open(md5_file, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split(None, 1)
            if len(parts) == 2:
                entries.append((parts[0], os.path.normpath(
                    os.path.join(batch_dir, parts[1].strip().lstrip("*")))))
    logger.info(f"md5sum.txt 存在，校验 {len(entries)} 个文件（{workers} 并行）")
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
                    failed_samples.add(os.path.basename(
                        os.path.dirname(p)) or os.path.basename(p))
                    logger.error(f"md5 校验失败: {p} 期望 {md5} 实际 {actual} ({e_})")
    logger.result(f"md5 校验完成: OK={n_ok} FAIL={n_bad}")
    return failed_samples, len(entries)


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
