#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mosdepth：运行（markdup 与 BQSR 两个 BAM 各一次）+ summary/regions 解析，
是矩阵 ./. 裁决的数据源。对应笔记《5-测序质量》代码块 0/1。"""

import config
from modules.bcftools import gzopen_text


def run_mosdepth(runner, bam, prefix_host, threads, logger):
    """mosdepth --by targets.sorted.bed --fast-mode -t N <prefix> <bam>"""
    summary = f"{prefix_host}.mosdepth.summary.txt"
    rc = runner.run(
        runner.tool("mosdepth",
                    f"mosdepth --by {runner.cpath(config.TARGETS_SORTED_BED)} "
                    f"--fast-mode -t {threads} "
                    f"{runner.cpath(prefix_host)} {runner.cpath(bam)}"),
        logger=logger, outputs=[summary, f"{prefix_host}.regions.bed.gz"])
    return rc == 0


def parse_summary(prefix):
    """summary（列：chrom length bases mean min max）→ {'mean': 靶区 total_region 均值, ...}。
    --by 模式下 total 行是全基因组口径（mean 被摊薄），靶区深度取 total_region 行。"""
    rows = {}
    try:
        with open(f"{prefix}.mosdepth.summary.txt", encoding="utf-8") as f:
            next(f)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    try:
                        rows[parts[0]] = (int(parts[1]), int(parts[2]),
                                          float(parts[3]))
                    except ValueError:
                        continue
    except (OSError, StopIteration):
        pass
    region = rows.get("total_region") or rows.get("total")
    return {"chroms": {k: v[2] for k, v in rows.items()},
            "mean": region[2] if region else None}


def load_regions(prefix):
    """regions.bed.gz → list[(chrom, start, end, mean)]（按基因组排序缓存）"""
    regions = []
    try:
        with gzopen_text(f"{prefix}.regions.bed.gz") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    try:
                        regions.append((parts[0], int(parts[1]), int(parts[2]),
                                        float(parts[3])))
                    except ValueError:
                        continue
    except OSError:
        pass
    regions.sort(key=lambda r: (r[0], r[1]))
    return regions


def region_depth(regions, chrom, pos):
    """定位 pos 落在哪个 region（线性扫描 + 缓存索引），返回 mean 或 None"""
    lo, hi = 0, len(regions)
    while lo < hi:
        mid = (lo + hi) // 2
        c, s, e, _ = regions[mid]
        if c < chrom or (c == chrom and e <= pos):
            lo = mid + 1
        else:
            hi = mid
    if lo < len(regions):
        c, s, e, m = regions[lo]
        if c == chrom and s <= pos < e:
            return m
    return None
