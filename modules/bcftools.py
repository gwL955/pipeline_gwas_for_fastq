#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bcftools：norm / view / concat / query / stats / index 及矩阵导出、
每样本 VCF 拆分与 GT 替换重建。对应笔记《4-变异检测》《5-测序质量》。"""

import gzip
import re
import shlex

import config


def _bc(runner, args):
    return runner.tool("bcftools", f"bcftools {args}")


# ── 基础操作 ──
def index_tbi(runner, vcf, logger):
    rc = runner.run(_bc(runner, f"index -t {runner.cpath(vcf)}"),
                    logger=logger, outputs=[str(vcf) + ".tbi"])
    return rc == 0


def norm_split(runner, in_vcf, out_vcf, logger):
    """bcftools norm -f genome.fa -m -any：拆多等位 + 左对齐 + REF 校验 + 建 tbi 索引
    （笔记 4-9：索引是 SelectVariants/多文件统计的前置）。
    返回 (ok, stats_dict{total/split/realigned/skipped})。
    norm 命令的 outputs 只写 vcf——.tbi 由随后的 index 命令生成，混进本命令的
    产物检查会在两命令之间必然误报"命令成功但产物缺失/为空"（RUN-29）。"""
    rc, out, err = runner.run(
        _bc(runner, f"norm -f {runner.cpath(config.GENOME_FA)} -m -any "
                    f"{runner.cpath(in_vcf)} -Oz -o {runner.cpath(out_vcf)}"),
        logger=logger, outputs=[out_vcf], capture=True)
    stats = _parse_norm_stats((out or "") + (err or ""))
    if rc == 0:
        index_tbi(runner, out_vcf, logger)
    return rc == 0, stats


def _parse_norm_stats(text):
    """解析 'Lines total/split/joined/realigned/.../skipped: v1/v2/...' 统计行（按表头名取值）"""
    m = re.search(r"Lines\s+([\w/]+):\s*([\d/]+)", text or "")
    if not m:
        return {}
    keys = m.group(1).split("/")
    vals = [int(v) for v in m.group(2).split("/") if v.lstrip("-").isdigit()]
    if len(keys) != len(vals):
        return {}
    return dict(zip(keys, vals))


def view_pass(runner, in_vcf, out_vcf, logger):
    rc = runner.run(_bc(runner,
                        f"view -f PASS -Oz -o {runner.cpath(out_vcf)} "
                        f"{runner.cpath(in_vcf)}"),
                    logger=logger, outputs=[out_vcf])
    return rc == 0


def concat(runner, vcfs, out_vcf, logger):
    rc = runner.run(_bc(runner,
                        f"concat -a {' '.join(runner.cpath(v) for v in vcfs)} "
                        f"-Oz -o {runner.cpath(out_vcf)}"),
                    logger=logger, outputs=[out_vcf])
    return rc == 0


def count_records(runner, vcf, logger=None):
    out = runner.out(
        f"{_bc(runner, f'view -H {runner.cpath(vcf)}')} | wc -l", logger=logger)
    try:
        return int(out.strip().split()[-1])
    except (ValueError, IndexError):
        return -1


def restrict_targets(runner, in_vcf, out_vcf, bed_host, pad, logger):
    """靶区提取加 ±100bp padding（防左对齐后 indel 移出区间）；
    生成 pad BED 后 bcftools view -T。返回 (ok, bed_path)"""
    bed_pad = str(bed_host) + f".pad{pad}.bed"
    lines = []
    try:
        with open(bed_host, encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    lines.append(f"{parts[0]}\t{max(0, int(parts[1]) - pad)}\t{int(parts[2]) + pad}")
    except OSError:
        return False, bed_pad
    with open(bed_pad, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))
    rc = runner.run(_bc(runner,
                        f"view -T {runner.cpath(bed_pad)} "
                        f"-Oz -o {runner.cpath(out_vcf)} {runner.cpath(in_vcf)}"),
                    logger=logger, outputs=[out_vcf])
    return rc == 0, bed_pad


# ── 矩阵 / 位点详情导出（双口径） ──
# 注：格式串经 shell 单引号原样传给 bcftools，\\t 由 bcftools 自行解释为制表符
FMT_MATRIX = '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'
FMT_DETAIL = '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT\\t%DP\\t%AD]\\n'
FMT_GT_ONLY = '%CHROM\\t%POS\\t%GT\\n'
FMT_FILTER = '%FILTER\\n'


def export_genotype_matrix(runner, in_vcf, out_tsv, logger):
    """全口径矩阵（hardfiltered × targets）：CHROM POS REF ALT [GT...]"""
    cmd = (_bc(runner, "query -R " + runner.cpath(config.TARGETS_SORTED_BED)
               + " -f " + shlex.quote(FMT_MATRIX) + " " + runner.cpath(in_vcf))
           + f" > {out_tsv}")
    rc = runner.run(cmd, logger=logger, outputs=[out_tsv])
    return rc == 0


def export_detail_pass(runner, in_vcf, out_tsv, logger):
    """PASS 位点 + 质控详情（GT:DP:AD——下游关联分析的正式输入）"""
    cmd = (_bc(runner, "query -f " + shlex.quote(FMT_DETAIL) + " "
               + runner.cpath(in_vcf)) + f" > {out_tsv}")
    rc = runner.run(cmd, logger=logger, outputs=[out_tsv])
    return rc == 0


def query_lines(runner, vcf, fmt, logger=None):
    """bcftools query -f <fmt>，返回行列表（内存中比对用）"""
    out = runner.out(
        _bc(runner, "query -f " + shlex.quote(fmt) + " " + runner.cpath(vcf)),
        logger=logger)
    return [l for l in out.splitlines() if l]


# ── stats 解析 ──
def run_stats(runner, vcf, out_host, logger):
    rc = runner.run(f"{_bc(runner, f'stats {runner.cpath(vcf)}')} > {out_host}",
                    logger=logger, outputs=[out_host])
    return rc == 0


def parse_stats(path):
    """bcftools stats → {records, snps, indels, ts, tv, titv...}
    SN 行格式：SN <fileID> <key>: <value>（键值按名匹配，不依赖列位置）"""
    d = {}
    mapping = {
        "number of records": "records",
        "number of no-ALTs": "noalts",
        "number of SNPs": "snps",
        "number of MNPs": "mnps",
        "number of indels": "indels",
        "number of others": "others",
        "number of multiallelic sites": "multiallelic",
    }
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("SN\t"):
                    parts = line.rstrip("\n").split("\t")
                    for i, p in enumerate(parts):
                        key = p.rstrip(":").strip()
                        if key in mapping and i + 1 < len(parts):
                            try:
                                d[mapping[key]] = int(parts[i + 1])
                            except ValueError:
                                pass
                elif line.startswith("TSTV\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 6:
                        try:
                            d["ts"], d["tv"], d["titv"] = \
                                int(parts[2]), int(parts[3]), float(parts[4])
                        except ValueError:
                            pass
    except OSError:
        pass
    return d


# ── FILTER 列质检 ──
def filter_column_check(runner, vcf, logger):
    """FILTER 列不得出现 '.'（出现 = 过滤漏跑，判错）"""
    out = runner.out(_bc(runner, "query -f " + shlex.quote(FMT_FILTER) + " "
                         + runner.cpath(vcf)))
    dist = {}
    for line in out.splitlines():
        if line:
            dist[line] = dist.get(line, 0) + 1
    if "." in dist:
        logger.error(f"FILTER 列出现 '.'（{dist['.']} 条）——过滤漏跑，判 FAIL")
        return False, dist
    return True, dist


# ── 每样本 VCF：拆分 + GT 替换重建（./. 裁决后） ──
def split_sample(runner, in_vcf, sample, out_vcf, logger):
    rc = runner.run(_bc(runner,
                        f"view -s {sample} --min-ac 0 -Oz "
                        f"-o {runner.cpath(out_vcf)} {runner.cpath(in_vcf)}"),
                    logger=logger, outputs=[out_vcf])
    return rc == 0


def adjudicate_matrix(matrix_tsv, regions_reader, samples, dp_min):
    """矩阵 ./. 裁决：mosdepth regions 深度 DP≥20 改判 0/0，不足保留 ./.（DP_MIN 可配）。
    regions_reader(sm) → iterable[(chrom, start, end, mean)]。纯 Python。"""
    n_total = n_filled = n_kept = 0
    n_cells = 0
    out_lines = []
    with open(matrix_tsv, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 + len(samples):
                out_lines.append(line.rstrip("\n"))
                continue
            n_cells += len(samples)
            chrom, pos = p[0], int(p[1])
            new_gts = []
            for sm, gt in zip(samples, p[4:]):
                if gt == "./.":
                    n_total += 1
                    mean = regions_reader(sm, chrom, pos)
                    if mean is not None and mean >= dp_min:
                        new_gts.append("0/0")
                        n_filled += 1
                    else:
                        new_gts.append("./.")
                        n_kept += 1
                else:
                    new_gts.append(gt)
            out_lines.append("\t".join(p[:4] + new_gts))
    return out_lines, {"dotdot_total": n_total, "filled_00": n_filled,
                       "kept_dotdot": n_kept, "cells_total": n_cells}


def rebuild_sample_vcf(in_vcf_plain, gt_tsv, out_plain, sample):
    """方案一（笔记 5-3）：裁决后矩阵列 → 该样本 VCF 只替换 GT（其余字段原样保留）。纯 Python。"""
    new_gt = {}
    with open(gt_tsv, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                new_gt[(p[0], int(p[1]))] = p[4]
            elif len(p) == 3:
                new_gt[(p[0], int(p[1]))] = p[2]
    n_total = n_changed = 0
    with open(in_vcf_plain, encoding="utf-8") as f, \
            open(out_plain, "w", encoding="utf-8") as fo:
        for line in f:
            if line.startswith("#"):
                fo.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            n_total += 1
            if len(p) >= 10:
                chrom, pos = p[0], int(p[1])
                fmt = p[8].split(":")
                cols = p[9].split(":")
                if fmt[0] == "GT" and (chrom, pos) in new_gt:
                    cols[0] = new_gt[(chrom, pos)]
                    p[9] = ":".join(cols)
                    n_changed += 1
            fo.write("\t".join(p) + "\n")
    return {"sample": sample, "records": n_total, "gt_changed": n_changed}


def gzopen_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") \
        else open(path, encoding="utf-8")
