#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""双运行产物一致性比对（RUN-50 验证工具，不入测试集）：
矩阵 TSV 逐字节 md5；VCF 只比记录体（bcftools/gzip 解压后去 # 头行）——
VCF 头含 GATKCommandLine/bcftools 命令行（路径/日期/-Xmx）属元数据，
跨运行/跨目录必不同，体一致即分析结果一致。
用法：python3 compare_runs.py <基线结果目录> <优化结果目录>"""

import gzip
import hashlib
import os
import sys


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_body(path):
    """VCF(.gz) 记录体 md5：解压后剔除 # 开头头行"""
    h = hashlib.md5()
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line)
    return h.hexdigest()


def main(base, opt):
    diff = []
    # ① 三矩阵：逐字节
    for rel in ("matrix/genotype_matrix.tsv",
                "matrix/genotype_detail_PASS.tsv",
                "matrix/genotype_matrix.adjudicated.tsv"):
        b, o = os.path.join(base, rel), os.path.join(opt, rel)
        if not (os.path.isfile(b) and os.path.isfile(o)):
            diff.append(f"[缺失] {rel}")
            continue
        hb, ho = md5_file(b), md5_file(o)
        if hb == ho:
            print(f"[OK] {rel}: {hb}")
        else:
            diff.append(f"[差异] {rel}: {hb} vs {ho}")
    # ② cohort VCF 体
    for rel in ("cohort/cohort.PASS.vcf.gz", "cohort/cohort.hardfiltered.vcf.gz"):
        b, o = os.path.join(base, rel), os.path.join(opt, rel)
        if not (os.path.isfile(b) and os.path.isfile(o)):
            diff.append(f"[缺失] {rel}")
            continue
        hb, ho = md5_body(b), md5_body(o)
        if hb == ho:
            print(f"[OK] {rel} 体: {hb}")
        else:
            diff.append(f"[差异] {rel} 体: {hb} vs {ho}")
    # ③ 每样本裁决 VCF 体（全样本点名）
    ps = os.path.join(base, "per_sample_vcf")
    for name in sorted(os.listdir(ps)):
        if not name.endswith(".PASS.adjudicated.vcf.gz"):
            continue
        b = os.path.join(ps, name)
        o = os.path.join(opt, "per_sample_vcf", name)
        if not os.path.isfile(o):
            diff.append(f"[缺失] per_sample_vcf/{name}")
            continue
        hb, ho = md5_body(b), md5_body(o)
        if hb == ho:
            print(f"[OK] {name} 体: {hb}")
        else:
            diff.append(f"[差异] {name} 体: {hb} vs {ho}")
    print()
    if diff:
        print("=== 不一致清单 ===")
        for d in diff:
            print(d)
        return 1
    print("=== 全部一致（矩阵逐字节 + 全部 VCF 记录体）===")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
