#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bwa-mem2：比对 + RG 构造（管道流式接 samtools sort，不落中间 sam）。
对应笔记《2-比对》代码块 3/16。"""

import config


def rg_line(sample):
    """单 RG（多 Lane 合并后 ID=SM=样本名即可）"""
    return f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"


def build_align_pipe(runner, sample, r1_host, r2_host, out_bam_host,
                     bwa_threads, sort_threads, sort_mem):
    """bwa-mem2 mem -t T -K 100000000 -Y -R ... | samtools sort -@ T -m M -O BAM -o ... -
    bwa 与 sort 同管道共存，峰值内存 = 索引 ~17G + 两者缓冲之和"""
    bwa = runner.tool(
        "bwa",
        f"bwa-mem2 mem -t {bwa_threads} -K {config.BWA_K} -Y "
        f'-R "{rg_line(sample)}" '
        f"{runner.cpath(config.GENOME_FA)} "
        f"{runner.cpath(r1_host)} {runner.cpath(r2_host)}")
    sort = runner.tool(
        "samtools",
        f"samtools sort -@ {sort_threads} -m {sort_mem} -O BAM "
        f"-o {runner.cpath(out_bam_host)} -")
    return f"{bwa} | {sort}"
