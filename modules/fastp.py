#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fastp：修剪运行 + json 指标解析（保留率等）。
对应笔记《1-原始数据 QC + 修剪》代码块 10。"""

import json

import config


def run_fastp(runner, sample, r1_host, r2_host, out_r1, out_r2, html, json_path,
              threads, logger):
    """fastp -i R1 -I R2 -o clean_R1 -O clean_R2 -h html -j json
    --length_required 36 --thread N（自动接头检测）"""
    rc = runner.run(
        runner.tool("fastp",
                    f"fastp -i {runner.cpath(r1_host)} -I {runner.cpath(r2_host)} "
                    f"-o {runner.cpath(out_r1)} -O {runner.cpath(out_r2)} "
                    f"-h {runner.cpath(html)} -j {runner.cpath(json_path)} "
                    f"--length_required {config.FASTP_LENGTH_REQUIRED} --thread {threads}"),
        logger=logger, outputs=[out_r1, out_r2, html, json_path])
    return rc == 0


def parse_json(json_path):
    """→ {before/after reads, retention_pct, q30_pct(百分数), q30_pct_r1/r2...}
    键名按 fastp 1.3.6 真实产物校准：after_filtering 只有整体 q30_rate（小数），
    R1/R2 分列在顶层 read{1,2}_after_filtering 且无 q30_rate——用 q30_bases/
    total_bases 自算（RUN-25：曾取不存在的 q30_rate_r1 致 Q30 恒为 None）"""
    d = {}
    try:
        with open(json_path, encoding="utf-8") as f:
            j = json.load(f)
        bf = j.get("summary", {}).get("before_filtering", {})
        af = j.get("summary", {}).get("after_filtering", {})
        d["before_reads"] = bf.get("total_reads")
        d["after_reads"] = af.get("total_reads")
        d["before_bases"] = bf.get("total_bases")
        d["after_bases"] = af.get("total_bases")
        if bf.get("total_reads"):
            d["retention_pct"] = round(af.get("total_reads", 0) * 100.0
                                       / bf["total_reads"], 2)
        q30 = af.get("q30_rate")
        if isinstance(q30, (int, float)):
            d["q30_pct"] = round(q30 * 100, 2)
        for tag, key in (("q30_pct_r1", "read1_after_filtering"),
                         ("q30_pct_r2", "read2_after_filtering")):
            ra = j.get(key) or {}
            if ra.get("total_bases"):
                d[tag] = round(ra.get("q30_bases", 0) * 100.0
                               / ra["total_bases"], 2)
        filt = j.get("filtering_result", {})
        d["low_quality_filtered"] = filt.get("low_quality_reads")
        d["too_short_filtered"] = filt.get("too_short_reads")
        ad = j.get("adapter_cutting", {}) or {}
        d["adapter_trimmed_reads"] = ad.get("adapter_trimmed_reads")
    except (OSError, ValueError):
        pass
    return d


def check_retention(metrics, logger):
    """保留率 <FASTP_RETENTION_WARN（TH-02）的样本告警"""
    pct = metrics.get("retention_pct")
    if pct is None:
        return
    if pct < config.FASTP_RETENTION_WARN:
        logger.warn(f"fastp 保留率 {pct}% < {config.FASTP_RETENTION_WARN}%，请检查原始数据质量")
    else:
        logger.result(f"fastp 保留率 {pct}% ≥ {config.FASTP_RETENTION_WARN}%")
