#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""samtools：sort 管道由 bwa_mem2 构造；此处提供 index / flagstat / stats
运行与解析（三套统计）。对应笔记《2-比对》代码块 5/7/14/16。"""

import re
import config

# 去重前 flagstat 的 duplicates 行无意义，不采集（提示词硬性要求）
FLAG_IGNORE_KEYS = ("duplicates", "primary duplicates")


def run_index(runner, bam_host, logger):
    rc = runner.run(runner.tool("samtools",
                                f"samtools index {runner.cpath(bam_host)}"),
                    logger=logger, outputs=[str(bam_host) + ".bai"])
    return rc == 0


def run_flagstat(runner, bam_host, out_host, logger):
    rc = runner.run(
        f"{runner.tool('samtools', f'samtools flagstat {runner.cpath(bam_host)}')} > {out_host}",
        logger=logger, outputs=[out_host], timeout=config.STATS_TIMEOUT_S)
    return rc == 0


def run_stats(runner, bam_host, out_host, logger):
    rc = runner.run(
        f"{runner.tool('samtools', f'samtools stats {runner.cpath(bam_host)}')} > {out_host}",
        logger=logger, outputs=[out_host], timeout=config.STATS_TIMEOUT_S)
    return rc == 0


def parse_flagstat(path):
    """→ dict(total, mapped, mapped_pct, paired, properly_paired, pp_pct,
    singletons, sgl_pct, diff_chr, diff_chr_mapq5)"""
    d = {}
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                m = re.match(r"^(\d+) \+ \d+ ", line)
                if not m:
                    continue
                n = int(m.group(1))
                rest = line[m.end():]
                pm = re.search(r"\(([\d.]+)%", rest)
                if pm:                       # 带 " (xx.xx%" 的行：键为括号前部分
                    key, pct = rest[:pm.start()].strip(), float(pm.group(1))
                elif "(mapQ>=5)" in rest:    # 该行无百分比，须与普通跨染色体行区分
                    d["diff_chr_mapq5"] = n
                    continue
                else:                        # in total 行无百分比，键为 " (" 前部分
                    key, pct = rest.split(" (")[0].strip(), None
                if key == "in total":
                    d["total"] = n
                elif key == "mapped" and pct is not None:
                    d["mapped"], d["mapped_pct"] = n, pct
                elif key == "paired in sequencing":
                    d["paired"] = n
                elif key == "properly paired" and pct is not None:
                    d["properly_paired"], d["pp_pct"] = n, pct
                elif key == "singletons" and pct is not None:
                    d["singletons"], d["sgl_pct"] = n, pct
                elif key == "with mate mapped to a different chr":
                    d["diff_chr"] = n
                elif key == "with mate mapped to a different chr (mapQ>=5)":
                    d["diff_chr_mapq5"] = n
    except OSError:
        pass
    if "total" in d and d.get("mapped") is not None:
        d["mapped_calc_pct"] = round(d["mapped"] * 100.0 / d["total"], 2) if d["total"] else 0.0
    return d


def parse_stats(path):
    """samtools stats SN 行 → dict（raw total / mapped / properly paired / insert size / qual）"""
    d = {}
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if not line.startswith("SN\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                key = parts[1].rstrip(":").strip()
                try:
                    d[key] = float(parts[2])
                except ValueError:
                    d[key] = parts[2]
    except OSError:
        pass
    return d


def flagstat_identical(path_a, path_b):
    """BQSR 前后 flagstat 必须逐行一致（断言检查）"""
    try:
        with open(path_a, encoding="utf-8") as fa, open(path_b, encoding="utf-8") as fb:
            return fa.read() == fb.read()
    except OSError:
        return False


def qc_judgement(flagstat, logger):
    """质检口径：mapped>90%；singletons 应很低；properly paired ≈97-98%
    （偏低且伴跨染色体配对升高 → 节段重复区正常现象，记录告警不报错）"""
    ok = True
    mp = flagstat.get("mapped_pct")
    if mp is not None and mp < 90.0:
        logger.error(f"mapped {mp}% < 90%，比对率不达标")
        ok = False
    elif mp is not None:
        logger.result(f"mapped {mp}% ≥ 90%")
    sp = flagstat.get("sgl_pct")
    if sp is not None and sp > 1.0:
        logger.warn(f"singletons {sp}% 偏高（应很低）")
    pp = flagstat.get("pp_pct")
    dc = flagstat.get("diff_chr")
    if pp is not None and pp < 96.0:
        msg = f"properly paired {pp}% 偏低（参考 97-98%）"
        if dc and flagstat.get("paired") and dc / flagstat["paired"] > 0.03:
            msg += "；伴跨染色体配对升高 → 节段重复区正常现象（MAPQ=0 为主），告警不报错"
        logger.warn(msg)
    return ok
