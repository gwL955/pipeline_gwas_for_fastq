#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分级告警与步骤里程碑通知（P0/P1/P2/OK，v2.9.0 三级体系 DEC-21）。

P0：阻断级——严重影响分析的错误（输入损坏/环境不可用/注入面），直接中断批次
P1：严重——执行失败（样本级隔离）/严重数据异常（NTC 污染、mapped<90 QC 口径），
    报错但不中断批次
P2：提示——质量阈值越界/口径存疑（保留率、Q30、dup、覆盖度、Ti/Tv、对账等），
    只报错记录，不影响运行
OK：全部正常，仅常规里程碑播报

样本点名规则（v2.20.0/DEC-30）：逐指标**越界样本全点名**，每样本一行（名字序）；
此前实现每指标只点名最差一个样本（RUN-43 复盘：三样本 on-target 0.58/0.59/0.59
全部越界、钉钉只报 0.58 一个，代表性误读为"只有一个样本坏"）。例外：fastp 保留率
80-95% 提示带为 OK 级，聚合一行点名（防提示刷屏）；批级指标（Ti/Tv、call rate、
深度 CV）本就单值无点名问题。逐样本完整数值以 run_summary.json 为准。

消息模板（钉钉 markdown 官方子集：标题/引用/加粗/列表）：
    [GWAS][P1] 20260720批次 · Step 2 比对完成
    样本: 4/4 成功 | Lane 合并 16/16
    指标: mapped 98.7% | proper pair 94.2%
    异常: L20260615001 mapped 91.3%（阈值 95%）← 需确认
    产物: bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22
    日志: tail -f logs/sample_L20260615001.log
"""

import glob
import os
from datetime import datetime

import config

LEVEL_ORDER = {"P0": 0, "P1": 1, "P2": 2, "OK": 3}


def worst_level(anomalies):
    """anomalies: [(level, msg), ...] → 'P0'/'P1'/'OK'"""
    lv = "OK"
    for level, _ in anomalies or []:
        if LEVEL_ORDER.get(level, 2) < LEVEL_ORDER[lv]:
            lv = level
    return lv


def _fmt_pct(v):
    return f"{v:g}%" if isinstance(v, (int, float)) else "?"


def _avg(d):
    vals = [v for v in (d or {}).values() if isinstance(v, (int, float))]
    return round(sum(vals) / len(vals), 1) if vals else None


def _min_item(d):
    """→ (样本, 最小值)。仅用于指标播报（如 Step3 ELS 最小值）；
    告警点名不走本函数（全点名见 _violating，DEC-30）"""
    vals = [(k, v) for k, v in (d or {}).items() if isinstance(v, (int, float))]
    return min(vals, key=lambda x: x[1]) if vals else (None, None)


def _violating(d, threshold, below=True):
    """→ 越界样本 [(样本, 值), ...]（名字序，全点名 DEC-30）；
    below=True 取 <阈值（越低越坏），False 取 >阈值（越高越坏）"""
    vals = sorted((k, v) for k, v in (d or {}).items()
                  if isinstance(v, (int, float)))
    return [(k, v) for k, v in vals
            if (v < threshold if below else v > threshold)]


# ── 各步阈值检查（返回 anomalies 列表） ────────────────────────────────
def check_fastp(retention, q30):
    """越界样本全点名（DEC-30）：保留率/Q30 逐样本一行；
    80-95% 提示带为 OK 级，聚合一行点名（不逐行刷屏）"""
    out = []
    for sm, v in _violating(retention, config.FASTP_RETENTION_P1):
        out.append(("P2", f"{sm} fastp 保留率 {v}%（阈值 {config.FASTP_RETENTION_P1}%）"))
    band = sorted((k, v) for k, v in (retention or {}).items()
                  if isinstance(v, (int, float))
                  and config.FASTP_RETENTION_P1 <= v < config.FASTP_RETENTION_WARN)
    if band:
        named = "、".join(f"{sm} {v}%" for sm, v in band)
        out.append(("OK", f"fastp 保留率 <{config.FASTP_RETENTION_WARN}% 提示: {named}"))
    for sm, v in _violating(q30, config.FASTP_Q30_P1):
        out.append(("P2", f"{sm} Q30 {v}%（阈值 {config.FASTP_Q30_P1}%）"))
    return out


def check_flagstat(mapped_pct, pp_pct):
    out = []
    for sm, v in _violating(mapped_pct, config.MAPPED_NOTIFY_P1):
        out.append(("P2", f"{sm} mapped {v}%（阈值 {config.MAPPED_NOTIFY_P1}%）"))
    for sm, v in _violating(pp_pct, config.PROPER_PAIR_P1):
        out.append(("P2", f"{sm} properly paired {v}%（阈值 {config.PROPER_PAIR_P1}%）"))
    return out


def check_dup(dup_pct):
    return [("P2", f"{sm} 重复率 {v}%（阈值 {config.DUP_P1}%，建库复杂度告急）")
            for sm, v in _violating(dup_pct, config.DUP_P1, below=False)]


def check_capture(mean_depth, pct20x, pct_selected):
    """捕获效率口径 v2.19.0/DEC-29：PCT_SELECTED_BASES（on+near bait 占比对碱基比）
    为告警指标；on-target 不再告警（1bp SNP panel 下为几何产物，见 run_summary 信息指标）。
    越界样本全点名（DEC-30）"""
    out = []
    for sm, v in _violating(mean_depth, config.MEAN_DEPTH_P1):
        out.append(("P2", f"{sm} mean depth {v}×（阈值 {config.MEAN_DEPTH_P1}×）"))
    for sm, v in _violating(pct20x, config.PCT_20X_P1):
        out.append(("P2", f"{sm} ≥20x 靶比例 {v}%（阈值 {config.PCT_20X_P1}%）"))
    for sm, v in _violating(pct_selected, config.PCT_SELECTED_P1):
        out.append(("P2", f"{sm} 捕获效率 PCT_SELECTED {v}%（阈值 {config.PCT_SELECTED_P1}%）"))
    return out


def check_variantqc(titv_pass, call_rate):
    out = []
    if titv_pass is not None and titv_pass < config.TITV_P1:
        out.append(("P2", f"Ti/Tv {titv_pass}（阈值 {config.TITV_P1}，结论口径存疑）"))
    if call_rate is not None and call_rate < config.CALL_RATE_P1:
        out.append(("P2", f"call rate {call_rate}%（阈值 {config.CALL_RATE_P1}%，非 ./. 基因型占比）"))
    return out


def check_ntc(ntc_depth):
    if ntc_depth is not None and ntc_depth > config.NTC_DEPTH_P0:
        return [("P1", f"NTC 靶区深度 {ntc_depth}×（阈值 {config.NTC_DEPTH_P0}×）——阴性对照出现真实覆盖，疑似污染")]
    return []


def check_reads_low(reads):
    """reads: {sm: after_reads}——绝对量过低 → P1（上样不足，报错不中断，RUN-34）"""
    out = []
    for sm, n in sorted((reads or {}).items()):
        if isinstance(n, (int, float)) and n < config.READS_MIN:
            out.append(("P1", f"{sm} reads {int(n)} < {config.READS_MIN}"
                               f"（上样不足，结果可信度存疑）"))
    return out


def check_ntc_reads(ntc_reads, median_reads):
    """NTC reads 占批次中位样本 reads 比例异常 → P2（污染维度之二，与深度互补）"""
    if not isinstance(ntc_reads, (int, float)) or not isinstance(median_reads, (int, float)) \
            or median_reads <= 0:
        return []
    pct = ntc_reads / median_reads * 100
    if pct > config.NTC_READS_PCT_P2:
        return [("P2", f"NTC reads 占批次中位样本 {pct:.2f}%"
                       f"（阈值 {config.NTC_READS_PCT_P2}%）——阴性对照出现可观数据量，疑似污染")]
    return []


def check_depth_cv(mean_depth):
    """批次内样本间 mean depth 变异系数 CV 过大 → P2（疑似混入异常样本）"""
    vals = [v for v in (mean_depth or {}).values() if isinstance(v, (int, float))]
    if len(vals) < 2:
        return []
    mean = sum(vals) / len(vals)
    if mean <= 0:
        return []
    cv = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5 / mean
    if cv > config.DEPTH_CV_P2:
        return [("P2", f"批次内 mean depth 离散度 CV={cv:.2f}（阈值 {config.DEPTH_CV_P2}，"
                       f"深度 {min(vals):.1f}-{max(vals):.1f}×）——疑似混入异常样本")]
    return []


def check_recal_low(obs):
    """recal: {sm: M 事件观测数}——known-sites 覆盖崩坏致校准不可信 → P2"""
    out = []
    for sm, n in sorted((obs or {}).items()):
        if isinstance(n, (int, float)) and n < config.RECAL_OBS_MIN_P2:
            out.append(("P2", f"{sm} BQSR recal 观测数 {int(n)} < {int(config.RECAL_OBS_MIN_P2)}"
                              f"（known-sites 覆盖异常，校准不可信）"))
    return out


# ── 产物摘要 ────────────────────────────────────────────────────────────
def artifact_summary(pattern, basedir=None):
    """glob 产物 → '×N  mtime …'（取最新 mtime）"""
    hits = glob.glob(pattern)
    if not hits:
        return "×0"
    latest = max(os.path.getmtime(p) for p in hits)
    return (f"×{len(hits)}  mtime "
            f"{datetime.fromtimestamp(latest).strftime('%Y-%m-%d %H:%M')}")


# ── 里程碑消息构造 ──────────────────────────────────────────────────────
def step_milestone(batch, step_no, step_name, *, samples="", metrics="",
                   anomalies=(), artifacts="", log_hint="", extra=()):
    """→ (title, text)。title 带 [GWAS][P级别]，正文含小写 gwas 兜底由 dingtalk 补。"""
    level = worst_level(anomalies)
    title = f"[GWAS][{level}] {batch}批次 · " \
            f"{f'Step {step_no} ' if step_no is not None else ''}{step_name}完成"
    text = f"#### {title}"
    if samples:
        text += f"\n\n样本: {samples}"
    if metrics:
        text += f"\n\n指标: {metrics}"
    for line in extra:
        if line:
            text += f"\n\n{line}"
    if anomalies:
        for lv, msg in anomalies:
            tag = {"P0": " ← 阻断级（中断分析）", "P1": " ← 需确认"}.get(lv, "")
            text += f"\n\n异常: [{lv}] {msg}{tag}"
    else:
        text += "\n\n异常: 无"
    if artifacts:
        text += f"\n\n产物: {artifacts}"
    if log_hint:
        text += f"\n\n日志: {log_hint}"
    return title, text
