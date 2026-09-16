#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分级告警与步骤里程碑通知（P0/P1/OK）。

P0：阻断级/污染级——执行失败、合并失败、依赖缺失、磁盘不足、NTC 污染
P1：需确认——各步质量阈值越界、对账数量/新鲜度、命名规范
OK：全部正常，仅常规里程碑播报

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

LEVEL_ORDER = {"P0": 0, "P1": 1, "OK": 2}


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
    """→ (样本, 最小值)"""
    vals = [(k, v) for k, v in (d or {}).items() if isinstance(v, (int, float))]
    return min(vals, key=lambda x: x[1]) if vals else (None, None)


def _max_item(d):
    vals = [(k, v) for k, v in (d or {}).items() if isinstance(v, (int, float))]
    return max(vals, key=lambda x: x[1]) if vals else (None, None)


# ── 各步阈值检查（返回 anomalies 列表） ────────────────────────────────
def check_fastp(retention, q30):
    out = []
    sm, v = _min_item(retention)
    if v is not None and v < config.FASTP_RETENTION_P1:
        out.append(("P1", f"{sm} fastp 保留率 {v}%（阈值 {config.FASTP_RETENTION_P1}%）"))
    elif v is not None and v < config.FASTP_RETENTION_WARN:
        out.append(("OK", f"{sm} fastp 保留率 {v}%（<{config.FASTP_RETENTION_WARN}% 提示）"))
    sm, v = _min_item(q30)
    if v is not None and v < config.FASTP_Q30_P1:
        out.append(("P1", f"{sm} Q30 {v}%（阈值 {config.FASTP_Q30_P1}%）"))
    return out


def check_flagstat(mapped_pct, pp_pct):
    out = []
    sm, v = _min_item(mapped_pct)
    if v is not None and v < config.MAPPED_NOTIFY_P1:
        out.append(("P1", f"{sm} mapped {v}%（阈值 {config.MAPPED_NOTIFY_P1}%）"))
    sm, v = _min_item(pp_pct)
    if v is not None and v < config.PROPER_PAIR_P1:
        out.append(("P1", f"{sm} properly paired {v}%（阈值 {config.PROPER_PAIR_P1}%）"))
    return out


def check_dup(dup_pct):
    sm, v = _max_item(dup_pct)
    if v is not None and v > config.DUP_P1:
        return [("P1", f"{sm} 重复率 {v}%（阈值 {config.DUP_P1}%，建库复杂度告急）")]
    return []


def check_capture(mean_depth, pct20x, on_target):
    out = []
    sm, v = _min_item(mean_depth)
    if v is not None and v < config.MEAN_DEPTH_P1:
        out.append(("P1", f"{sm} mean depth {v}×（阈值 {config.MEAN_DEPTH_P1}×）"))
    sm, v = _min_item(pct20x)
    if v is not None and v < config.PCT_20X_P1:
        out.append(("P1", f"{sm} ≥20x 靶比例 {v}%（阈值 {config.PCT_20X_P1}%）"))
    sm, v = _min_item(on_target)
    if v is not None and v < config.ON_TARGET_P1:
        out.append(("P1", f"{sm} on-target {v}%（阈值 {config.ON_TARGET_P1}%）"))
    return out


def check_variantqc(titv_pass, call_rate):
    out = []
    if titv_pass is not None and titv_pass < config.TITV_P1:
        out.append(("P1", f"Ti/Tv {titv_pass}（阈值 {config.TITV_P1}，结论口径存疑）"))
    if call_rate is not None and call_rate < config.CALL_RATE_P1:
        out.append(("P1", f"call rate {call_rate}%（阈值 {config.CALL_RATE_P1}%，非 ./. 基因型占比）"))
    return out


def check_ntc(ntc_depth):
    if ntc_depth is not None and ntc_depth > config.NTC_DEPTH_P0:
        return [("P0", f"NTC 靶区深度 {ntc_depth}×（阈值 {config.NTC_DEPTH_P0}×）——阴性对照出现真实覆盖，疑似污染")]
    return []


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
            tag = " ← 需确认" if lv == "P1" else (" ← 阻断/污染级" if lv == "P0" else "")
            text += f"\n\n异常: [{lv}] {msg}{tag}"
    else:
        text += "\n\n异常: 无"
    if artifacts:
        text += f"\n\n产物: {artifacts}"
    if log_hint:
        text += f"\n\n日志: {log_hint}"
    return title, text
