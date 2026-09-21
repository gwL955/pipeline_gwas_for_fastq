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
80-90% 提示带（TH-02~03，v2.21.0 由 80-95 下调）为 OK 级，聚合一行点名（防提示
刷屏）；批级指标（Ti/Tv、call rate、深度 CV）本就单值无点名问题。逐样本完整数值
以 run_summary.json 为准。

对照样本豁免（v2.21.0/DEC-31）：excluded（--exclude-samples，默认 NTC）命中的
对照样本豁免样本级阈值——reads 低属阴性对照正常态，check_reads_low 降级为 OK 级
提示行播报实际数值；其余样本级指标（保留率/Q30/mapped/dup/捕获/recal/CV）对对照
无统计意义，直接跳过。污染监控不在此列，仍由 check_ntc（靶区深度）与
check_ntc_reads（占批次中位比例）专属口径负责。

消息模板（钉钉 markdown 官方子集：标题/引用/加粗/列表）：
    [GWAS][P1] 20260720批次 · Step 2 比对完成（共 6 步）
    样本: 4/4 成功 | Lane 合并 16/16
    指标: mapped 98.7% | proper pair 94.2%
    异常: L20260615001 mapped 91.3%（阈值 95%）← 需确认
    产物: bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22
    日志: tail -f logs/sample_L20260615001.log

指标播报对照豁免（v2.22.0/DEC-32）：DEC-31 只豁免了"告警点名"，"指标:"行的
批均值与 Step3 ELS 统计此前仍含对照——NTC reads 近 0 使 dup/ELS/深度/mapped
均值失真、ELS 最低值被 NTC 占据（误读为文库复杂度不足）。v2.22.0 起里程碑与
全流程汇总的样本级指标均值、ELS 统计一律排除 excluded 对照（Step1 里程碑的
fastp 保留率/Q30 均值保留原口径：修剪口径对对照同样成立，reads 已由 DEC-31
OK 提示行播报）。
"""

import glob
import os
import statistics
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


def els_summary(els, excluded=()):
    """ELS 播报串（v2.22.0/DEC-32；最低值带 Z 值 v2.31.0/DEC-41）：
    "均值 … · 方差 … · 最低 <样本> …（Z=±x.x）"（科学计数法；Z=最低值偏离
    均值的 SD 倍数，直观判断该样本是否离群）——排除对照样本：NTC reads 近 0，
    其 ELS 无统计意义且必然占据最低值（误读为文库复杂度不足）。方差取总体
    方差（除以 n，与 check_depth_cv 同口径）；无有效数值/SD=0（无法定标）时
    不附 Z；无有效数值返回 None"""
    vals = [(sm, v) for sm, v in (els or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
    if not vals:
        return None
    mean = sum(v for _, v in vals) / len(vals)
    var = sum((v - mean) ** 2 for _, v in vals) / len(vals)
    lo_sm, lo_v = min(vals, key=lambda x: x[1])
    z_txt = ""
    sd = var ** 0.5
    if sd > 0:
        z_txt = f"（Z={(lo_v - mean) / sd:+.1f}）"
    return f"均值 {mean:.2e} · 方差 {var:.2e} · 最低 {lo_sm} {lo_v:.2e}{z_txt}"


def dup_stats_summary(dup_pct, excluded=()):
    """Step3 重复率分布播报（DEC-41）："min-max …%｜SD …%｜P25/P50/P75 …"——
    与批内变异 CV 共同刻画批次重复率形态（分位=statistics.quantiles 线性插值
    inclusive 口径）；对照排除；有效样本 <2 返回 None（调用方显示 n/a）"""
    vals = sorted(v for sm, v in (dup_pct or {}).items()
                  if isinstance(v, (int, float)) and sm not in excluded)
    if len(vals) < 2:
        return None
    mean = sum(vals) / len(vals)
    sd = (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5
    q1, q2, q3 = statistics.quantiles(vals, n=4, method="inclusive")
    return (f"min-max {vals[0]:.1f}-{vals[-1]:.1f}%｜SD {sd:.1f}%｜"
            f"P25/P50/P75 {q1:.1f}/{q2:.1f}/{q3:.1f}%")


def _violating(d, threshold, below=True):
    """→ 越界样本 [(样本, 值), ...]（名字序，全点名 DEC-30）；
    below=True 取 <阈值（越低越坏），False 取 >阈值（越高越坏）"""
    vals = sorted((k, v) for k, v in (d or {}).items()
                  if isinstance(v, (int, float)))
    return [(k, v) for k, v in vals
            if (v < threshold if below else v > threshold)]


def _only_samples(d, excluded):
    """剔除对照样本后的指标 dict（DEC-31）：样本级实验口径指标（保留率/mapped/
    dup/捕获/recal/深度 CV）对对照样本无统计意义，一律跳过——污染监控由
    check_ntc/check_ntc_reads 专属口径负责，不走这些检查"""
    if not excluded:
        return d or {}
    return {sm: v for sm, v in (d or {}).items() if sm not in excluded}


# ── 各步阈值检查（返回 anomalies 列表） ────────────────────────────────
def check_fastp(retention, q30, excluded=()):
    """越界样本全点名（DEC-30）：保留率/Q30 逐样本一行；
    80-90% 提示带为 OK 级，聚合一行点名（不逐行刷屏）；对照样本跳过（DEC-31）"""
    retention = _only_samples(retention, excluded)
    q30 = _only_samples(q30, excluded)
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


def check_flagstat(mapped_pct, pp_pct, excluded=()):
    mapped_pct = _only_samples(mapped_pct, excluded)
    pp_pct = _only_samples(pp_pct, excluded)
    out = []
    for sm, v in _violating(mapped_pct, config.MAPPED_NOTIFY_P1):
        out.append(("P2", f"{sm} mapped {v}%（阈值 {config.MAPPED_NOTIFY_P1}%）"))
    for sm, v in _violating(pp_pct, config.PROPER_PAIR_P1):
        out.append(("P2", f"{sm} properly paired {v}%（阈值 {config.PROPER_PAIR_P1}%）"))
    return out


def check_dup(dup_pct, excluded=()):
    return [("P2", f"{sm} 重复率 {v}%（阈值 {config.DUP_P1}%，建库复杂度告急）")
            for sm, v in _violating(_only_samples(dup_pct, excluded),
                                    config.DUP_P1, below=False)]


def check_capture(mean_depth, pct20x, pct_selected, excluded=()):
    """捕获效率口径 v2.19.0/DEC-29：PCT_SELECTED_BASES（on+near bait 占比对碱基比）
    为告警指标；on-target 不再告警（1bp SNP panel 下为几何产物，见 run_summary 信息指标）。
    越界样本全点名（DEC-30）；对照样本跳过（DEC-31）"""
    mean_depth = _only_samples(mean_depth, excluded)
    pct20x = _only_samples(pct20x, excluded)
    pct_selected = _only_samples(pct_selected, excluded)
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


def check_reads_low(reads, excluded=()):
    """reads: {sm: after_reads}——绝对量过低 → P1（上样不足，报错不中断，RUN-34）。
    对照样本（excluded，如 NTC）reads 接近 0 是正常状态：不按实验样本口径告警，
    降级为 OK 级提示行逐个播报实际数值（DEC-31）——reads 偏高属污染，由
    check_ntc_reads 占批次中位口径负责；数值缺失（None/非数值）不报"""
    out = []
    for sm, n in sorted((reads or {}).items()):
        if not isinstance(n, (int, float)) or n >= config.READS_MIN:
            continue
        if sm in excluded:
            out.append(("OK", f"{sm} reads {int(n)}（阴性对照，低 reads 属正常；"
                              f"reads 偏高属污染，由 NTC reads 占中位口径负责）"))
        else:
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


def cv_of(d, excluded=()):
    """样本间变异系数 CV（总体方差口径，与 check_depth_cv 一致）→ float|None
    （有效样本 <2 或无值 → None；对照排除，DEC-31）"""
    vals = [v for sm, v in (d or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
    if len(vals) < 2:
        return None
    mean = sum(vals) / len(vals)
    if mean <= 0:
        return None
    return (sum((v - mean) ** 2 for v in vals) / len(vals)) ** 0.5 / mean


def check_els_min_z(els, excluded=()):
    """ELS 最低样本 Z ≤ TH-39（-3 SD）→ P1：偏离批均值 3 个标准差以上，文库
    复杂度离群（该样本建库/上样环节异常）。总体 SD 口径下单点 Z 极值=-√(n-1)，
    批内 <10 样本数学上不可达（不属漏报）；SD=0 无法定标不判；对照排除"""
    vals = [(sm, v) for sm, v in (els or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
    if len(vals) < 2:
        return []
    mean = sum(v for _, v in vals) / len(vals)
    sd = (sum((v - mean) ** 2 for _, v in vals) / len(vals)) ** 0.5
    if sd <= 0:
        return []
    lo_sm, lo_v = min(vals, key=lambda x: x[1])
    z = (lo_v - mean) / sd
    if z <= config.ELS_MIN_Z_P1:
        return [("P1", f"ELS 最低样本 {lo_sm} Z={z:.1f}（阈值 "
                       f"{config.ELS_MIN_Z_P1:.0f} SD，ELS {lo_v:.2e} vs 批均值 "
                       f"{mean:.2e}）——文库复杂度离群，建议核查该样本建库/上样环节")]
    return []


def check_dup_cv(dup_pct, excluded=()):
    """批内重复率变异系数 CV 过大 → P1（TH-38=20%，DEC-40）：重复率是建库/上样
    环节的指纹，批内离散常指向文库质量或上样量差异——比 Step6 深度 CV（TH-35）
    早三步暴露批次异质；对照样本排除（NTC 重复率无统计意义）"""
    cv = cv_of(dup_pct, excluded)
    if cv is None or cv <= config.DUP_CV_P1:
        return []
    vals = [v for sm, v in (dup_pct or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
    return [("P1", f"批内重复率离散度 CV={cv * 100:.1f}%（阈值 "
                   f"{config.DUP_CV_P1 * 100:.0f}%，重复率 {min(vals):.1f}-{max(vals):.1f}%）"
                   f"——疑似批次异质（文库质量/上样量差异）")]


def check_depth_cv(mean_depth, excluded=()):
    """批次内样本间 mean depth 变异系数 CV 过大 → P2（疑似混入异常样本）；
    对照样本深度近 0 会拉爆 CV，不参与统计（DEC-31）"""
    vals = [v for sm, v in (mean_depth or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
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


def check_recal_low(obs, excluded=()):
    """recal: {sm: M 事件观测数}——known-sites 覆盖崩坏致校准不可信 → P2；
    对照样本跳过（DEC-31）"""
    out = []
    for sm, n in sorted(_only_samples(obs, excluded).items()):
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
                   anomalies=(), artifacts="", log_hint="", extra=(),
                   total_steps=None):
    """→ (title, text)。title 带 [GWAS][P级别]；total_steps 给定时标题尾部
    追加"（共 X 步）"（X=--step：Step 0 为清点预备步不计入，全流程=6——
    RUN-47 修正，首版误用 --step+1 显示 7）；正文含小写 gwas 兜底由 dingtalk 补。"""
    level = worst_level(anomalies)
    title = f"[GWAS][{level}] {batch}批次 · " \
            f"{f'Step {step_no} ' if step_no is not None else ''}{step_name}完成"
    if step_no is not None and total_steps:
        title += f"（共 {total_steps} 步）"
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
