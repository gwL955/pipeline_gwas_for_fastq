#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""运行报告（markdown）与 run_summary.json 生成。仅标准库。"""

import json
import os
from datetime import datetime

import config


# ── run_summary.json ────────────────────────────────────────────────────
def write_run_summary(summary, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    summary["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)


# ── 运行报告（results/<batch>_<date>/run_report.md） ────────────────────
def render_run_report(batch, bdata, out_path):
    m = bdata.get("metrics", {}) or {}
    lines = [f"# {batch} 运行报告", "",
             f"- 状态：{bdata.get('status')}",
             f"- 有效样本：{', '.join(bdata.get('samples', {}).get('valid', []))}"
             f"（无效/跳过：{bdata.get('samples', {}).get('invalid') or '无'}）",
             f"- 联合分型样本（排除对照）："
             f"{', '.join(bdata.get('samples', {}).get('calling', []))}", "",
             "## 关键指标", "",
             "| 指标 | 值 |", "| --- | --- |"]
    for key, label in (("fastp_retention", "fastp 保留率 %"),
                       ("mapped_pct", "mapped %"), ("dup_pct", "dup %"),
                       ("mean_target_coverage", "MEAN 靶深度 ×"),
                       ("pct_20x", "20X 覆盖 %"),
                       ("pct_selected", "捕获效率 %（PCT_SELECTED）"),
                       ("on_target_pct", "on-target %（信息指标）"),
                       ("snp_raw", "raw SNP"), ("indel_raw", "raw INDEL"),
                       ("snp_pass", "PASS SNP"), ("indel_pass", "PASS INDEL"),
                       ("titv_raw", "Ti/Tv raw"), ("titv_pass", "Ti/Tv PASS")):
        lines.append(f"| {label} | {m.get(key)} |")
    art = bdata.get("artifacts", {})
    if art:
        lines += ["", "## 主要产物", ""]
        for k, v in art.items():
            lines.append(f"- {k}: `{v}`")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
