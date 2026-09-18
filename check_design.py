#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""设计文档一致性校验（防漂移）：
  1) design_doc/DESIGN.md 的 TH 阈值表必须与 config.py 实际值一致
  2) REQ 表不允许存在未回写状态的行（人新增需求后机器必须实现并回写）
用法：python3 pipeline/check_design.py   （run_tests.sh 末尾自动执行）
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config                                       # noqa: E402

DESIGN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "design_doc", "DESIGN.md")


def _parse_value(raw):
    raw = raw.strip().strip("`")
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
            return raw


def main():
    if not os.path.isfile(DESIGN):
        print(f"[FAIL] 设计文档不存在: {DESIGN}")
        return 1
    text = open(DESIGN, encoding="utf-8").read()
    errors, warnings = [], []

    # ① TH 表 ↔ config.py
    th_rows = re.findall(r"^\|\s*(TH-\d+)\s*\|\s*([A-Z_0-9]+)\s*\|\s*([^|]+)\|", text,
                         re.MULTILINE)
    if not th_rows:
        errors.append("DESIGN.md 中未找到 TH 阈值表（格式：| TH-xx | CONFIG_KEY | value |）")
    for th_id, key, val in th_rows:
        if not hasattr(config, key):
            errors.append(f"{th_id}: config.py 无属性 {key}")
            continue
        actual, declared = getattr(config, key), _parse_value(val)
        same = (abs(actual - declared) < 1e-9
                if isinstance(actual, (int, float)) and isinstance(declared, (int, float))
                else str(actual) == str(declared))
        if not same:
            errors.append(f"{th_id}: config.{key}={actual!r} ≠ 文档值 {declared!r}")

    # ② config.py 中可镜像的标量阈值是否有 TH 行覆盖（防一侧新增未同步）
    mirrored = {"FASTP_LENGTH_REQUIRED", "FASTP_RETENTION_WARN", "FASTP_RETENTION_P1",
                "FASTP_Q30_P1", "MAPPED_MIN_PCT", "MAPPED_NOTIFY_P1", "PROPER_PAIR_P1",
                "DUP_WARN_PCT", "DUP_P1", "MEAN_COV_MIN", "MEAN_DEPTH_P1", "PCT_20X_MIN",
                "PCT_20X_P1", "PCT_SELECTED_P1", "DP_MIN", "HC_INTERVAL_PADDING",
                "NTC_DEPTH_P0", "TITV_P1", "CALL_RATE_P1", "DISK_MIN_FREE_GB",
                "DISK_PER_SAMPLE_GB", "BWA_INDEX_MEM_GB", "SORT_MEM_MIN", "SORT_MEM_MAX",
                "GATK_MEM_MIN_GB", "GATK_MEM_MAX_GB", "COHORT_MEM_MIN_GB",
                "COHORT_MEM_MAX_GB", "READS_MIN", "NTC_READS_PCT_P2",
                "DEPTH_CV_P2", "RECAL_OBS_MIN_P2"}
    documented = {key for _, key, _ in th_rows}
    for key in sorted(mirrored - documented):
        warnings.append(f"config.{key} 未在 DESIGN.md TH 表登记（若属新阈值请补 TH-xx 行")

    # ③ REQ 状态完整性（人新增目标 → 机器必须实现并回写状态）
    for line in text.splitlines():
        if not re.match(r"^\|\s*REQ-\d+\s*\|", line):
            continue
        cols = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cols) >= 3 and cols[0].startswith("REQ-") and not cols[2]:
            errors.append(f"{cols[0]}: 需求已提出但状态列为空（机器实现后须回写 §2 状态与证据）")

    # ④ 版本锚点存在 + 与 config.PIPELINE_VERSION 同步（交付/run_summary 溯源依赖）
    m_ver = re.search(r"^version:\s*(\d+\.\d+\.\d+)", text, re.MULTILINE)
    if not m_ver:
        errors.append("design-meta 缺少 version 字段（人机共写协议要求）")
    elif m_ver.group(1) != config.PIPELINE_VERSION:
        errors.append(f"版本漂移：DESIGN.md version={m_ver.group(1)} ≠ "
                      f"config.PIPELINE_VERSION={config.PIPELINE_VERSION}（须两处同步）")

    for w in warnings:
        print(f"[WARN] {w}")
    if errors:
        for e in errors:
            print(f"[FAIL] {e}")
        print(f"check_design: {len(errors)} 项不一致 —— 请先同步 DESIGN.md 与 config.py 再提交")
        return 1
    print(f"check_design: OK（TH×{len(th_rows)} 与 config.py 一致，REQ 状态完整，版本锚点在）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
