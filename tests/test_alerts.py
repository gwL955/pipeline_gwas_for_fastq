#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""alerts.py 单测（P0/P1 阈值判定、最差级别、里程碑模板字段、产物摘要）"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import alerts                                     # noqa: E402
import config                                     # noqa: E402


class TestLevels(unittest.TestCase):

    def test_worst_level(self):
        self.assertEqual(alerts.worst_level([]), "OK")
        self.assertEqual(alerts.worst_level([("P1", "x")]), "P1")
        self.assertEqual(alerts.worst_level([("P1", "x"), ("P0", "y")]), "P0")
        self.assertEqual(alerts.worst_level([("OK", "x")]), "OK")


class TestChecks(unittest.TestCase):

    def test_fastp(self):
        a = alerts.check_fastp({"A": 75.0}, {"A": 90.0})
        self.assertEqual(a[0][0], "P1")
        self.assertIn("保留率", a[0][1])
        a = alerts.check_fastp({"A": 99.0}, {"A": 84.0})
        self.assertIn(("P1", a[-1][1]), a)          # Q30 <85 → P1（85-95 保留率提示行可能排前）
        self.assertEqual(alerts.check_fastp({"A": 99.0}, {"A": 95.0}), [])
        levels = [lv for lv, _ in alerts.check_fastp({"A": 90.0}, {"A": 95.0})]
        self.assertEqual(levels, ["OK"])            # 80-95 区间为提示行不升级

    def test_flagstat(self):
        a = alerts.check_flagstat({"A": 91.0}, {"A": 99.0})
        self.assertEqual(a[0][0], "P1")
        a = alerts.check_flagstat({"A": 99.0}, {"A": 84.0})
        self.assertEqual(a[0][0], "P1")            # pp <85
        self.assertEqual(alerts.check_flagstat({"A": 99.0}, {"A": 99.0}), [])

    def test_dup(self):
        self.assertEqual(alerts.check_dup({"A": 31.0})[0][0], "P1")
        self.assertEqual(alerts.check_dup({"A": 5.0}), [])

    def test_capture(self):
        # on-target 阈值 8%（小 panel 正常 8-10%，60% 为误用 off-bait 口径的历史值）
        a = alerts.check_capture({"A": 45.0}, {"A": 96.0}, {"A": 9.0})
        self.assertEqual(a[0][0], "P1")
        self.assertIn("mean depth", a[0][1])
        a = alerts.check_capture({"A": 60.0}, {"A": 94.0}, {"A": 9.0})
        self.assertEqual(a[0][0], "P1")            # 20x <95
        a = alerts.check_capture({"A": 60.0}, {"A": 96.0}, {"A": 7.0})
        self.assertEqual(a[0][0], "P1")            # on-target <8
        self.assertIn("on-target", a[0][1])
        self.assertEqual(alerts.check_capture({"A": 80.0}, {"A": 97.0}, {"A": 8.5}), [])
        self.assertEqual(alerts.check_capture({"A": 80.0}, {"A": 97.0}, {"A": 70.0}), [])

    def test_variantqc(self):
        self.assertEqual(alerts.check_variantqc(1.9, 99.0)[0][0], "P1")
        self.assertEqual(alerts.check_variantqc(2.2, 90.0)[0][0], "P1")
        self.assertEqual(alerts.check_variantqc(2.2, 99.0), [])

    def test_ntc_p0(self):
        a = alerts.check_ntc(15.0)
        self.assertEqual(a[0][0], "P0")
        self.assertIn("污染", a[0][1])
        self.assertEqual(alerts.check_ntc(0.01), [])
        self.assertEqual(alerts.check_ntc(None), [])


class TestMilestone(unittest.TestCase):

    def test_template_fields(self):
        title, text = alerts.step_milestone(
            "20260720", 2, "比对",
            samples="4/4 成功 | Lane 合并 16/16",
            metrics="mapped 98.7% | proper pair 94.2%",
            anomalies=[("P1", "L20260615001 mapped 91.3%（阈值 95%）")],
            artifacts="bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22",
            log_hint="tail -f logs/sample_L20260615001.log")
        self.assertEqual(title, "[GWAS][P1] 20260720批次 · Step 2 比对完成")
        for kw in ("样本: 4/4 成功", "指标: mapped 98.7%",
                   "异常: [P1] L20260615001", "← 需确认",
                   "产物: bam/*/*.sort.bam ×4", "日志: tail -f"):
            self.assertIn(kw, text)

    def test_level_in_title(self):
        _, text_p0 = alerts.step_milestone("B", 6, "质量汇总",
                                           anomalies=[("P0", "NTC 污染")])
        self.assertIn("← 阻断/污染级", text_p0)
        t_ok, text_ok = alerts.step_milestone("B", 3, "去重", anomalies=[])
        self.assertIn("[OK]", t_ok)
        self.assertIn("异常: 无", text_ok)

    def test_artifact_summary(self):
        with tempfile.TemporaryDirectory() as td:
            for i in range(3):
                with open(os.path.join(td, f"f{i}.bam"), "wb") as f:
                    f.write(b"x")
            s = alerts.artifact_summary(os.path.join(td, "*.bam"))
            self.assertTrue(s.startswith("×3  mtime "))
            self.assertEqual(alerts.artifact_summary(os.path.join(td, "*.none")), "×0")


if __name__ == "__main__":
    unittest.main()
