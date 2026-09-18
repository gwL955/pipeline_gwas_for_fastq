#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""alerts.py 单测（P0/P1/P2 三级阈值判定、最差级别、里程碑模板字段、产物摘要）
分级语义（v2.9.0 DEC-21）：P0=阻断级（严重影响分析→中断批次）；
P1=严重（执行失败隔离/NTC 污染/mapped<90 QC 口径，报错不中断）；P2=质量提示"""

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
        self.assertEqual(alerts.worst_level([("P2", "x")]), "P2")
        self.assertEqual(alerts.worst_level([("P2", "x"), ("P1", "y")]), "P1")
        self.assertEqual(alerts.worst_level([("P1", "x"), ("P0", "y")]), "P0")
        self.assertEqual(alerts.worst_level([("OK", "x")]), "OK")


class TestChecks(unittest.TestCase):

    def test_fastp(self):
        a = alerts.check_fastp({"A": 75.0}, {"A": 90.0})
        self.assertEqual(a[0][0], "P2")
        self.assertIn("保留率", a[0][1])
        a = alerts.check_fastp({"A": 99.0}, {"A": 84.0})
        self.assertEqual(a[-1][0], "P2")           # Q30 <85 → P2（保留率提示行可能排前）
        self.assertIn("Q30", a[-1][1])
        self.assertEqual(alerts.check_fastp({"A": 99.0}, {"A": 95.0}), [])
        levels = [lv for lv, _ in alerts.check_fastp({"A": 85.0}, {"A": 95.0})]
        self.assertEqual(levels, ["OK"])            # 80-90 区间为提示行不升级（TH-02=90）

    def test_fastp_band_threshold_90(self):
        """★ RUN-45/DEC-31：提示带上沿 TH-02 由 95 调 90——93%（旧带内）不再
        进提示带，85% 仍在带内（260918 实测全批 92.18-94.47%，95 线全批点名
        失去区分度）"""
        self.assertEqual(alerts.check_fastp({"A": 93.0}, {"A": 95.0}), [])
        self.assertEqual(alerts.check_fastp({"A": 90.0}, {"A": 95.0}), [])   # 带不含 90 本身
        a = alerts.check_fastp({"A": 85.0}, {"A": 95.0})
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0][0], "OK")
        self.assertIn("A 85.0%", a[0][1])

    def test_fastp_all_violators_named(self):
        # 全点名（DEC-30）：多个越界样本逐个一行，不再只报最差一个
        a = alerts.check_fastp({"A": 75.0, "B": 79.5, "C": 99.0}, {"A": 95.0})
        self.assertEqual([lv for lv, _ in a], ["P2", "P2"])
        self.assertIn("A", a[0][1])
        self.assertIn("B", a[1][1])
        self.assertNotIn("C", a[0][1] + a[1][1])
        # 80-90 提示带：OK 级聚合一行点名全部带内样本
        a = alerts.check_fastp({"A": 82.0, "B": 83.5, "C": 99.0}, {"A": 95.0})
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0][0], "OK")
        self.assertIn("A 82.0%", a[0][1])
        self.assertIn("B 83.5%", a[0][1])

    def test_flagstat(self):
        a = alerts.check_flagstat({"A": 91.0}, {"A": 99.0})
        self.assertEqual(a[0][0], "P2")
        a = alerts.check_flagstat({"A": 99.0}, {"A": 84.0})
        self.assertEqual(a[0][0], "P2")            # pp <85
        self.assertEqual(alerts.check_flagstat({"A": 99.0}, {"A": 99.0}), [])
        # 全点名：mapped 与 pp 各自逐越界样本一行
        a = alerts.check_flagstat({"A": 93.0, "B": 91.0, "C": 99.0},
                                  {"A": 80.0, "B": 99.0})
        self.assertEqual(len(a), 3)
        self.assertIn("mapped", a[0][1])
        self.assertIn("properly paired", a[2][1])

    def test_dup(self):
        self.assertEqual(alerts.check_dup({"A": 31.0})[0][0], "P2")
        self.assertEqual(alerts.check_dup({"A": 5.0}), [])
        a = alerts.check_dup({"A": 31.0, "B": 45.0, "C": 5.0})
        self.assertEqual(len(a), 2)
        self.assertIn("A", a[0][1])
        self.assertIn("B", a[1][1])

    def test_capture(self):
        # 捕获效率口径 v2.19.0/DEC-29：PCT_SELECTED ≥85% 告警（新 panel 实测 89-91%，
        # 外送 pct_selected_bases 交叉验证）；on-target 已降为信息指标不再告警
        a = alerts.check_capture({"A": 45.0}, {"A": 96.0}, {"A": 90.0})
        self.assertEqual(a[0][0], "P2")
        self.assertIn("mean depth", a[0][1])
        a = alerts.check_capture({"A": 60.0}, {"A": 94.0}, {"A": 90.0})
        self.assertEqual(a[0][0], "P2")            # 20x <95
        a = alerts.check_capture({"A": 60.0}, {"A": 96.0}, {"A": 84.0})
        self.assertEqual(a[0][0], "P2")            # PCT_SELECTED <85
        self.assertIn("捕获效率", a[0][1])
        self.assertEqual(alerts.check_capture({"A": 80.0}, {"A": 97.0}, {"A": 85.0}), [])
        # on-target 任意低值（1bp SNP panel 下 ≈0.6%）不触发告警
        self.assertEqual(alerts.check_capture({"A": 80.0}, {"A": 97.0}, {"A": 90.0}), [])

    def test_capture_all_violators_named(self):
        # RUN-43 复盘锚：三样本 on-target 0.58/0.59/0.59 全部越界只报一个的历史行为
        # 不再复现——PCT_SELECTED 同场景三样本全点名
        a = alerts.check_capture(
            {"TG017": 136.2, "TG018": 194.7, "TG019": 193.1},
            {"TG017": 99.56, "TG018": 99.63, "TG019": 99.65},
            {"TG017": 40.0, "TG018": 41.5, "TG019": 89.42})
        self.assertEqual(len(a), 2)
        self.assertIn("TG017", a[0][1])
        self.assertIn("TG018", a[1][1])
        self.assertNotIn("TG019", a[0][1] + a[1][1])

    def test_variantqc(self):
        self.assertEqual(alerts.check_variantqc(1.9, 99.0)[0][0], "P2")
        self.assertEqual(alerts.check_variantqc(2.2, 90.0)[0][0], "P2")
        self.assertEqual(alerts.check_variantqc(2.2, 99.0), [])

    def test_ntc_p1(self):
        """NTC 污染：严重数据异常——P1 报错不中断（v2.9.0 由 P0 降级，DEC-21）"""
        a = alerts.check_ntc(15.0)
        self.assertEqual(a[0][0], "P1")
        self.assertIn("污染", a[0][1])
        self.assertEqual(alerts.check_ntc(0.01), [])
        self.assertEqual(alerts.check_ntc(None), [])


class TestNewChecks(unittest.TestCase):
    """RUN-34 新增检查：P1 reads 不足/空 PASS VCF；P2 NTC reads、深度 CV、recal 观测"""

    def test_reads_low_p1(self):
        a = alerts.check_reads_low({"A": 500000, "B": 5000000})
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0][0], "P1")
        self.assertIn("A", a[0][1])
        self.assertEqual(alerts.check_reads_low({"A": 5000000}), [])
        self.assertEqual(alerts.check_reads_low({}), [])

    def test_ntc_reads_p2(self):
        self.assertEqual(alerts.check_ntc_reads(1000, 1000000), [])   # 0.1% 正常
        a = alerts.check_ntc_reads(50000, 1000000)                    # 5% 异常
        self.assertEqual(a[0][0], "P2")
        self.assertIn("污染", a[0][1])
        self.assertEqual(alerts.check_ntc_reads(None, 1000), [])
        self.assertEqual(alerts.check_ntc_reads(50000, 0), [])        # 中位数无效

    def test_depth_cv_p2(self):
        self.assertEqual(alerts.check_depth_cv({"A": 100.0, "B": 102.0}), [])  # 均匀
        a = alerts.check_depth_cv({"A": 50.0, "B": 300.0})            # CV≈0.71
        self.assertEqual(a[0][0], "P2")
        self.assertIn("CV", a[0][1])
        self.assertEqual(alerts.check_depth_cv({"A": 100.0}), [])     # 单样本不判

    def test_recal_low_p2(self):
        a = alerts.check_recal_low({"A": 12345.0, "B": 45782114.0})
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0][0], "P2")
        self.assertIn("校准不可信", a[0][1])
        self.assertEqual(alerts.check_recal_low({"B": 45782114.0}), [])


class TestControlExemption(unittest.TestCase):
    """★ RUN-45/DEC-31：对照样本（excluded）豁免样本级阈值——reads 低降级 OK
    级提示行播报实际数值；其余指标直接跳过；污染仍由 check_ntc/check_ntc_reads
    专属口径负责（NTC reads 32 曾被"上样不足"P1 误报、保留率 73.37% 被 P2 误报）"""

    def test_reads_low_control_downgraded_to_ok(self):
        # 低 reads 的 NTC 不再 P1，而是 OK 级提示行且文案含实际数值；普通样本越界仍 P1
        a = alerts.check_reads_low({"NTC": 32, "A": 500000, "B": 5000000},
                                   excluded={"NTC"})
        p1 = [x for x in a if x[0] == "P1"]
        ok = [x for x in a if x[0] == "OK"]
        self.assertEqual(len(p1), 1)
        self.assertIn("A", p1[0][1])                # 普通样本照常 P1 全点名
        self.assertEqual(len(ok), 1)
        self.assertIn("NTC", ok[0][1])
        self.assertIn("32", ok[0][1])               # 播报实际数值
        self.assertIn("阴性对照", ok[0][1])
        # 不传 excluded（默认）→ 行为与旧版一致，NTC 仍 P1
        a = alerts.check_reads_low({"NTC": 32})
        self.assertEqual(a[0][0], "P1")
        # reads 数值缺失（None/非数值）的对照不报
        self.assertEqual(alerts.check_reads_low({"NTC": None}, excluded={"NTC"}), [])

    def test_reads_low_multiple_controls_each_named(self):
        # 多个对照样本逐个点名（DEC-30 名字序全点名风格）
        a = alerts.check_reads_low({"NTC1": 10, "NTC2": 20},
                                   excluded={"NTC1", "NTC2"})
        self.assertEqual([lv for lv, _ in a], ["OK", "OK"])
        self.assertIn("NTC1", a[0][1])
        self.assertIn("NTC2", a[1][1])

    def test_fastp_excluded_control_skipped(self):
        # 对照保留率/Q30 不进 P2 越界点名、不进提示带统计；普通样本带内照常提示
        a = alerts.check_fastp({"NTC": 73.37, "A": 85.0}, {"NTC": 70.0, "A": 95.0},
                               excluded={"NTC"})
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0][0], "OK")             # 仅 A 85% 提示带聚合行
        self.assertNotIn("NTC", a[0][1])
        self.assertEqual(alerts.check_fastp({"NTC": 73.37}, {"NTC": 70.0},
                                            excluded={"NTC"}), [])

    def test_depth_cv_excluded_control_skipped(self):
        # NTC 深度近 0 会拉爆批次 CV，对照不参与统计；真实离散仍报
        self.assertEqual(alerts.check_depth_cv({"NTC": 0.5, "A": 100.0, "B": 102.0},
                                               excluded={"NTC"}), [])
        a = alerts.check_depth_cv({"NTC": 0.5, "A": 50.0, "B": 300.0},
                                  excluded={"NTC"})
        self.assertEqual(a[0][0], "P2")
        # 不传 excluded（默认）→ NTC 参与统计拉爆 CV（旧行为保持）
        self.assertEqual(alerts.check_depth_cv({"NTC": 0.5, "A": 100.0, "B": 102.0})[0][0],
                         "P2")

class TestElsSummary(unittest.TestCase):
    """★ RUN-46/DEC-32：Step3 ELS 播报排除对照并扩为 均值/方差/最低值+样本——
    NTC reads 近 0，其 ELS 无统计意义且必然占据最低值（曾误读为文库复杂度不足）"""

    def test_control_excluded_from_els_stats(self):
        # 不排除：NTC(3) 占据最低且拉低均值（旧误判形态）；排除后恢复实验样本口径
        els = {"TG017": 8305698, "TG018": 9100000, "TG019": 8500000, "NTC": 3}
        with_ctl = alerts.els_summary(els)
        self.assertIn("NTC", with_ctl)
        s = alerts.els_summary(els, excluded={"NTC"})
        self.assertNotIn("NTC", s)
        self.assertIn("均值", s)
        self.assertIn("方差", s)
        self.assertIn("最低 TG017", s)               # 最低值带对应样本名
        self.assertIn("8.31e+06", s)                 # TG017 最低值本体

    def test_els_population_variance(self):
        # 方差=总体方差（÷n，与 check_depth_cv 同口径）：{100,200} → 均值150、方差2500
        s = alerts.els_summary({"A": 100, "B": 200})
        self.assertIn("1.50e+02", s)
        self.assertIn("2.50e+03", s)
        self.assertIn("最低 A", s)

    def test_els_no_valid_values(self):
        # 空/全对照/非数值（None 是 GATK "?" 的解析产物）→ None（调用方显示 n/a）
        self.assertIsNone(alerts.els_summary(None))
        self.assertIsNone(alerts.els_summary({}))
        self.assertIsNone(alerts.els_summary({"NTC": 3}, excluded={"NTC"}))
        self.assertIsNone(alerts.els_summary({"A": None, "B": "?"}))


class TestBroadcastExemption(unittest.TestCase):
    """★ RUN-46/DEC-32：指标播报均值排除对照（run_pipeline._avg）——DEC-31 只
    豁免了告警点名，"指标:"行均值此前仍含 NTC（mapped/深度等被 reads 近 0 拉低）"""

    def test_avg_excludes_control(self):
        import run_pipeline
        d = {"A": 98.0, "B": 99.0, "NTC": 2.0}
        self.assertEqual(run_pipeline._avg(d), round((98.0 + 99.0 + 2.0) / 3, 2))
        self.assertEqual(run_pipeline._avg(d, {"NTC"}), 98.5)
        self.assertIsNone(run_pipeline._avg({}, {"NTC"}))


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

    def test_total_steps_in_title(self):
        """★ RUN-46/DEC-32 + RUN-47：里程碑标题增"（共 X 步）"，X=--step——
        Step 0 为清点预备步不计入，全流程共 6 步（首版误用 +1 显示 7）；
        --step 0（只跑清点）时 0 为假值不追加后缀；不传保持旧标题（兼容）"""
        title, text = alerts.step_milestone("260918", 3, "去重", total_steps=6)
        self.assertEqual(title, "[GWAS][OK] 260918批次 · Step 3 去重完成（共 6 步）")
        self.assertIn("#### [GWAS][OK] 260918批次 · Step 3 去重完成（共 6 步）", text)
        t_zero, _ = alerts.step_milestone("260918", 0, "清点与Lane合并", total_steps=0)
        self.assertEqual(t_zero, "[GWAS][OK] 260918批次 · Step 0 清点与Lane合并完成")
        t_no, _ = alerts.step_milestone("260918", 3, "去重")
        self.assertEqual(t_no, "[GWAS][OK] 260918批次 · Step 3 去重完成")

    def test_level_in_title(self):
        _, text_p0 = alerts.step_milestone("B", 6, "质量汇总",
                                           anomalies=[("P0", "样本名非法")])
        self.assertIn("← 阻断级（中断分析）", text_p0)
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
