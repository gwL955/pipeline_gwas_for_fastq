#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""resource.py 单测（对应踩坑：samtools sort -m 只接受整数+单位——'2.0G' 被解析为
2 字节；以及低配快速失败、档位推导、计划透明）"""

import os
import re
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import resource                                  # noqa: E402


class TestSortMemFormat(unittest.TestCase):
    """★ 核心踩坑：sort -m 必须是 '128M'/'2048M' 这类整数+单位"""

    def _plan(self, profile="auto", **kw):
        """auto 档用真实内存推导：小内存机器（CI runner ~16G）会因单样本峰值
        19.5G 触发快速失败 SystemExit（--max-memory 是上限压不高它）——统一
        mock 探测为 32 线程/64G，测试只验证推导逻辑本身"""
        with mock.patch.object(resource, "detect_cpu", return_value=32), \
                mock.patch.object(resource, "detect_mem_gb", return_value=64.0):
            return resource.plan(profile=profile, **kw)

    def test_sort_mem_integer_megabytes(self):
        for profile in ("auto", "low", "high"):
            with self.subTest(profile=profile):
                p = self._plan(profile=profile)
                self.assertRegex(p.sort_mem, r"^\d+M$",
                                 f"{profile} 档 sort_mem={p.sort_mem} 含小数点/其他单位")

    def test_gatk_mem_format(self):
        p = self._plan()
        self.assertRegex(p.gatk_mem, r"^\d+g$")
        self.assertRegex(p.cohort_mem, r"^\d+g$")

    def test_plan_invariants(self):
        p = self._plan()
        self.assertGreaterEqual(p.workers, 1)
        self.assertGreaterEqual(p.threads, 4)
        self.assertIn(p.threads, (4, 8, 12, 24))       # 分档表
        self.assertGreater(p.peak_per_sample_gb, 0)
        self.assertEqual(p.workers * p.threads, p.workers * p.threads)  # int


class TestProfiles(unittest.TestCase):

    def test_low_profile_conservative(self):
        """低配档（16 线程/20G）：workers=1、sort 128M、GATK 1g（提示词硬性口径）"""
        p = resource.plan(profile="low")
        self.assertEqual(p.workers, 1)
        self.assertEqual(p.threads, 4)
        self.assertEqual(p.sort_mem, "128M")
        self.assertEqual(p.gatk_mem, "1g")

    def test_high_profile(self):
        p = resource.plan(profile="high")
        self.assertEqual(p.workers, 4)
        self.assertEqual(p.threads, 24)
        self.assertEqual(p.gatk_mem, "8g")
        self.assertLessEqual(int(p.sort_mem.rstrip("M")), 2048)

    def test_fastfail_when_index_does_not_fit(self):
        """快速失败：可用内存装不下 bwa-mem2 索引（17G）时启动即报错退出"""
        with mock.patch.object(resource, "detect_cpu", return_value=8), \
                mock.patch.object(resource, "detect_mem_gb", return_value=10.0):
            with self.assertRaises(SystemExit):
                resource.plan(profile="auto")

    def test_plan_table_transparent(self):
        p = resource.plan(profile="low")
        t = p.table()
        for kw in ("探测", "预留", "workers", "sort -m", "GATK -Xmx", "峰值"):
            self.assertIn(kw, t)

    def test_workers_override_wins(self):
        p = resource.plan(profile="low", workers_override=1)
        self.assertEqual(p.workers, 1)

    def test_to_dict_json_serializable(self):
        import json
        p = resource.plan(profile="low")
        json.dumps(p.to_dict())


class TestWorkerClasses(unittest.TestCase):
    """★ 步骤类型分化 workers（DEC-35/RUN-50）：GATK 单线程工具以比对类
    workers（每路含 17G bwa 索引峰值）并行时整机 CPU ~10% 严重空转——
    48 样本批次 6h53m 中 Step3/4/5/6 全被 workers=4 钉死；分化后 50 核机
    GATK 类 12 路（12×hmm4=48 核满载不超订），比对类维持 bwa 饱和设计"""

    def _plan50(self, **kw):
        """事故机规格（50 线程/305.7GB，auto）"""
        with mock.patch.object(resource, "detect_cpu", return_value=50), \
                mock.patch.object(resource, "detect_mem_gb", return_value=305.7):
            return resource.plan(**kw)

    def test_50core_gatk_class_12way(self):
        p = self._plan50()
        self.assertEqual(p.workers, 4)            # 比对类不变（索引峰值口径）
        self.assertEqual(p.gatk_mem, "8g")
        self.assertEqual(p.workers_gatk, 12)      # min(48//4, 290//10.4)
        self.assertEqual(p.workers_gatk_cpu, 12)
        self.assertEqual(p.workers_gatk_mem, 27)
        self.assertEqual(p.workers_io, 12)        # 内存宽裕 → 与 GATK 类同路上限
        self.assertEqual(p.hc_hmm_threads, 4)     # 与现行档一致（单样本耗时不变）

    def test_no_oversubscription(self):
        """HC 防超订阅不变式：workers_gatk × hmm ≤ usable_cores（全档位）"""
        for p in (self._plan50(), resource.plan(profile="low"),
                  resource.plan(profile="high"), self._plan50(workers_override=20)):
            with self.subTest(profile=p.profile, ov=p.workers_gatk):
                self.assertLessEqual(p.workers_gatk * p.hc_hmm_threads,
                                     p.usable_cores)

    def test_low_profile_gatk_class(self):
        """低配档（16 线程/20G）：比对类仍 workers=1（索引口径），GATK 类放宽
        到 7 路（14//2，内存 18//1.3=13 不设限）——GATK 阶段索引非工作集"""
        p = resource.plan(profile="low")
        self.assertEqual(p.workers, 1)
        self.assertEqual(p.workers_gatk, 7)
        self.assertEqual(p.workers_io, 7)
        self.assertEqual(p.hc_hmm_threads, 2)     # min(2 档, 14//7)

    def test_high_profile_gatk_class(self):
        p = resource.plan(profile="high")
        self.assertEqual(p.workers, 4)
        self.assertEqual(p.workers_gatk, 24)      # min(96//4, 855//10.4)
        self.assertEqual(p.hc_hmm_threads, 4)

    def test_workers_override_applies_to_all_classes(self):
        """--workers 最高优先级：三类 workers 同步覆盖，hmm 随之反推防超订阅"""
        p = self._plan50(workers_override=20)
        self.assertEqual(p.workers, 20)
        self.assertEqual(p.workers_gatk, 20)
        self.assertEqual(p.workers_io, 20)
        self.assertEqual(p.hc_hmm_threads, 2)     # min(4 档, 48//20)

    def test_table_shows_class_rows(self):
        t = self._plan50().table()
        for kw in ("workers 比对类", "workers GATK 类", "workers IO 类",
                   "防超订阅", "不含 bwa 索引"):
            self.assertIn(kw, t)

    def test_to_dict_carries_class_keys(self):
        import json
        d = self._plan50().to_dict()
        for k in ("workers_gatk", "workers_io", "workers_gatk_cpu",
                  "workers_gatk_mem", "gatk_peak_gb"):
            self.assertIn(k, d)
        json.dumps(d)


if __name__ == "__main__":
    unittest.main()
