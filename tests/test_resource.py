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


if __name__ == "__main__":
    unittest.main()
