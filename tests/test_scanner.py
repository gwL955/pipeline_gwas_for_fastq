#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scanner.py 单测（三种输入布局解析、无效样本标记、布局冲突、md5 校验、Lane 合并与 dry-run 行为）"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scanner                                   # noqa: E402
from runner import Runner                       # noqa: E402


def _touch(path, data=b"@read\nACGT\n+\nIIII\n"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(data)


class TestScanLayouts(unittest.TestCase):

    def test_illumina_flat_layout(self):
        with tempfile.TemporaryDirectory() as td:
            for lane in ("L001", "L002"):
                for r in ("R1", "R2"):
                    _touch(os.path.join(td, f"NA12878_S46_{lane}_{r}_001.fastq.gz"))
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(list(valid), ["NA12878"])
            self.assertEqual(valid["NA12878"].snum, "S46")
            self.assertEqual(valid["NA12878"].lanes, ["001", "002"])   # 存数字串
            self.assertEqual(len(valid["NA12878"].r1), 2)
            self.assertEqual(valid["NA12878"].layout, "illumina")
            self.assertEqual(invalid, {})

    def test_outsourced_subdir_layout(self):
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "L200100435", "L200100435_R1.fastq.gz"))
            _touch(os.path.join(td, "L200100435", "L200100435_R2.fastq.gz"))
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(list(valid), ["L200100435"])
            self.assertEqual(valid["L200100435"].layout, "outsourced")

    def test_outsourced_flat_layout(self):
        """★ 外送平铺布局（RUN-36/DEC-23）：<样本>_R{1,2}.fastq.gz 直接放批次目录。
        真实事故：外送交付平铺数据两种识别器都不认 → 整批"无有效样本"静默跳过"""
        with tempfile.TemporaryDirectory() as td:
            sm = "SURFS45391-260810-HS02-U241-20260811cap18ZFD-TG018"   # 真实交付名截取
            for r in ("R1", "R2"):
                _touch(os.path.join(td, f"{sm}_{r}.fastq.gz"))
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(list(valid), [sm])
            self.assertEqual(valid[sm].layout, "outsourced_flat")
            self.assertEqual(len(valid[sm].r1), 1)
            self.assertEqual(len(valid[sm].r2), 1)
            self.assertIn("外送平铺", valid[sm].note())
            self.assertEqual(invalid, {})

    def test_outsourced_flat_r1_only_invalid(self):
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM_R1.fastq.gz"))          # 缺 R2
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(valid, {})
            self.assertIn("不一致", invalid["SM"])

    def test_illumina_name_not_mistaken_as_flat(self):
        """互斥锚：Illumina 命名（尾 _001.fastq.gz）不得同时被外送平铺识别器认出"""
        with tempfile.TemporaryDirectory() as td:
            for lane in ("L001", "L002"):
                for r in ("R1", "R2"):
                    _touch(os.path.join(td, f"SM_S1_{lane}_{r}_001.fastq.gz"))
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(list(valid), ["SM"])               # 唯一样本，无第二识别结果
            self.assertEqual(valid["SM"].layout, "illumina")
            self.assertEqual(invalid, {})

    def test_layout_conflict_marks_invalid(self):
        """同一样本名被多种布局认出 → 先认者优先，后到者记冲突无效（历史行为保持）"""
        with tempfile.TemporaryDirectory() as td:
            for r in ("R1", "R2"):
                _touch(os.path.join(td, f"SM_S1_L001_{r}_001.fastq.gz"))
                _touch(os.path.join(td, f"SM_{r}.fastq.gz"))    # 外送平铺同文名
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(list(valid), ["SM"])               # illumina 先认
            self.assertIn("冲突", invalid["SM"])

    def test_invalid_r1r2_mismatch(self):
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM_S1_L001_R1_001.fastq.gz"))   # 只有 R1
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(valid, {})
            self.assertIn("不一致", invalid["SM"])

    def test_invalid_zero_byte(self):
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM_S1_L001_R1_001.fastq.gz"), b"")
            _touch(os.path.join(td, "SM_S1_L001_R2_001.fastq.gz"))
            valid, invalid = scanner.scan_batch(td)
            self.assertEqual(valid, {})
            self.assertIn("0 字节", invalid["SM"])

    def test_undetermined_ignored(self):
        """★ RUN-45/DEC-31：Illumina 下机自带 Undetermined（BCLConvert 未匹配
        index reads）不得作为样本进入分析——曾一路进联合分型/基因型矩阵/Output
        污染整批。剔除后进 ignored 清单：不算 invalid、不触发告警；
        既有 `valid, invalid = scan_batch(...)` 二元组解包不受影响"""
        with tempfile.TemporaryDirectory() as td:
            for r in ("R1", "R2"):
                _touch(os.path.join(td, f"SM_S1_L001_{r}_001.fastq.gz"))
                _touch(os.path.join(td, f"Undetermined_S0_L001_{r}_001.fastq.gz"))
            res = scanner.scan_batch(td)
            valid, invalid = res                    # 二元组解包兼容锚
            self.assertEqual(list(valid), ["SM"])
            self.assertNotIn("Undetermined", valid)
            self.assertNotIn("Undetermined", invalid)   # 不算 invalid
            self.assertIn("Undetermined", res.ignored)  # 进独立 ignored 清单

    def test_undetermined_match_case_insensitive(self):
        """剔除规则不区分大小写（IGNORED_SAMPLES 常量，便于日后扩充）"""
        with tempfile.TemporaryDirectory() as td:
            for r in ("R1", "R2"):
                _touch(os.path.join(td, f"undetermined_S0_L001_{r}_001.fastq.gz"))
            res = scanner.scan_batch(td)
            self.assertEqual(res[0], {})
            self.assertEqual(res.ignored, ["undetermined"])


class TestMd5(unittest.TestCase):

    def test_verify_md5_pass_and_fail(self):
        import hashlib
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM1", "SM1_R1.fastq.gz"), b"aaaa")
            _touch(os.path.join(td, "SM1", "SM1_R2.fastq.gz"), b"bbbb")
            good = hashlib.md5(b"aaaa").hexdigest()
            bad = "0" * 32
            with open(os.path.join(td, "md5sum.txt"), "w", encoding="utf-8") as f:
                f.write(f"{good}  SM1/SM1_R1.fastq.gz\n{bad}  SM1/SM1_R2.fastq.gz\n")

            class _L:
                def info(self, *a): pass
                def error(self, *a): pass
                def result(self, *a): pass
                def warn(self, *a): pass
            failed, n = scanner.verify_md5(td, _L(), workers=2)
            self.assertEqual(n, 2)
            self.assertIn("SM1", failed)      # 坏的那条 → 样本列入失败集

    def test_verify_md5_flat_outsourced_sample_name(self):
        """★ md5 失败样本名推导覆盖外送平铺（RUN-36）：平铺文件样本名=去 _R# 尾——
        曾误取整个文件名，致 md5_failed 与 valid 无交集（P0 失效，同 RUN-33 病根）"""
        import hashlib
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM1_R1.fastq.gz"), b"aaaa")
            _touch(os.path.join(td, "SM1_R2.fastq.gz"), b"bbbb")
            bad, good = "0" * 32, hashlib.md5(b"bbbb").hexdigest()
            with open(os.path.join(td, "md5sum.txt"), "w", encoding="utf-8") as f:
                f.write(f"{bad}  SM1_R1.fastq.gz\n{good}  SM1_R2.fastq.gz\n")

            class _L:
                def info(self, *a): pass
                def error(self, *a): pass
                def result(self, *a): pass
                def warn(self, *a): pass
            failed, n = scanner.verify_md5(td, _L(), workers=2)
            self.assertEqual(n, 2)
            self.assertIn("SM1", failed)
            self.assertNotIn("SM1_R1.fastq.gz", failed)         # 不得是整个文件名


class TestMerge(unittest.TestCase):

    def test_merge_all_real(self):
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM_S1_L001_R1_001.fastq.gz"), b"a")
            _touch(os.path.join(td, "SM_S1_L001_R2_001.fastq.gz"), b"b")
            valid, _ = scanner.scan_batch(td)
            out = os.path.join(td, "merged")
            ok, failed = scanner.merge_all(valid, out, Runner(), lambda sm: None, workers=2)
            self.assertEqual(failed, {})
            r1, r2 = scanner.merged_paths(out, "SM")
            self.assertTrue(os.path.isfile(r1) and os.path.getsize(r1) > 0)

    def test_merge_dry_run_returns_ok_without_files(self):
        """踩坑回归：dry-run 下合并不落文件，曾被误判'合并失败'导致批次失败"""
        with tempfile.TemporaryDirectory() as td:
            _touch(os.path.join(td, "SM_S1_L001_R1_001.fastq.gz"), b"a")
            _touch(os.path.join(td, "SM_S1_L001_R2_001.fastq.gz"), b"b")
            valid, _ = scanner.scan_batch(td)
            out = os.path.join(td, "merged")
            ok, failed = scanner.merge_all(valid, out, Runner(dry_run=True),
                                           lambda sm: None, workers=1)
            self.assertIn("SM", ok)
            self.assertEqual(failed, {})

    def test_samples_tsv_columns(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "samples.tsv")
            merged = {"SM": ("/x/SM_R1.fastq.gz", "/x/SM_R2.fastq.gz")}
            info = scanner.SampleInfo("SM", "illumina")
            info.snum, info.lanes = "S1", ["L001"]
            scanner.write_samples_tsv(p, merged, {"SM": info})
            with open(p, encoding="utf-8") as f:
                lines = f.read().strip().split("\n")
            self.assertEqual(lines[0].split("\t")[:3], ["sample", "R1", "R2"])
            self.assertEqual(len(lines[1].split("\t")), 4)   # +note 追溯列


if __name__ == "__main__":
    unittest.main()
