#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scanner.py 单测（两种输入布局解析、无效样本标记、md5 校验、Lane 合并与 dry-run 行为）"""

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
