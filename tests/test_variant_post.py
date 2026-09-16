#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""变异后处理核心逻辑单测（对应踩坑：GT 列序必须与 VCF 样本序一致、
rebuild 只换 GT 其余字段保留、矩阵 ./. 裁决阈值）
（比对相关用例已随比对模块移除，v2.0.0）"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules import bcftools as mbc               # noqa: E402


class TestAdjudicate(unittest.TestCase):
    """矩阵 ./. 裁决：DP≥20 改判 0/0，不足保留；返回总格子数供 call rate"""

    MATRIX = ("chr1\t100\tA\tG\t0/1\t./.\t1/1\n"
              "chr1\t200\tC\tT\t./.\t0/0\t./.\n"
              "chr2\t300\tG\tA\t0/0\t0/1\t./.\n")

    def _run(self, depth):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "m.tsv")
            open(p, "w").write(self.MATRIX)
            out, stats = mbc.adjudicate_matrix(
                p, lambda sm, c, pos: depth, ["S1", "S2", "S3"], 20)
            return out, stats

    def test_adjudicate_deep_becomes_homref(self):
        out, stats = self._run(depth=35.0)
        self.assertEqual(out[0], "chr1\t100\tA\tG\t0/1\t0/0\t1/1")
        self.assertEqual(stats["cells_total"], 9)
        self.assertEqual(stats["dotdot_total"], 4)
        self.assertEqual(stats["filled_00"], 4)

    def test_adjudicate_shallow_keeps_nocall(self):
        out, stats = self._run(depth=5.0)
        self.assertEqual(out[0], "chr1\t100\tA\tG\t0/1\t./.\t1/1")
        self.assertEqual(stats["filled_00"], 0)
        self.assertEqual(stats["kept_dotdot"], 4)
        self.assertEqual(out[1], "chr1\t200\tC\tT\t./.\t0/0\t./.")


VCF_BODY = ("##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/1:40:20,20\n"
            "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP:AD\t1/1:30:0,30\n")


class TestRebuild(unittest.TestCase):
    """★ GT 列序踩坑：重建必须按 (CHROM,POS) 显式映射取 GT（防样本列错位），
    且只替换 GT、其余 FORMAT 字段原样保留"""

    def test_rebuild_replaces_gt_only(self):
        with tempfile.TemporaryDirectory() as td:
            vcf = os.path.join(td, "in.vcf")
            tsv = os.path.join(td, "gt.tsv")
            out = os.path.join(td, "out.vcf")
            open(vcf, "w").write(VCF_BODY)
            open(tsv, "w").write("chr1\t100\t0/0\n")       # 只改第一个位点
            st = mbc.rebuild_sample_vcf(vcf, tsv, out, "S1")
            self.assertEqual(st["records"], 2)
            self.assertEqual(st["gt_changed"], 1)
            lines = [l for l in open(out).read().splitlines() if not l.startswith("#")]
            self.assertEqual(lines[0].split("\t")[9], "0/0:40:20,20")   # DP/AD 保留
            self.assertEqual(lines[1].split("\t")[9], "1/1:30:0,30")    # 未映射不动

    def test_rebuild_column_map_explicit(self):
        """裁决矩阵列 → gt.tsv 必须按样本名取列（曾因线程池顺序错位互换两样本 GT）"""
        adj = ("chr1\t100\tA\tG\t0/0\t1/1\n")   # 列4=S1(0/0)、列5=S2(1/1)
        with tempfile.TemporaryDirectory() as td:
            adj_p = os.path.join(td, "adj.tsv")
            open(adj_p, "w").write(adj)
            for sm, want_col in (("S1", "0/0"), ("S2", "1/1")):
                gt_tsv = os.path.join(td, f"{sm}.gt.tsv")
                with open(adj_p) as fa, open(gt_tsv, "w") as fo:
                    idx = ["S1", "S2"].index(sm) + 4
                    for line in fa:
                        p = line.rstrip("\n").split("\t")
                        fo.write("\t".join([p[0], p[1], p[idx]]) + "\n")
                self.assertEqual(open(gt_tsv).read().strip(), f"chr1\t100\t{want_col}")


if __name__ == "__main__":
    unittest.main()
