#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""各软件模块的解析器单测（fixture 取自真实输出格式；对应踩坑：
fastqc 状态为小写、markdup/hsmetrics 按表头名解析、bcftools stats SN 行带文件 ID 列、
mosdepth summary 为 6 列且靶区深度应取 total_region 行、norm 统计行 7 字段动态表头）"""

import io
import os
import sys
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules import fastqc as mfastqc                # noqa: E402
from modules import fastp as mfastp                  # noqa: E402
from modules import samtools as msam                 # noqa: E402
from modules import gatk as mgatk                    # noqa: E402
from modules import bcftools as mbc                  # noqa: E402
from modules import mosdepth as mmos                 # noqa: E402


class _L:
    def __getattr__(self, name):
        return lambda *a, **k: None


def _zip_with(data_text, name="SM_R1_fastqc"):
    import tempfile
    fd, path = tempfile.mkstemp(suffix=".zip")
    os.close(fd)
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{name}/fastqc_data.txt", data_text)
    return path


FASTQC_DATA = """##FastQC 0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\tSM_R1.fastq.gz
Total Sequences\t2387762
Sequence length\t150
>>Per base sequence quality\twarn
>>Adapter Content\tfail
#Position\tAdapter
>>END_MODULE
"""


class TestFastqc(unittest.TestCase):

    def test_parse_statuses_lowercase(self):
        """踩坑回归：fastqc 模块状态是小写 pass/warn/fail，曾按大写 PASS 比较误报"""
        z = _zip_with(FASTQC_DATA)
        try:
            modules, basic = mfastqc.parse_fastqc_zip(z)
            self.assertEqual(modules["Basic Statistics"], "pass")
            self.assertEqual(modules["Adapter Content"], "fail")
            self.assertEqual(basic["Total Sequences"], "2387762")
        finally:
            os.unlink(z)

    def test_adapter_cleared_case_insensitive(self):
        z_raw = _zip_with(FASTQC_DATA)
        z_trim = _zip_with(FASTQC_DATA.replace(">>Adapter Content\tfail",
                                               ">>Adapter Content\tpass"))
        try:
            mfastqc.check_adapter_cleared(z_raw, z_trim, _L())   # 不应抛异常
        finally:
            os.unlink(z_raw), os.unlink(z_trim)


class TestFastp(unittest.TestCase):

    def test_parse_json(self):
        import tempfile
        # 结构按 fastp 1.3.6 真实产物校准（RUN-25）：after_filtering 只有整体
        # q30_rate（小数）；R1/R2 分列在顶层 read{1,2}_after_filtering 且无
        # q30_rate——曾凭想象造 q30_rate_r1 键，测试全绿但真实数据恒 None
        j = {"summary": {"before_filtering": {"total_reads": 1000, "total_bases": 150000},
                         "after_filtering": {"total_reads": 920, "total_bases": 138000,
                                             "q30_rate": 0.85764}},
             "read1_after_filtering": {"total_bases": 69000, "q30_bases": 60720},
             "read2_after_filtering": {"total_bases": 69000, "q30_bases": 58650}}
        fd, p = tempfile.mkstemp(suffix=".json")
        with os.fdopen(fd, "w") as f:
            import json
            json.dump(j, f)
        met = mfastp.parse_json(p)
        os.unlink(p)
        self.assertEqual(met["retention_pct"], 92.0)
        self.assertEqual(met["q30_pct"], 85.76)        # 整体：小数→百分数
        self.assertEqual(met["q30_pct_r1"], 88.0)      # bases 自算
        self.assertEqual(met["q30_pct_r2"], 85.0)


MARKDUP_METRICS = """## htsjdk.samtools.metrics.StringHeader
# MarkDuplicates INPUT=...
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
Unknown Library\t0\t2193438\t0\t0\t0\t44236\t10920\t0.020183\t54890174
"""


class TestGatkParsers(unittest.TestCase):

    def test_markdup_by_header_name(self):
        import tempfile
        fd, p = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(MARKDUP_METRICS)
        d = mgatk.parse_markdup_metrics(p)
        os.unlink(p)
        self.assertEqual(d["READ_PAIRS_EXAMINED"], 2193438)
        self.assertAlmostEqual(d["PERCENT_DUPLICATION"], 0.020183)
        self.assertEqual(d["ESTIMATED_LIBRARY_SIZE"], 54890174)

    def test_hsmetrics(self):
        import tempfile
        cols = ["BAIT_SET", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE",
                "PCT_TARGET_BASES_20X", "PCT_SELECTED_BASES", "ON_BAIT_BASES",
                "ON_TARGET_BASES", "PF_UQ_BASES_ALIGNED"]
        vals = ["targets", "93.519497", "93.0", "0.953103", "0.395441",
                "1000", "600", "4000"]
        text = "## METRICS CLASS\tpicard.analysis.directed.HsMetrics\n" \
            + "\t".join(cols) + "\n" + "\t".join(vals) + "\n"
        fd, p = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(text)
        d = mgatk.parse_hsmetrics(p)
        os.unlink(p)
        self.assertEqual(d["MEAN_TARGET_COVERAGE"], 93.519497)
        self.assertEqual(d["PCT_TARGET_BASES_20X"], 0.953103)
        self.assertEqual(d["PCT_SELECTED_BASES"], 0.395441)   # 捕获效率口径（DEC-29）
        self.assertEqual(d["ON_BAIT_PCT"], 25.0)    # 1000/4000×100
        self.assertEqual(d["ON_TARGET_PCT"], 15.0)


FLAGSTAT = """6813859 + 0 in total (QC-passed reads + QC-failed reads)
6733928 + 0 primary
0 + 0 duplicates
6813270 + 0 mapped (99.99% : N/A)
6733928 + 0 paired in sequencing
6450400 + 0 properly paired (95.79% : N/A)
489 + 0 singletons (0.01% : N/A)
242774 + 0 with mate mapped to a different chr
217338 + 0 with mate mapped to a different chr (mapQ>=5)
"""

STATS = """# This file was produced by samtools stats
SN\traw total sequences:\t6813859
SN\treads mapped:\t6813270
SN\tinsert size average:\t290.5
SN\taverage quality:\t36.1
"""


class TestSamtools(unittest.TestCase):

    def test_flagstat(self):
        import tempfile
        fd, p = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(FLAGSTAT)
        d = msam.parse_flagstat(p)
        os.unlink(p)
        self.assertEqual(d["total"], 6813859)
        self.assertEqual(d["mapped"], 6813270)
        self.assertEqual(d["mapped_pct"], 99.99)
        self.assertEqual(d["pp_pct"], 95.79)
        self.assertEqual(d["sgl_pct"], 0.01)
        self.assertEqual(d["diff_chr"], 242774)

    def test_stats(self):
        import tempfile
        fd, p = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(STATS)
        d = msam.parse_stats(p)
        os.unlink(p)
        self.assertEqual(d["raw total sequences"], 6813859.0)
        self.assertEqual(d["insert size average"], 290.5)

    def test_flagstat_identical(self):
        import tempfile
        fd, p1 = tempfile.mkstemp()
        os.close(fd)
        fd, p2 = tempfile.mkstemp()
        os.close(fd)
        open(p1, "w").write(FLAGSTAT)
        open(p2, "w").write(FLAGSTAT.replace("95.79", "95.80"))
        self.assertTrue(msam.flagstat_identical(p1, p1))
        self.assertFalse(msam.flagstat_identical(p1, p2))
        os.unlink(p1), os.unlink(p2)


BCF_STATS = """# This file was produced by bcftools stats
SN\t0\tnumber of samples:\t2
SN\t0\tnumber of records:\t11865
SN\t0\tnumber of no-ALTs:\t0
SN\t0\tnumber of SNPs:\t10279
SN\t0\tnumber of MNPs:\t0
SN\t0\tnumber of indels:\t1599
SN\t0\tnumber of others:\t0
SN\t0\tnumber of multiallelic sites:\t210
TSTV\t0\t6882\t3435\t2.00\t6865\t3395\t2.02
"""

NORM_LINE = ("Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:"
             "\t11865/210/0/117/0/0/0\n")

MOSDEPTH_SUMMARY = """chrom\tlength\tbases\tmean\tmin\tmax
chr1\t248956422\t59460088\t0.24\t0\t489
chr1_region\t104196\t16813882\t161.37\t0\t453
total\t3099519240\t615070121\t0.20\t0\t915
total_region\t1006786\t156384830\t155.33\t0\t662
"""


class TestBcftoolsMosdepth(unittest.TestCase):

    def test_bcftools_stats_with_file_id_column(self):
        """踩坑回归：SN 行带文件 ID 列（SN 0 number of SNPs: N），按列位置解析曾得 None"""
        import tempfile
        fd, p = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(BCF_STATS)
        d = mbc.parse_stats(p)
        os.unlink(p)
        self.assertEqual(d["records"], 11865)
        self.assertEqual(d["snps"], 10279)
        self.assertEqual(d["indels"], 1599)
        self.assertEqual(d["ts"], 6882)
        self.assertEqual(d["tv"], 3435)
        self.assertEqual(d["titv"], 2.00)

    def test_norm_stats_dynamic_header(self):
        """踩坑回归：bcftools 1.24 的 norm 统计行是 7 字段（含 joined/…_removed），
        按旧 4 字段正则解析曾得空 dict"""
        d = mbc._parse_norm_stats(NORM_LINE)
        self.assertEqual(d["total"], 11865)
        self.assertEqual(d["split"], 210)
        self.assertEqual(d["realigned"], 117)
        self.assertEqual(d["skipped"], 0)
        self.assertEqual(mbc._parse_norm_stats("no match here"), {})

    def test_mosdepth_six_columns_total_region(self):
        """踩坑回归：summary 是 6 列（chrom length bases mean min max），mean 在第 4 列；
        且 --by 模式下 total 行被摊薄，靶区深度必须取 total_region 行"""
        import tempfile
        fd, p = tempfile.mkstemp(suffix=".mosdepth.summary.txt")
        os.close(fd)
        prefix = p[:-len(".mosdepth.summary.txt")]
        with open(p, "w") as f:
            f.write(MOSDEPTH_SUMMARY)
        d = mmos.parse_summary(prefix)
        os.unlink(p)
        self.assertEqual(d["mean"], 155.33)          # total_region 而非 total 的 0.20
        self.assertEqual(d["chroms"]["total"], 0.20)

    def test_region_depth_binary_search(self):
        regions = [("chr1", 100, 200, 30.0), ("chr1", 200, 300, 25.0),
                   ("chr2", 1, 50, 12.0)]
        regions.sort(key=lambda r: (r[0], r[1]))
        self.assertEqual(mmos.region_depth(regions, "chr1", 150), 30.0)
        self.assertEqual(mmos.region_depth(regions, "chr1", 250), 25.0)
        self.assertEqual(mmos.region_depth(regions, "chr2", 10), 12.0)
        self.assertIsNone(mmos.region_depth(regions, "chr1", 350))


class TestRecalObservations(unittest.TestCase):
    """RecalTable1 M 事件观测数解析（RUN-34）：按表头名取列、只加 M 行、
    下一 # 表即止；known-sites 覆盖崩坏时该值骤降 → P2 信号"""

    FIXTURE = (
        "#:GATKReport.v1.1:5\n"
        "#:GATKTable:2:17:%s:%s:;\n"
        "#:GATKTable:Arguments:...\n"
        "covariate  ReadGroupCovariate\n"
        "#:GATKTable:6:3:%s:%s:%.4g:%.4g:%.4g:%.4g:;\n"
        "#:GATKTable:RecalTable1:\n"
        "ReadGroup  QualityScore  EventType  EmpiricalQuality  Observations  Errors\n"
        "NA12878  14  M  13.0000  45782114  2230327.00\n"
        "NA12878  14  I   9.0000       1234       45.00\n"
        "NA12878  21  M  21.0000   6627250    52626.00\n"
        "#:GATKTable:RecalTable2:\n"
        "ReadGroup  QualityScore  EventType  Observations\n"
        "NA12878  40  M  999999999\n")   # 后续表不得计入

    def test_sums_m_event_observations(self):
        from modules import gatk as mg
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".table", delete=False) as f:
            f.write(self.FIXTURE)
            path = f.name
        try:
            self.assertEqual(mg.parse_recal_observations(path),
                             45782114 + 6627250)
        finally:
            os.unlink(path)

    def test_missing_file_returns_none(self):
        from modules import gatk as mg
        self.assertIsNone(mg.parse_recal_observations("/nonexistent/x.table"))


if __name__ == "__main__":
    unittest.main()
