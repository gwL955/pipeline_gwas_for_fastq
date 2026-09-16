#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""杂项单测：日志 tee 自动落盘（Logger(None) 控制台模式）、目录命名
（results/<批次>_<执行日期>）、交付导出（md5sum/MANIFEST/幂等）"""

import io
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config                                      # noqa: E402
from logger import Logger, capture_stdio           # noqa: E402


class TestLoggerTee(unittest.TestCase):

    def test_capture_stdio_mirrors_to_file(self):
        """自动日志：无 shell 重定向时 stdout 内容自动进入运行日志"""
        with tempfile.TemporaryDirectory() as td:
            logp = os.path.join(td, "run.log")
            old_out, old_err = sys.stdout, sys.stderr
            capture_stdio(logp)
            try:
                print("TEE_LINE_MARK")
            finally:
                sys.stdout, sys.stderr = old_out, old_err
            self.assertIn("TEE_LINE_MARK", open(logp, encoding="utf-8").read())

    def test_capture_stdio_buffered_then_flush(self):
        """缓冲模式（RUN-24）：capture_stdio(None) 先内存缓冲，目标确定后
        （tee_set_log_path）缓冲内容一并落盘、后续持续镜像——运行日志归属
        批次目录，但目标目录批次发现后才知道"""
        from logger import tee_set_log_path
        with tempfile.TemporaryDirectory() as td:
            logp = os.path.join(td, "batch_20xx", "logs", "run.log")
            old_out, old_err = sys.stdout, sys.stderr
            capture_stdio(None)
            try:
                print("EARLY_LINE")          # 目标确定前（资源计划阶段）
                tee_set_log_path(logp)
                print("LATER_LINE")          # 目标确定后
            finally:
                sys.stdout, sys.stderr = old_out, old_err
            body = open(logp, encoding="utf-8").read()
            self.assertIn("EARLY_LINE", body)   # 缓冲被 flush
            self.assertIn("LATER_LINE", body)
        # 从不设路径（dry-run/早退）→ 零落盘
        with tempfile.TemporaryDirectory() as td:
            before = os.listdir(td)
            old_out, old_err = sys.stdout, sys.stderr
            capture_stdio(None)
            try:
                print("DISCARDED_LINE")
            finally:
                sys.stdout, sys.stderr = old_out, old_err
            self.assertEqual(os.listdir(td), before)

    def test_results_root_only_batch_dirs(self):
        """目录规约（RUN-24）：results/ 根下只允许 批次_日期 目录——
        全局 run_summary.json 与 logs/ 不得再产生（dry-run 端到端验证）"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "DEMO")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(bdir, f"S1_L001_{r}_001.fastq.gz"), "wb") as f:
                    f.write(b"x")
            env = {**os.environ,
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low",
                 # ↑ low 档按 16 线程/20G 口径规划：CI runner 内存 ~16G，
                 #   auto 档探测会因单样本峰值 19.8G 触发快速失败（dry-run
                 #   不执行命令，规划口径不影响断言语义）
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 0, r.stderr[-500:])
            res = os.path.join(td, "results")
            names = os.listdir(res) if os.path.isdir(res) else []
            for n in names:                        # 只允许目录，且命名批次_日期
                self.assertTrue(os.path.isdir(os.path.join(res, n)), n)
                self.assertRegex(n, r"^DEMO_\d{8}$")
            self.assertNotIn("logs", names)        # 无全局日志目录
            self.assertNotIn("run_summary.json", names)   # 无全局散文件
    def test_full_flow_dry_run_zero_writes(self):
        """★ 全流程 dry-run（有效样本走到 Step 5/6）必须 success 且零落盘。
        gvcf.list 写入与 fastq_merged 目录创建曾无 dry-run 守卫——同日实跑过
        目录已存在而被掩盖，跨日新日期目录直接 FileNotFoundError（RUN-27 修复）"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "260422")
            os.makedirs(bdir)
            for lane in ("L001", "L002"):
                for r in ("R1", "R2"):
                    with open(os.path.join(
                            bdir, f"NA12878_S1_{lane}_{r}_001.fastq.gz"), "wb") as f:
                        f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low",
                 # ↑ low 档按 16 线程/20G 口径规划：CI runner 内存 ~16G，
                 #   auto 档探测会因单样本峰值 19.8G 触发快速失败（dry-run
                 #   不执行命令，规划口径不影响断言语义）
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 0, r.stdout[-800:])
            self.assertIn("批次 260422: success", r.stdout)
            res = os.path.join(td, "results")
            self.assertFalse(os.path.exists(res),
                             "dry-run 不得在 results/ 落任何目录/文件")

    def test_console_only_logger_no_file(self):
        lg = Logger(None, prefix="X")
        lg.info("不落盘")          # 不应抛异常
        self.assertEqual(lg.tail(), "")
        lg.close()

    def test_file_logger_tail(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "b.log")
            lg = Logger(p, prefix="B", quiet=True)
            lg.error("boom")
            self.assertIn("boom", lg.tail(5))
            self.assertIn("ERROR", lg.tail(5))
            lg.close()


class TestDirNaming(unittest.TestCase):

    def test_batch_result_dir_with_date(self):
        p = config.batch_result_dir("260422", "20260914")
        self.assertTrue(p.endswith(os.path.join("results", "260422_20260914")))

    def test_batch_result_dir_without_date(self):
        self.assertTrue(config.batch_result_dir("260422").endswith("results/260422")
                        or config.batch_result_dir("260422").endswith("results\\260422"))

    def test_batch_delivery_dir(self):
        p = config.batch_delivery_dir("260422", "20260914")
        self.assertTrue(p.endswith(os.path.join("Output", "260422_20260914")))

    def test_run_date_fixed_at_startup(self):
        """跨 0 点防回归锚（DEC-01）：执行日期必须取自启动时间戳 ts 一次性派生
        并传参到各批次（run_pipeline.py: run_date = ts.split("_")[0]）；
        禁止任何人把 batch_result_dir 的调用改为运行中现取 datetime.now()
        ——否则 23:50 启动的运行跨午夜后会把产物写进两个日期目录"""
        import inspect
        import run_pipeline
        main_src = inspect.getsource(run_pipeline.main)
        self.assertIn('run_date = ts.split("_")[0]', main_src)
        # process_batch 只接受外部传入的 run_date，函数体内不再现取日期
        pb_src = inspect.getsource(run_pipeline.process_batch)
        self.assertNotIn("strftime(\"%Y%m%d\")", pb_src.replace("'", '"'))

    def test_multiqc_ordering_after_all_steps_before_delivery(self):
        """★ MultiQC 时序锚（v2.3.0 用户口径）：其他全部流程（含矩阵裁决与
        每样本 VCF 重建）结束、QC 文件生成完全后才生成最终汇总报告，随后才
        交付导出——防止回归成"裁决前跑 MultiQC"或"交付后才跑"（报告缺最后产物
        或进不了交付目录）"""
        import inspect
        import run_pipeline
        src = inspect.getsource(run_pipeline.process_batch)
        i_rebuild = src.index('bdata["_adj_vcfs"] = adj_vcfs')
        i_multiqc = src.index("run_multiqc")
        i_export = src.index("export_delivery(")
        self.assertLess(i_rebuild, i_multiqc,
                        "MultiQC 必须在矩阵裁决与每样本 VCF 重建之后")
        self.assertLess(i_multiqc, i_export,
                        "MultiQC 报告生成后才交付导出（报告随交付拷贝）")


class TestHostPythonGuard(unittest.TestCase):
    """★ 嵌套容器防回归（RUN-22/23）：本机 `python` 是容器别名，容器内运行
    曾致全部子进程 `singularity: not found`——必须在启动即失败并给出指引"""

    def test_guard_blocks_inside_container(self):
        from run_pipeline import guard_host_python
        for marker in ({"SINGULARITY_NAME": "mamba_2.6.2.sif"},
                       {"SINGULARITY_CONTAINER": "/x.sif"}):
            with mock.patch.dict(os.environ, marker):
                with self.assertRaises(SystemExit):
                    guard_host_python()

    def test_guard_passes_on_host(self):
        from run_pipeline import guard_host_python
        host_env = {k: v for k, v in os.environ.items()
                    if not k.startswith("SINGULARITY_")}
        with mock.patch.dict(os.environ, host_env, clear=True):
            self.assertIsNone(guard_host_python())   # 宿主环境直接放行


class _FakeRunner:
    dry_run = True          # 跳过 bcftools 计数（无容器环境）

    def cpath(self, p):
        return p


class _FakeLog:
    def __getattr__(self, name):
        return lambda *a, **k: None


class TestExportDelivery(unittest.TestCase):

    def _src(self, td):
        import hashlib
        os.makedirs(td, exist_ok=True)
        files = {}
        for sm, payload in (("NA12878", b"VCF-A"), ("NA18544", b"VCF-B")):
            v = os.path.join(td, f"{sm}.PASS.adjudicated.vcf.gz")
            open(v, "wb").write(payload)
            open(v + ".tbi", "wb").write(b"TBI")
            files[sm] = v
        return files

    def test_export_delivery_layout_and_md5(self):
        from run_pipeline import export_delivery
        with tempfile.TemporaryDirectory() as src, \
                tempfile.TemporaryDirectory() as dst:
            files = self._src(src)
            mq = os.path.join(src, "GWAS-Panel_multiqc_report.html")
            open(mq, "wb").write(b"<html>MQ</html>")
            with mock.patch.object(config, "DELIVERY_DIR", dst):
                out = export_delivery(_FakeRunner(), "260422_20260914", files,
                                      _FakeLog(), batch_note="｜批次 260422",
                                      extra_files=[mq])
            self.assertTrue(out.endswith("260422_20260914"))
            names = os.listdir(out)
            for f in ("NA12878.PASS.adjudicated.vcf.gz",
                      "NA12878.PASS.adjudicated.vcf.gz.tbi",
                      "GWAS-Panel_multiqc_report.html",
                      "md5sum.txt", "MANIFEST.tsv", "README.md"):
                self.assertIn(f, names)
            import hashlib
            expect = hashlib.md5(b"VCF-A").hexdigest()
            md5s = open(os.path.join(out, "md5sum.txt")).read()
            self.assertIn(f"{expect}  NA12878.PASS.adjudicated.vcf.gz", md5s)
            self.assertIn(hashlib.md5(b"<html>MQ</html>").hexdigest()
                          + "  GWAS-Panel_multiqc_report.html", md5s)
            man = open(os.path.join(out, "MANIFEST.tsv")).read()
            self.assertIn("sample\tfile\tsize_bytes\tmd5\trecords\tsource", man)
            self.assertIn(expect, man)
            self.assertIn("multiqc", man)          # 附加文件登记（sample 列=multiqc）
            rdm = open(os.path.join(out, "README.md")).read()
            self.assertIn("DP≥" + str(config.DP_MIN), rdm)
            self.assertIn("MultiQC 汇总报告", rdm)  # 交付说明含 MultiQC 行

    def test_export_delivery_without_extra_files(self):
        """不传 extra_files 时行为与旧版一致（无 MultiQC 行/登记）"""
        from run_pipeline import export_delivery
        with tempfile.TemporaryDirectory() as src, \
                tempfile.TemporaryDirectory() as dst:
            files = self._src(src)
            with mock.patch.object(config, "DELIVERY_DIR", dst):
                out = export_delivery(_FakeRunner(), "B_20260914", files, _FakeLog())
            names = os.listdir(out)
            self.assertNotIn("GWAS-Panel_multiqc_report.html", names)
            self.assertNotIn("MultiQC 汇总报告",
                             open(os.path.join(out, "README.md")).read())

    def test_export_delivery_idempotent(self):
        from run_pipeline import export_delivery
        with tempfile.TemporaryDirectory() as src, \
                tempfile.TemporaryDirectory() as dst:
            files = self._src(src)
            mq = os.path.join(src, "GWAS-Panel_multiqc_report.html")
            open(mq, "wb").write(b"<html>MQ</html>")
            with mock.patch.object(config, "DELIVERY_DIR", dst):
                out1 = export_delivery(_FakeRunner(), "B_20260914", files, _FakeLog(),
                                       extra_files=[mq])
                mt1 = os.path.getmtime(os.path.join(out1, "NA12878.PASS.adjudicated.vcf.gz"))
                mtq1 = os.path.getmtime(os.path.join(out1, os.path.basename(mq)))
                out2 = export_delivery(_FakeRunner(), "B_20260914", files, _FakeLog(),
                                       extra_files=[mq])
                mt2 = os.path.getmtime(os.path.join(out2, "NA12878.PASS.adjudicated.vcf.gz"))
                mtq2 = os.path.getmtime(os.path.join(out2, os.path.basename(mq)))
            self.assertEqual(out1, out2)
            self.assertEqual(mt1, mt2)      # 源未更新不重拷（mtime 保留）
            self.assertEqual(mtq1, mtq2)    # 附加文件（MultiQC）同样幂等


class TestDeliveryIndex(unittest.TestCase):

    def test_index_accumulates_history(self):
        """★ INDEX.md 累积锚（RUN-29）：合并 Output/ 全部历史交付目录 ∪ 本次运行——
        曾只写本次运行批次，跨日新运行会把历史交付从索引中挤掉
        （260422_20260914 在 260422/260720_20260916 运行后从索引消失）"""
        from run_pipeline import write_delivery_index
        with tempfile.TemporaryDirectory() as td:
            for b in ("260422_20260914", "260720_20260916"):   # 历史交付（已在磁盘）
                d = os.path.join(td, b)
                os.makedirs(d)
                open(os.path.join(d, "x.PASS.adjudicated.vcf.gz"), "wb").write(b"v")
            new = os.path.join(td, "260529_20260916")           # 本次运行交付
            os.makedirs(new)
            open(os.path.join(new, "y.PASS.adjudicated.vcf.gz"), "wb").write(b"v")
            with mock.patch.object(config, "DELIVERY_DIR", td):
                idx = write_delivery_index({new})
            body = open(idx, encoding="utf-8").read()
            for b in ("260422_20260914", "260529_20260916", "260720_20260916"):
                self.assertIn(b, body)          # 历史 + 本次全部在索引中


if __name__ == "__main__":
    unittest.main()
