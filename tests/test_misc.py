#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""杂项单测：日志 tee 自动落盘（Logger(None) 控制台模式）、目录命名
（results/<批次>_<执行日期>）、交付导出（md5sum/MANIFEST/幂等）"""

import io
import os
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config                                      # noqa: E402
from logger import Logger, capture_stdio           # noqa: E402


def _dep_env(td):
    """E2E dry-run 的哑依赖环境：开跑前 P0 依赖检查（reference/SIF 缺失即
    RuntimeError 中止——实跑快速失败语义）要求非空文件存在；CI/干净克隆没有
    真实数据，经 GWAS_REFERENCE_DIR/GWAS_SIF_DIR 重定向到临时哑树即可全流程
    dry-run（dry-run 不执行命令，哑文件内容无关）"""
    ref = os.path.join(td, "deps", "reference")
    sif = os.path.join(td, "deps", "singularity")
    os.makedirs(os.path.join(ref, "genome"), exist_ok=True)
    os.makedirs(sif, exist_ok=True)
    for n in ("genome.fa", "genome.fa.fai", "genome.dict"):
        open(os.path.join(ref, "genome", n), "wb").write(b"x")
    for n in ("targets.bed", "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
              "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"):
        open(os.path.join(ref, n), "wb").write(b"x")
    for v in config.SIF.values():
        open(os.path.join(sif, os.path.basename(v)), "wb").write(b"x")
    return {"GWAS_REFERENCE_DIR": ref, "GWAS_SIF_DIR": sif,
            "GWAS_DISK_MIN_FREE_GB": "0"}   # 磁盘 P0 现为阻断级（v2.9.0），测试环境豁免


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

    def test_tee_set_echo_silences_console_not_file(self):
        """★ 实跑外显关闭锚（DEC-09 v2.13.0/RUN-37）：tee_set_echo(False) 后
        控制台不再回显（nohup 后台不向 nohup.out 倾倒运行日志），
        文件镜像不受影响（含关闭前的缓冲前缀）"""
        from logger import tee_set_log_path, tee_set_echo
        with tempfile.TemporaryDirectory() as td:
            logp = os.path.join(td, "batch_x", "logs", "run.log")
            con = io.StringIO()
            old_out, old_err = sys.stdout, sys.stderr
            capture_stdio(None)
            sys.stdout._stream = con          # 控制台端换成可断言的内存流
            try:
                print("BEFORE_ECHO_OFF")      # 外显期：控制台+缓冲
                tee_set_log_path(logp)
                tee_set_echo(False)
                print("AFTER_ECHO_OFF")       # 静默期：只落文件
            finally:
                sys.stdout, sys.stderr = old_out, old_err
            body = open(logp, encoding="utf-8").read()
            self.assertIn("BEFORE_ECHO_OFF", body)   # 缓冲前缀完整落盘
            self.assertIn("AFTER_ECHO_OFF", body)    # 文件镜像不受 echo 开关影响
            self.assertIn("BEFORE_ECHO_OFF", con.getvalue())
            self.assertNotIn("AFTER_ECHO_OFF", con.getvalue())   # 控制台已静默

    def test_real_run_silent_after_log_path_line(self):
        """★ 实跑端到端静默锚（RUN-37）：真跑（非 dry-run）控制台止于"运行日志:
        <路径>"提示行，批次处理日志只进 run_<ts>.log；dry-run 行为不变（见相邻
        dry-run 用例对 stdout 的断言）"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "in", "EMPTY"))   # 空批次 → skipped，无需容器
            env = {**os.environ, "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--resource-profile", "low", "--notify", "off",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr[-500:])
            self.assertIn("运行日志", r.stdout)              # 路径提示是控制台最后一截
            self.assertNotIn("批次 EMPTY 结束", r.stdout)     # 处理日志不再外显
            self.assertNotIn("skipped", r.stdout)
            # run_<ts>.log 含全量日志（含静默后的批次处理段）；目录名日期动态匹配
            res = os.path.join(td, "results")
            bdirs = [n for n in os.listdir(res)
                     if n.startswith("EMPTY_") and os.path.isdir(os.path.join(res, n))]
            self.assertTrue(bdirs, "批次结果目录必须生成")
            logs_dir = os.path.join(res, bdirs[0], "logs")
            run_logs = [n for n in os.listdir(logs_dir) if n.startswith("run_")]
            self.assertTrue(run_logs, "run_<ts>.log 必须落盘")
            body = open(os.path.join(logs_dir, run_logs[0]), encoding="utf-8").read()
            self.assertIn("批次 EMPTY", body)
            self.assertIn("skipped", body)

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
            env = {**os.environ, **_dep_env(td),
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
            env = {**os.environ, **_dep_env(td),
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
            self.assertEqual(r.returncode, 0, r.stdout[-800:])
            self.assertIn("批次 260422: success", r.stdout)
            res = os.path.join(td, "results")
            self.assertFalse(os.path.exists(res),
                             "dry-run 不得在 results/ 落任何目录/文件")

    def test_full_flow_dry_run_outsourced_flat(self):
        """★ 外送平铺布局端到端锚（RUN-36/DEC-23）：<样本>_R{1,2}.fastq.gz 直接放
        批次目录须走完全流程 success——此前两种识别器都不认，整批静默 skipped"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "260917")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(bdir, f"NA12878_{r}.fastq.gz"), "wb") as f:
                    f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 0, r.stdout[-800:])
            self.assertIn("批次 260917: success", r.stdout)
            self.assertIn("有效样本 1 个", r.stdout)
            self.assertFalse(os.path.exists(os.path.join(td, "results")))

    def test_bad_sample_name_rejected(self):
        """★ 样本名白名单 P0 阻断锚（DEC-19/RUN-32）：非常规字符（如 ;）→
        P0 级错误直接中断分析（退出码非零、批次 failed）——样本名进入 shell
        命令拼接与 @RG 头属注入面；v2.7.0 曾为剔除继续，v2.8.0 收紧为阻断"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "B")
            os.makedirs(bdir)
            for sm in ("bad;name_S1", "NA12878_S2"):
                for r in ("R1", "R2"):
                    with open(os.path.join(
                            bdir, f"{sm}_L001_{r}_001.fastq.gz"), "wb") as f:
                        f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 1, r.stdout[-500:])   # 批次 failed
            self.assertIn("开跑前-P0", r.stdout)                  # P0 异常项落日志
            self.assertIn("样本名含中文或非常规字符", r.stdout)   # RUN-38 起消息点名中文
            self.assertIn("bad;name", r.stdout)                   # 指名道姓
            self.assertNotIn("有效样本 1 个", r.stdout)            # 不再剔除继续

    def test_chinese_batch_name_p0(self):
        """★ 批次名中文 P0 阻断锚（RUN-38）：中文批次名进入容器 bind/命令路径
        在 singularity 内因 locale 报错——开跑前直接中断（退出码 1、指名批次名）"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "测试批次")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(
                        bdir, f"NA12878_S1_L001_{r}_001.fastq.gz"), "wb") as f:
                    f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low", "--notify", "off",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 1, r.stdout[-500:])   # 批次 failed
            self.assertIn("开跑前-P0", r.stdout)
            self.assertIn("批次名 测试批次 含中文或非常规字符", r.stdout)   # 指名道姓
            self.assertNotIn("Step 0 完成", r.stdout)             # 中断于合并之前

    def test_gitignore_blocks_private_bed(self):
        """★ 靶区 bed 防泄露锚（RUN-38）：panel 靶区设计属私密内容，
        .gitignore 必须拦截 *.bed 与派生 *.interval_list（改址经 .env GWAS_TARGETS_BED）"""
        gi = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), ".gitignore")
        text = open(gi, encoding="utf-8").read()
        self.assertIn("*.bed", text)
        self.assertIn("*.interval_list", text)

    def test_relative_input_respects_raw_data_override(self):
        """★ --input 相对路径解析锚（RUN-32 事故）：相对名优先在 RAW_DATA_DIR
        （含 GWAS_RAW_DATA 覆盖）下解析——曾直接回落 $WORK/0_raw_data，导致
        验证脚本带着覆盖变量却扫到真实批次发起实跑"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "raw", "0_raw_data", "DEMO"))
            env = {**os.environ,
                   "GWAS_RAW_DATA": os.path.join(td, "raw", "0_raw_data"),
                   "GWAS_RESULTS": os.path.join(td, "res")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low",   # CI 内存 <auto 档峰值
                 "--input", "0_raw_data"],
                env=env, capture_output=True, text=True, timeout=60)
            self.assertIn("待处理批次: ['DEMO']", r.stdout)

    def test_cli_version_flag(self):
        """--version：输出版本并退出码 0（与 DESIGN.md 同步，check_design 校验）"""
        import subprocess
        r = subprocess.run(
            [sys.executable,
             os.path.join(os.path.dirname(os.path.dirname(
                 os.path.abspath(__file__))), "run_pipeline.py"),
             "--version"], capture_output=True, text=True, timeout=30)
        self.assertEqual(r.returncode, 0)
        self.assertIn(config.PIPELINE_VERSION, r.stdout)

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
        # v2.7.0 起步骤实现拆分为 BatchCtx 方法，锚定 step6 方法源码
        src = inspect.getsource(run_pipeline.BatchCtx.step6_summary_delivery)
        i_rebuild = src.index('self.bdata["_adj_vcfs"] = adj_vcfs')
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
            self.assertIn(f"v{config.PIPELINE_VERSION}", rdm)   # 溯源：交付物带版本
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


class TestDiskGuard(unittest.TestCase):

    def test_disk_guard_between_steps_p0(self):
        """★ P0 Step 间磁盘复查锚（RUN-34）：低于阈值立即 RuntimeError 终止批次；
        dry-run 零落盘跳过检查"""
        import argparse
        import shutil as _shutil
        from run_pipeline import BatchCtx
        with tempfile.TemporaryDirectory() as td:
            args = argparse.Namespace(dry_run=True)
            ctx = BatchCtx("B", td, os.path.join(td, "w"), args, None, None)
            ctx.disk_guard("step0")            # dry-run：不检查不抛
            ctx.runner.dry_run = False         # 切换为实跑语义再验阈值
            with mock.patch.object(config, "DISK_MIN_FREE_GB", 1e9), \
                    mock.patch("shutil.disk_usage",
                               return_value=_shutil.disk_usage("/")):
                with self.assertRaises(RuntimeError):   # 阈值 1e9GB 必低于真实剩余
                    ctx.disk_guard("step1")
            with mock.patch.object(config, "DISK_MIN_FREE_GB", 0):
                ctx.disk_guard("step2")        # 阈值 0：正常放行


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


class TestIntervalListPrep(unittest.TestCase):
    """★ BedToIntervalList 参数锚（RUN-40，DEC-26）：--UNIQUE true 去重合并 +
    --DROP_MISSING_CONTIGS true 丢字典外 contig——新 panel bed 按 GRCh38 完整版
    （含 ALT contig）制定，genome.dict（194 序列无 ALT）遇 `chr22_KI270879v1_alt`
    即 PicardException 中断；重叠/相邻探针区间须合并为唯一区间（靶区碱基唯一口径）"""

    def _prep(self, td, sorted_exists=False, interval_exists=False):
        from modules import gatk as mgatk
        bed = os.path.join(td, "targets.bed")
        open(bed, "wb").write(b"chr1\t100\t200\n")
        sorted_bed = os.path.join(td, "targets.sorted.bed")
        ilist = os.path.join(td, "targets.sorted.interval_list")
        if sorted_exists:
            open(sorted_bed, "wb").write(b"chr1\t100\t200\n")
        if interval_exists:
            open(ilist, "wb").write(b"@HD\tVN:1.6\n")
        cmds = []

        class R:
            dry_run = True      # prep 的新鲜度重建 post-check 只在实跑校验产物

            def run(self, cmd, logger=None, outputs=(), timeout=None, capture=False):
                cmds.append(cmd)
                return 0

            def tool(self, sif_key, args, binds=None):
                return f"singularity:{sif_key} {args}"

            def cpath(self, path):
                return path

        with mock.patch.object(config, "TARGETS_BED", bed), \
                mock.patch.object(config, "TARGETS_SORTED_BED", sorted_bed), \
                mock.patch.object(config, "TARGETS_INTERVAL_LIST", ilist), \
                mock.patch.object(config, "GENOME_DICT", os.path.join(td, "genome.dict")):
            ok = mgatk.prep_interval_list(R(), _FakeLog())
        return ok, cmds

    def test_bedtointervallist_flags(self):
        """两参数必须同时在命令行上：漏 --DROP_MISSING_CONTIGS 遇 ALT contig 必中断；
        漏 --UNIQUE 重叠探针虚增靶区碱基（HsMetrics 口径失真）。
        awk 侧必须带 genome.dict contig 白名单（DEC-28）：漏过滤则 bed 中字典外
        contig（ALT）进 sorted.bed，坑留给 bcftools/mosdepth"""
        with tempfile.TemporaryDirectory() as td:
            ok, cmds = self._prep(td)
            self.assertTrue(ok)
            self.assertEqual(len(cmds), 2)          # awk 排序去重 + BedToIntervalList
            self.assertIn("sort -k1,1V -k2,2n -u", cmds[0])   # bed 级去完全重复行
            self.assertIn("NR==FNR", cmds[0])         # 字典白名单装载（DEC-28）
            self.assertIn("$1 in c", cmds[0])         # bed 行 contig 必须在白名单内
            self.assertIn("genome.dict", cmds[0])     # 字典作为 awk 首输入
            gatk_cmd = cmds[1]
            self.assertIn("BedToIntervalList", gatk_cmd)
            self.assertIn("--UNIQUE true", gatk_cmd)
            self.assertIn("--DROP_MISSING_CONTIGS true", gatk_cmd)
            self.assertIn("-SD", gatk_cmd)

    def test_interval_list_idempotent_skip(self):
        """sorted.bed 与 interval_list 均已存在（非空且不旧）→ 零命令直接放行（REQ-04 幂等）"""
        with tempfile.TemporaryDirectory() as td:
            ok, cmds = self._prep(td, sorted_exists=True, interval_exists=True)
            self.assertTrue(ok)
            self.assertEqual(cmds, [])

    def test_derived_targets_stale_rebuild(self):
        """★ 派生靶区新鲜度锚（DEC-28/RUN-42）：源 bed 更新（mtime 更新）后，
        已存在且非空的 sorted.bed 必须重建——曾只看非空即跳过，换 bed 后旧派生
        文件（含 ALT contig 行）一路带进下游；interval_list 同理随 sorted.bed 失效"""
        with tempfile.TemporaryDirectory() as td:
            ok, cmds = self._prep(td, sorted_exists=True, interval_exists=True)
            self.assertEqual(cmds, [])                       # 新鲜 → 全跳过
            os.utime(os.path.join(td, "targets.bed"),
                     (time.time() + 5, time.time() + 5))     # 源 bed 更新
            ok, cmds = self._prep(td)
            self.assertTrue(ok)
            self.assertEqual(len(cmds), 1)                   # 只重建 sorted.bed
            self.assertTrue(cmds[0].startswith("awk"))
            os.utime(os.path.join(td, "targets.sorted.bed"),
                     (time.time() + 5, time.time() + 5))     # sorted.bed 更新
            ok, cmds = self._prep(td)
            self.assertTrue(ok)
            self.assertEqual(len(cmds), 1)                   # 只重建 interval_list
            self.assertIn("BedToIntervalList", cmds[0])

    def test_haplotypecaller_uses_interval_list(self):
        """★ HC -L 口径锚（DEC-28/RUN-42）：-L 必须用 interval_list（字典口径+
        DEC-26 唯一合并），不得用 sorted.bed——GATK 引擎对 -L 区间 contig 严格
        校验字典，bed 含字典外 contig 即 USER ERROR（实跑三样本 HC 全灭）"""
        from modules import gatk as mgatk
        cmds = []

        class R:
            def run(self, cmd, logger=None, outputs=(), timeout=None, capture=False):
                cmds.append(cmd)
                return 0

            def tool(self, sif_key, args, binds=None):
                return f"singularity:{sif_key} {args}"

            def cpath(self, path):
                return path

        with mock.patch.object(config, "TARGETS_INTERVAL_LIST",
                               "/x/targets.sorted.interval_list"), \
                mock.patch.object(config, "TARGETS_SORTED_BED",
                                  "/x/targets.sorted.bed"):
            ok = mgatk.haplotypecaller(R(), "/x/a.bam", "/x/a.g.vcf.gz",
                                       "2g", 2, _FakeLog())
        self.assertTrue(ok)
        self.assertIn("-L /x/targets.sorted.interval_list", cmds[0])
        self.assertNotIn("targets.sorted.bed", cmds[0])


class TestBatchExceptionPath(unittest.TestCase):

    def test_batch_exception_path_no_unbound_notify(self):
        """★ 失败路径锚（RUN-42）：step 中途异常时 except 分支不得再抛
        UnboundLocalError——曾引用仅在成功汇总段才赋值的局部 notify_on，导致
        钉钉失败通知发不出、且裸 traceback 炸穿 main()（nohup 下控制台止于
        日志路径行，报错只埋在 run log，形同"无通知无报错"静默死亡）"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "测试批次")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(
                        bdir, f"NA12878_S1_L001_{r}_001.fastq.gz"), "wb") as f:
                    f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low", "--notify", "off",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            combined = r.stdout + r.stderr
            self.assertEqual(r.returncode, 1, combined[-500:])   # 批次 P0 失败
            self.assertIn("批次 测试批次 失败", combined)   # except 分支完整走完
            self.assertNotIn("UnboundLocalError", combined)   # 修复锚：不再二次崩

    def test_unrecognized_extension_skip_is_loud(self):
        """★ 跳出不静默锚（RUN-52/DEC-37）：全部输入文件未被任何布局识别（如
        .fqq.gz 这类扩展名）→ 批次跳过（退出码 0 属设计行为），但日志/通知必须
        点名未识别文件与正确扩展名——曾未匹配文件不进 invalid、跳过原因不可见，
        nohup 下终端形同静默失败"""
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "260921")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(bdir, f"SM1_{r}.fqq.gz"), "wb") as f:
                    f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low", "--notify", "off",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 0, r.stdout[-500:])   # 跳过≠失败
            self.assertIn("无任何有效样本", r.stdout)
            self.assertIn("SM1_R1.fqq.gz", r.stdout)     # 点名未识别文件
            self.assertIn("fastq.gz/.fq.gz", r.stdout)   # 指引正确扩展名


class TestExcludeMatching(unittest.TestCase):

    def test_match_excluded_token_level(self):
        """★ 对照排除 token 级匹配锚（DEC-39/RUN-55）：整串相等或 [_\-.] token
        命中（不区分大小写）——外送长前缀对照名（..._NTC_combined）曾因仅整串
        比对而漏排：排除集为空 → NTC 进联合分型且 TH-21/34 污染监控同时失效；
        部分子串不得误中（MNTCX 不中 NTC）"""
        from run_pipeline import _match_excluded
        ntc_long = ("202609182149_AE01-231101003_4P260813139US293229A2_B_"
                    "ZM20260918D_NTC_combined")
        samples = {ntc_long, "SM1", "SM2_ntc", "MNTCX", "SM3"}
        self.assertEqual(_match_excluded(samples, ["NTC"]), {ntc_long, "SM2_ntc"})
        self.assertEqual(_match_excluded(samples, ["SM1", "SM3"]), {"SM1", "SM3"})
        self.assertEqual(_match_excluded(samples, ["SM"]), set())      # 部分子串不中
        self.assertEqual(_match_excluded(samples, []), set())
        self.assertEqual(_match_excluded({"ntc"}, ["NTC"]), {"ntc"})   # 整串不区分大小写

    def test_long_prefix_ntc_excluded_e2e(self):
        """★ 长前缀 NTC 端到端锚（RUN-55）：外送平铺布局下对照名带长前缀——
        默认 --exclude-samples NTC 须 token 命中将其排除出 HaplotypeCaller
        （dry-run 命令行中不得出现该样本名）并播报对照命中行"""
        import subprocess
        ntc = ("202609182149_AE01-231101003_4P260813139US293229A2_B_"
               "ZM20260918D_NTC_combined")
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "260921")
            os.makedirs(bdir)
            for sm in ("SM1", ntc):
                for r in ("R1", "R2"):
                    with open(os.path.join(bdir, f"{sm}_{r}.fastq.gz"), "wb") as f:
                        f.write(b"@x\nACGT\n+\nIIII\n")
            env = {**os.environ, **_dep_env(td),
                   "GWAS_RESULTS": os.path.join(td, "results")}
            r = subprocess.run(
                [sys.executable,
                 os.path.join(os.path.dirname(os.path.dirname(
                     os.path.abspath(__file__))), "run_pipeline.py"),
                 "--dry-run", "--resource-profile", "low", "--notify", "off",
                 "--input", os.path.join(td, "in")],
                env=env, capture_output=True, text=True, timeout=90)
            self.assertEqual(r.returncode, 0, r.stdout[-500:])
            self.assertIn(f"对照样本（排除出联合分型，QC 照跑）: {ntc}", r.stdout)
            hc_lines = [ln for ln in r.stdout.splitlines() if "HaplotypeCaller" in ln]
            self.assertTrue(hc_lines)
            self.assertFalse(any(ntc in ln for ln in hc_lines))   # 不进 HC/联合分型


class TestInputSizeFormatting(unittest.TestCase):

    def test_fmt_gib_binary_unit(self):
        """★ 输入体量二进制口径锚（RUN-54）：_fmt_gib 按 1024³（GiB）计并带
        IEC 后缀——曾十进制 1e9（GB），与 ls -lh/du 等工具的二进制口径不一致"""
        from run_pipeline import _fmt_gib
        self.assertEqual(_fmt_gib(1 << 30), "1.00GiB")
        self.assertEqual(_fmt_gib(3 * (1 << 30) + 512 * (1 << 20)), "3.50GiB")
        self.assertEqual(_fmt_gib(5 * (1 << 20)), "0.00GiB")
        self.assertNotIn("GB", _fmt_gib(1 << 30))


class TestStartupNotifyOrder(unittest.TestCase):

    def test_startup_notify_before_md5(self):
        """★ 启动通知时序锚（DEC-38/RUN-53）：启动信息必须是第一条——先于
        md5 校验发出（md5 哈希 GB 级输入可达数分钟，期间群里应已收到启动信息）；
        曾置于 md5 之后，第一条消息被 md5 耗时阻塞且 md5 异常态提前外显冗余"""
        import argparse
        import types
        from run_pipeline import BatchCtx
        calls = []
        with tempfile.TemporaryDirectory() as td:
            bdir = os.path.join(td, "in", "B")
            os.makedirs(bdir)
            for r in ("R1", "R2"):
                with open(os.path.join(bdir, f"SM1_{r}.fastq.gz"), "wb") as f:
                    f.write(b"@x\nACGT\n+\nIIII\n")
            res = os.path.join(td, "res")
            os.makedirs(os.path.join(res, "logs"))     # Logger 落盘目标
            plan = types.SimpleNamespace(
                cpu_detected=8, mem_detected=32, profile="auto", reserve_cores=2,
                reserve_mem_gb=2, workers=1, threads=4, workers_gatk=2, workers_io=2,
                sort_mem="128M", gatk_mem="1g", cohort_mem="8g",
                peak_per_sample_gb=20)
            args = argparse.Namespace(dry_run=False, samples=None, step=6,
                                      exclude_samples="", notify="on")
            ctx = BatchCtx("B", bdir, res, args, plan, None)
            ctx.bdata = {"samples": {}}               # 平时由 process_batch 挂载
            ctx.notify_on = True

            def fake_notify(title, body, logger=None):
                calls.append(("notify", title))

            def fake_md5(*a, **k):
                calls.append(("md5", ""))
                return set(), 0, "OK"

            with mock.patch("dingtalk.notify", side_effect=fake_notify), \
                    mock.patch("scanner.verify_md5", side_effect=fake_md5), \
                    mock.patch("scanner.merge_all", return_value=({}, {})), \
                    mock.patch("run_pipeline.nonempty", return_value=True):
                try:
                    ctx.step0_scan()
                except RuntimeError:
                    pass          # merge 空结果的收尾 raise，不影响顺序断言
        self.assertTrue(calls)
        self.assertEqual(calls[0][0], "notify")        # 第一件事：启动通知
        self.assertIn("启动", calls[0][1])
        md5_i = [i for i, c in enumerate(calls) if c[0] == "md5"]
        self.assertTrue(md5_i, "md5 校验应被调用")
        self.assertLess(0, md5_i[0])                   # notify 在 md5 之前


class TestSighupHardening(unittest.TestCase):
    """★ SIGHUP 防护锚（DEC-33/RUN-48）：nohup 只让 Python 忽略 SIGHUP，
    SIG_IGN 经 fork/exec 继承本可覆盖 sh/singularity/bcftools 等全部子进程，
    但 JVM 启动时会安装自己的 SIGHUP 处理器覆盖继承位（除非 -Xrs）——
    外部服务器实跑 Step 5 四个 HC JVM 同瞬 "Hangup"（exit=129=128+SIGHUP）"""

    def test_sighup_ignored_after_guard(self):
        import signal
        from run_pipeline import ignore_sighup
        old = signal.getsignal(signal.SIGHUP)
        try:
            ignore_sighup()
            self.assertIs(signal.getsignal(signal.SIGHUP), signal.SIG_IGN)
        finally:
            signal.signal(signal.SIGHUP, old)

    def test_gatk_java_options_carry_xrs(self):
        """GATK 全部 10 类命令必须带 -Xrs（JVM 不装信号处理器，继承的
        忽略位得以保留）；旧形态 `--java-options -Xmx`（无 -Xrs）不得残留"""
        from modules import gatk as mgatk
        src = open(mgatk.__file__, encoding="utf-8").read()
        self.assertEqual(src.count('--java-options "-Xrs -Xmx'), 10)
        self.assertNotIn("--java-options -Xmx", src)

    def test_haplotypecaller_command_shape(self):
        """代表性命令形态锚：单参数引号包裹多 JVM 选项（GATK 启动器文档口径
        --java-options 'OPTION1[ OPTION2 ...]'；实测容器内跑通）"""
        from modules import gatk as mgatk
        cmds = []

        class R:
            def run(self, cmd, logger=None, outputs=(), timeout=None, capture=False):
                cmds.append(cmd)
                return 0

            def tool(self, sif_key, args, binds=None):
                return f"singularity:{sif_key} {args}"

            def cpath(self, path):
                return path

        with mock.patch.object(config, "TARGETS_INTERVAL_LIST", "/x/i.list"):
            ok = mgatk.haplotypecaller(R(), "/x/a.bam", "/x/a.g.vcf.gz",
                                       "2g", 2, _FakeLog())
        self.assertTrue(ok)
        self.assertIn('--java-options "-Xrs -Xmx2g" HaplotypeCaller', cmds[0])

    def test_fastqc_env_injects_xrs(self):
        """fastqc 同为 JVM 但启动器不收 --java-options → 命令前缀
        _JAVA_OPTIONS=-Xrs 注入（singularity 默认透传宿主环境进容器）"""
        from modules import fastqc as mfastqc
        cmds = []

        class R:
            def run(self, cmd, logger=None, outputs=(), timeout=None, capture=False):
                cmds.append(cmd)
                return 0

            def tool(self, sif_key, args, binds=None):
                return f"rt exec {args}"

            def cpath(self, path):
                return path

        self.assertTrue(mfastqc.run_fastqc(R(), ["/x/a_R1.fastq.gz"],
                                           "/x/qc", 2, _FakeLog()))
        self.assertTrue(cmds[0].startswith("_JAVA_OPTIONS=-Xrs rt exec"))
        self.assertIn("fastqc -t 2", cmds[0])

    def test_step5_split_uses_hcok(self):
        """★ 拆分名单锚（DEC-34/RUN-48）：每样本 hardfiltered/PASS 拆分必须遍历
        实际进入 cohort 的 hc_ok——曾遍历 HC 失败前的原始 calling，失败样本不在
        cohort VCF 头，bcftools view -s 连锁报"样本不存在"（exit=255）"""
        import inspect
        from run_pipeline import BatchCtx
        src = inspect.getsource(BatchCtx.step5_variant_calling)
        self.assertIn("for sm in hc_ok:", src)


class TestCohortRerunGuard(unittest.TestCase):
    """★ cohort 复跑守卫锚（DEC-34/RUN-48）：同日断点续跑补回 HC 失败样本后
    gvcf.list 样本集变化，旧 cohort 链产物非空仍被幂等 SKIP 沿用旧口径——
    补回样本 view -s 静默失败、交付残缺还报 success；守卫按现存 VCF header
    样本清单（bcftools query -l）比对，不一致即作废 cohort/matrix/
    per_sample_vcf 全部派生产物重算"""

    class _R:
        dry_run = False

        def __init__(self, samples_out):
            self._out = samples_out

        def out(self, cmd, logger=None, timeout=None):
            return self._out

        def tool(self, sif_key, args, binds=None):
            return f"rt:{sif_key} {args}"

        def cpath(self, p):
            return p

    def _mk(self, td):
        work = os.path.join(td, "results", "B_20990909")
        for d in ("cohort", "matrix", "per_sample_vcf",
                  "qc" + os.sep + "bcftools_stats", "qc" + os.sep + "multiqc"):
            os.makedirs(os.path.join(work, *d.split(os.sep)))
        paths = [os.path.join(work, "cohort", "cohort.g.vcf.gz"),
                 os.path.join(work, "cohort", "cohort.PASS.vcf.gz"),
                 os.path.join(work, "matrix", "genotype_matrix.tsv"),
                 os.path.join(work, "per_sample_vcf", "S1.PASS.adjudicated.vcf.gz"),
                 # RUN-49：cohort 级 stats 与 MultiQC 自带幂等 SKIP，须随守卫作废
                 os.path.join(work, "qc", "bcftools_stats", "cohort.raw.stats"),
                 os.path.join(work, "qc", "multiqc", "x_multiqc_report.html")]
        for p in paths:
            open(p, "wb").write(b"x")
        mdata = os.path.join(work, "qc", "multiqc", "multiqc_data")
        os.makedirs(mdata)
        open(os.path.join(mdata, "multiqc_data.json"), "wb").write(b"x")
        return work, paths

    def _guard(self, td, samples_out, hc_ok):
        from run_pipeline import cohort_rerun_guard
        work, paths = self._mk(td)
        fired = cohort_rerun_guard(
            self._R(samples_out), work, hc_ok, _FakeLog(),
            checks=(("cohort.g.vcf.gz", paths[0]),))
        return work, paths, fired

    def test_invalidate_on_sample_set_change(self):
        """旧 cohort 缺 S3（4 样本 HC 被杀后补跑场景）→ 五个派生位置清空重算
        （含 RUN-49 补漏：cohort 级 stats 与 MultiQC——二者自带幂等 SKIP，
        不作废则陈旧统计/报告被沿用进通知与交付；multiqc_data/ 子目录整树删）"""
        with tempfile.TemporaryDirectory() as td:
            work, paths, fired = self._guard(td, "S1\nS2\n", ["S1", "S2", "S3"])
            self.assertTrue(fired)
            for d in ("cohort", "matrix", "per_sample_vcf",
                      "qc" + os.sep + "bcftools_stats", "qc" + os.sep + "multiqc"):
                self.assertEqual(os.listdir(os.path.join(work, *d.split(os.sep))), [])

    def test_keep_when_same_set(self):
        """样本集一致（正常同日续跑）→ 零动作，幂等语义不变（REQ-04）；
        stats/MultiQC 同样保留（含 multiqc_data/ 子目录）"""
        with tempfile.TemporaryDirectory() as td:
            _, paths, fired = self._guard(td, "S1\nS2\n", ["S1", "S2"])
            self.assertFalse(fired)
            for p in paths:
                self.assertGreater(os.path.getsize(p), 0)

    def test_keep_when_header_unreadable(self):
        """header 读不出（dry-run / 容器异常）→ 不判不作废，零副作用"""
        with tempfile.TemporaryDirectory() as td:
            _, paths, fired = self._guard(td, "", ["S1", "S2", "S3"])
            self.assertFalse(fired)
            for p in paths:
                self.assertGreater(os.path.getsize(p), 0)


if __name__ == "__main__":
    unittest.main()
