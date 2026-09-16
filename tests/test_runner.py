#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""runner.py 单测（对应踩坑：binds 变量遮蔽导致 tool() 崩溃、幂等 SKIP 语义、
dry-run 不执行、超时 kill、capture 通道）"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from runner import Runner, nonempty          # noqa: E402
import config                                 # noqa: E402


class TestRunner(unittest.TestCase):

    def setUp(self):
        self.r = Runner()

    def test_tool_no_shadowing_crash(self):
        """踩坑回归：tool() 内 binds 变量遮蔽曾致 ValueError: too many values to unpack"""
        cmd = self.r.tool("bcftools", "bcftools view -h /data/x.vcf.gz")
        self.assertIn("exec", cmd)
        self.assertIn(f"--bind {config.WORK_DIR}:/data", cmd)
        self.assertIn("bcftools_1.24.sif", cmd)
        self.assertIn("bcftools view -h /data/x.vcf.gz", cmd)

    def test_rt_resolved_absolute(self):
        """踩坑回归：PATH 受限环境启动曾致 `singularity: not found` exit=127
        （实测 /usr/local/bin 不在调用方 PATH）——rt 必须自解析为绝对路径"""
        self.assertTrue(os.path.isabs(self.r.rt))
        self.assertTrue(os.access(self.r.rt, os.X_OK))
        self.assertIn("exec", self.r.tool("bcftools", "bcftools view -h x"))

    def test_env_path_hardened(self):
        """PATH 兜底：标准系统目录必须在子进程 PATH 里（精简 env 不应再炸）"""
        dirs = self.r._env["PATH"].split(":")
        for d in ("/usr/local/bin", "/usr/bin", "/bin"):
            self.assertIn(d, dirs)
        self.assertEqual(self.r._env["LC_ALL"], "C")

    def test_tool_extra_binds(self):
        r = Runner(extra_binds=[("/mnt/data", "/illumina")])
        cmd = r.tool("bcftools", "bcftools view /illumina/Output/x.vcf")
        self.assertIn("--bind /mnt/data:/illumina", cmd)
        self.assertEqual(cmd.count("--bind"), 2)

    def test_cpath(self):
        self.assertEqual(self.r.cpath(config.GENOME_FA), "/data/reference/genome/genome.fa")
        r = Runner(extra_binds=[("/mnt/data", "/illumina")])
        self.assertEqual(r.cpath("/mnt/data/Output/x"), "/illumina/Output/x")
        self.assertEqual(r.cpath("/etc/hosts"), "/etc/hosts")   # 未挂载路径原样

    def test_nonempty(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "f")
            self.assertFalse(nonempty(p))
            open(p, "wb").close()
            self.assertFalse(nonempty(p))          # 0 字节 = 无效
            with open(p, "wb") as f:
                f.write(b"x")
            self.assertTrue(nonempty(p))

    def test_run_skip_when_outputs_exist(self):
        """幂等 SKIP：产物已存在时不执行命令（exit 7 也不跑，直接返回 0）"""
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "done.txt")
            with open(out, "wb") as f:
                f.write(b"data")
            rc = self.r.run("exit 7", outputs=[out])
            self.assertEqual(rc, 0)

    def test_run_executes_and_fails_without_outputs(self):
        self.assertEqual(self.r.run("exit 7"), 7)
        self.assertEqual(self.r.run("true"), 0)

    def test_dry_run_never_executes(self):
        r = Runner(dry_run=True)
        self.assertEqual(r.run("exit 7"), 0)
        self.assertEqual(r.run("rm -rf /nonexistent-xyz"), 0)

    def test_run_timeout_kills(self):
        self.assertEqual(self.r.run("sleep 5", timeout=1), -1)

    def test_capture(self):
        rc, out, err = self.r.run("echo hi", capture=True)
        self.assertEqual(rc, 0)
        self.assertEqual(out.strip(), "hi")
        self.assertEqual(self.r.out("echo yo").strip(), "yo")


if __name__ == "__main__":
    unittest.main()
