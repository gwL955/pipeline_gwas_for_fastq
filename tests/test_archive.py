#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""archive_results.py 单测（超期批次目录筛选 / 7z 命令口径 / 幂等 / -sdel 真往返）。
命令构造用例 mock subprocess（CI 无 7z 也能跑）；真实往返 skipUnless 本机有 7z。"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from datetime import date
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import archive_results                              # noqa: E402


class TestSelectBatchDirs(unittest.TestCase):

    def test_selects_only_stale_dated_dirs(self):
        """只有 `<名>_<YYYYMMDD>` 且执行日期超期的目录入选；新目录/散文件/非法日期跳过"""
        with tempfile.TemporaryDirectory() as td:
            old = os.path.join(td, "260422_20260101")     # 远超 30 天
            fresh = os.path.join(td, "260917_" + date.today().strftime("%Y%m%d"))
            nodate = os.path.join(td, "plain_dir")
            for d in (old, fresh, nodate):
                os.makedirs(d)
            open(os.path.join(td, "260422_20991301"), "w").close()   # 非法日期+是文件
            open(os.path.join(td, "readme.md"), "w").close()
            hits = archive_results.select_batch_dirs(td, date.today(), 30)
            self.assertEqual([os.path.basename(p) for p, _ in hits], ["260422_20260101"])

    def test_days_boundary_exclusive(self):
        """恰好 31 天 >30 入选、30 天不入选（>days 判定锚）"""
        with tempfile.TemporaryDirectory() as td:
            today = date(2026, 9, 18)
            for name, d in (("A_20260819", date(2026, 8, 19)),    # 30 天 → 不选
                            ("B_20260818", date(2026, 8, 18))):   # 31 天 → 选
                os.makedirs(os.path.join(td, name))
            hits = archive_results.select_batch_dirs(td, today, 30)
            self.assertEqual([os.path.basename(p) for p, _ in hits], ["B_20260818"])


class TestArchiveOne(unittest.TestCase):

    def test_command_flags_exact(self):
        """★ 压缩口径锚（用户指定）：7z a -t7z -mx=9 -mfb=192 -ms=on -md=256m
        -snl -mmt -sdel <dest.7z> <src>（含 -sdel 压缩成功后删源）"""
        with tempfile.TemporaryDirectory() as td:
            src = os.path.join(td, "B_20260101")
            os.makedirs(src)
            out = os.path.join(td, "archive")
            argv_seen = {}

            def fake_run(argv, **kw):
                argv_seen["argv"] = argv
                return subprocess.CompletedProcess(argv, 0, "", "")

            with mock.patch.object(archive_results.subprocess, "run",
                                   side_effect=fake_run):
                status, _ = archive_results.archive_one("/usr/bin/7z", src, out)
            self.assertEqual(status, "ok")
            argv = argv_seen["argv"]
            self.assertEqual(argv[0], "/usr/bin/7z")
            self.assertEqual(argv[1], "a")
            self.assertEqual(argv[2:10],
                             ["-t7z", "-mx=9", "-mfb=192", "-ms=on",
                              "-md=256m", "-snl", "-mmt", "-sdel"])
            self.assertEqual(argv[10], os.path.join(out, "B_20260101.7z"))
            self.assertEqual(argv[11], src)

    def test_idempotent_skip_existing_archive(self):
        """同名 .7z 已存在 → 跳过（不重压不覆盖）"""
        with tempfile.TemporaryDirectory() as td:
            src = os.path.join(td, "B_20260101")
            os.makedirs(src)
            out = os.path.join(td, "archive")
            os.makedirs(out)
            open(os.path.join(out, "B_20260101.7z"), "wb").write(b"x")
            with mock.patch.object(archive_results.subprocess, "run") as run:
                status, msg = archive_results.archive_one("/usr/bin/7z", src, out)
            run.assert_not_called()
            self.assertEqual(status, "skip")

    def test_fail_keeps_source_semantics(self):
        """7z rc≠0 → fail（源目录保留由 7z -sdel 语义保证，脚本不再自删）"""
        with tempfile.TemporaryDirectory() as td:
            src = os.path.join(td, "B_20260101")
            os.makedirs(src)
            bad = subprocess.CompletedProcess([], 2, "", "disk full")
            with mock.patch.object(archive_results.subprocess, "run",
                                   return_value=bad):
                status, msg = archive_results.archive_one("/usr/bin/7z", src,
                                                          os.path.join(td, "out"))
            self.assertEqual(status, "fail")
            self.assertIn("源目录保留", msg)
            self.assertTrue(os.path.isdir(src))

    @unittest.skipUnless(shutil.which("7z"), "本机无 7z（CI 可跳过）")
    def test_real_7z_roundtrip_sdel(self):
        """真实 7z 往返：压缩成功 → 源目录被 -sdel 删除、包可列出内容（平台口径验证）"""
        with tempfile.TemporaryDirectory() as td:
            src = os.path.join(td, "B_20260101")
            os.makedirs(src)
            open(os.path.join(src, "a.txt"), "w").write("payload")
            out = os.path.join(td, "archive")
            status, msg = archive_results.archive_one(shutil.which("7z"), src, out)
            self.assertEqual(status, "ok", msg)
            self.assertFalse(os.path.exists(src), "-sdel 应删除源目录")
            zp = os.path.join(out, "B_20260101.7z")
            self.assertTrue(os.path.isfile(zp))
            lst = subprocess.run(["7z", "l", zp], capture_output=True, text=True)
            self.assertEqual(lst.returncode, 0)
            self.assertIn("a.txt", lst.stdout)


if __name__ == "__main__":
    unittest.main()
