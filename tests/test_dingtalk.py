#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dingtalk.py 单测（消息结构规范化纯函数 + 企业机器人发送链路 mock）——
官方 markdown 子集不支持表格且换行需 \\n\\n；发送用例 mock _request 不发网络）"""

import json
import os
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dingtalk                                   # noqa: E402
from dingtalk import md_table_to_list, normalize_markdown, MAX_TEXT_BYTES  # noqa: E402


class TestTableDowngrade(unittest.TestCase):

    def test_two_col_table_to_list(self):
        src = "| 指标 | 值 |\n| --- | --- |\n| mapped | 99.9% |\n| dup | 2.0% |"
        out = md_table_to_list(src)
        self.assertNotIn("|", out.replace("-", ""))
        self.assertNotIn("---", out)
        self.assertIn("- **指标**: mapped｜**值**: 99.9%", out)
        self.assertIn("- **指标**: dup｜**值**: 2.0%", out)

    def test_table_inside_text(self):
        src = "#### 标题\n| a | b |\n| --- | --- |\n| 1 | 2 |\n尾巴"
        out = md_table_to_list(src)
        self.assertIn("#### 标题", out)
        self.assertTrue(any(l.startswith("- ") for l in out.splitlines()))
        self.assertIn("尾巴", out)

    def test_no_table_passthrough(self):
        src = "普通文本\n没有表格"
        self.assertEqual(md_table_to_list(src), src)


class TestNormalize(unittest.TestCase):

    def test_single_newline_becomes_double(self):
        out = normalize_markdown("第一行\n第二行")
        self.assertIn("第一行\n\n第二行", out)
        self.assertNotIn("第一行\n第二行", out)

    def test_double_newline_untouched(self):
        src = "段落一\n\n段落二"
        self.assertEqual(normalize_markdown(src), src)

    def test_table_then_newline_combined(self):
        out = normalize_markdown("| a | b |\n| --- | --- |\n| 1 | 2 |\n结语")
        self.assertNotIn("| 1 | 2 |", out)
        self.assertIn("- **a**: 1", out)
        self.assertIn("结语", out)

    def test_truncation_guard(self):
        out = normalize_markdown("x" * (MAX_TEXT_BYTES + 5000))
        self.assertLess(len(out.encode("utf-8")), MAX_TEXT_BYTES + 500)
        self.assertIn("已截断", out)

    def test_list_items_renderable(self):
        """列表项之间经 \\n\\n 分隔后仍以 '- ' 开头（钉钉可渲染为列表）"""
        out = normalize_markdown("- 甲\n- 乙\n- 丙")
        items = [l for l in out.splitlines() if l.strip()]
        self.assertTrue(all(l.strip().startswith("- ") for l in items))
        self.assertEqual(len(items), 3)


class _Log:
    def __init__(self):
        self.lines = []

    def info(self, msg):
        self.lines.append(("info", msg))

    def warn(self, msg):
        self.lines.append(("warn", msg))


def _configured_patch():
    """把 config 打成'已配置企业机器人'的样子"""
    return [mock.patch.object(dingtalk.config, "DINGTALK_CLIENT_ID", "cid_test"),
            mock.patch.object(dingtalk.config, "DINGTALK_CLIENT_SECRET", "sec_test"),
            mock.patch.object(dingtalk.config, "DINGTALK_ROBOT_CODE", ""),
            mock.patch.object(dingtalk.config, "DINGTALK_CONVERSATION_ID", "cidGROUP==")]


class TestEnterpriseSend(unittest.TestCase):
    """★ 企业机器人链路（RUN-39/DEC-24，mock _request 不发网络）：
    v1.0 groupMessages/send + msgParam JSON 字符串 + token 进程内缓存"""

    def setUp(self):
        dingtalk._reset_token_cache()
        dingtalk._no_config_warned = False

    def tearDown(self):
        # mock 的 _request 会让真 access_token() 把假 token 写进模块级缓存，
        # 不清掉会泄漏给后续用例（曾致 test_envfile 带假 token 发出真实网络请求）
        dingtalk._reset_token_cache()

    def _requests(self, responses):
        """按顺序回放 (status, body)；记录 (url, body) 供断言"""
        calls = []

        def fake(url, method="GET", data=None, headers=None, timeout=60):
            calls.append({"url": url, "method": method,
                          "data": data, "headers": headers or {}})
            return responses.pop(0)

        return calls, fake

    def test_send_markdown_group_message_and_log(self):
        calls, fake = self._requests([
            (200, json.dumps({"accessToken": "T1", "expireIn": 7200})),
            (200, json.dumps({"processQueryKey": "k"})),
        ])
        log = _Log()
        with mock.patch.object(dingtalk, "_request", side_effect=fake):
            for p in _configured_patch():
                p.start()
            try:
                ok, err = dingtalk.send_markdown("[GWAS][OK] 标题X", "正文", logger=log)
            finally:
                mock.patch.stopall()
        self.assertTrue(ok, err)
        self.assertEqual(calls[0]["url"], "https://api.dingtalk.com/v1.0/oauth2/accessToken")
        self.assertTrue(calls[1]["url"].endswith("/v1.0/robot/groupMessages/send"))
        self.assertEqual(calls[1]["headers"].get("x-acs-dingtalk-access-token"), "T1")
        sent = json.loads(calls[1]["data"].decode())
        self.assertEqual(sent["msgKey"], "sampleMarkdown")
        self.assertEqual(sent["openConversationId"], "cidGROUP==")
        # msgParam 必须是 JSON 字符串（传对象钉钉报 invalidParameter）
        mp = json.loads(sent["msgParam"])
        self.assertEqual(mp["title"], "[GWAS][OK] 标题X")
        self.assertIn("正文", mp["text"])
        self.assertTrue(any(lv == "info" and "标题X" in m for lv, m in log.lines))

    def test_notify_success_logs_info(self):
        """RUN-29 语义保持：notify 成功落 INFO（notify=on 可事后确认）"""
        calls, fake = self._requests([
            (200, json.dumps({"accessToken": "T1", "expireIn": 7200})),
            (200, json.dumps({"processQueryKey": "k"})),
        ])
        log = _Log()
        with mock.patch.object(dingtalk, "_request", side_effect=fake):
            for p in _configured_patch():
                p.start()
            try:
                dingtalk.notify("[GWAS] 标题Y", "正文", logger=log)
            finally:
                mock.patch.stopall()
        self.assertTrue(any(lv == "info" and "标题Y" in m for lv, m in log.lines))

    def test_token_cached_across_sends(self):
        """accessToken 进程内缓存：两次发送只取一次 token（频繁取会被限流）"""
        calls, fake = self._requests([
            (200, json.dumps({"accessToken": "T1", "expireIn": 7200})),
            (200, json.dumps({"processQueryKey": "k"})),
            (200, json.dumps({"processQueryKey": "k2"})),
        ])
        with mock.patch.object(dingtalk, "_request", side_effect=fake):
            for p in _configured_patch():
                p.start()
            try:
                self.assertTrue(dingtalk.send_markdown("t1", "x")[0])
                self.assertTrue(dingtalk.send_markdown("t2", "x")[0])
            finally:
                mock.patch.stopall()
        token_calls = [c for c in calls if "accessToken" in c["url"]]
        self.assertEqual(len(token_calls), 1)

    def test_notify_unconfigured_no_network(self):
        """凭证未配置 → 不发任何网络请求、不抛异常、只 WARN 一次"""
        calls, fake = self._requests([])
        log = _Log()
        with mock.patch.object(dingtalk, "_request", side_effect=fake), \
                mock.patch.object(dingtalk.config, "DINGTALK_CLIENT_ID", ""), \
                mock.patch.object(dingtalk.config, "DINGTALK_CLIENT_SECRET", ""), \
                mock.patch.object(dingtalk.config, "DINGTALK_CONVERSATION_ID", ""):
            dingtalk.notify("t", "x", logger=log)
            dingtalk.notify("t2", "x", logger=log)
        self.assertEqual(calls, [])
        warns = [m for lv, m in log.lines if lv == "warn"]
        self.assertEqual(len(warns), 1)
        self.assertIn("DINGTALK_CLIENT_ID", warns[0])

    def test_file_reject_reason_ext_and_size(self):
        """文件硬限制：后缀白名单 + 20MB 上限（交付目录打包 zip 即为此）"""
        with tempfile.TemporaryDirectory() as td:
            zp = os.path.join(td, "d.zip")
            open(zp, "wb").write(b"x")
            self.assertIsNone(dingtalk.file_reject_reason(zp))
            fq = os.path.join(td, "s.fastq.gz")
            open(fq, "wb").write(b"x")
            self.assertIn("不支持", dingtalk.file_reject_reason(fq))
            self.assertIn("文件不存在", dingtalk.file_reject_reason(td + "/nope.zip"))
            open(zp, "wb").write(b"x" * 100)          # 超过 mock 上限
            with mock.patch.object(dingtalk, "FILE_SIZE_LIMIT", 1):
                self.assertIn("20MB", dingtalk.file_reject_reason(zp))

    def test_send_file_upload_then_card(self):
        """文件卡片两步：oapi media/upload（multipart 含文件名）→ sampleFile"""
        calls, fake = self._requests([
            (200, json.dumps({"access_token": "OT", "expires_in": 7200})),
            (200, json.dumps({"errcode": 0, "media_id": "M1"})),
            (200, json.dumps({"accessToken": "T1", "expireIn": 7200})),
            (200, json.dumps({"processQueryKey": "k"})),
        ])
        with tempfile.TemporaryDirectory() as td:
            zp = os.path.join(td, "交付.zip")
            open(zp, "wb").write(b"zip-content")
            with mock.patch.object(dingtalk, "_request", side_effect=fake):
                for p in _configured_patch():
                    p.start()
                try:
                    ok, err = dingtalk.send_file(zp, logger=_Log())
                finally:
                    mock.patch.stopall()
            self.assertTrue(ok, err)
        up = calls[1]
        self.assertIn("oapi.dingtalk.com/media/upload", up["url"])
        self.assertIn("multipart/form-data", up["headers"]["Content-Type"])
        self.assertIn('name="media"; filename="交付.zip"'.encode(), up["data"])
        self.assertIn(b"zip-content", up["data"])
        sent = json.loads(calls[3]["data"].decode())
        self.assertEqual(sent["msgKey"], "sampleFile")
        self.assertEqual(json.loads(sent["msgParam"])["mediaId"], "M1")

    def test_send_zip_dir_packages_and_sends(self):
        """★ 交付推送入口（DEC-24）：目录 → 临时 zip → 说明消息 + 文件卡片 → 清理"""
        sends = []

        def fake_send_markdown(title, text, logger=None):
            sends.append(("md", title))
            return True, None

        def fake_send_file(path, logger=None):
            sends.append(("file", os.path.basename(path)))
            return True, None

        with tempfile.TemporaryDirectory() as td:
            ddir = os.path.join(td, "260422_20260914")
            os.makedirs(ddir)
            open(os.path.join(ddir, "x.PASS.adjudicated.vcf.gz"), "wb").write(b"vcf")
            made = []

            def fake_archive(base, fmt, root_dir, base_dir):
                zp = base + "." + fmt
                with zipfile.ZipFile(zp, "w") as zf:
                    zf.write(os.path.join(root_dir, base_dir, "x.PASS.adjudicated.vcf.gz"),
                             arcname=os.path.join(base_dir, "x.PASS.adjudicated.vcf.gz"))
                made.append(zp)
                return zp

            with mock.patch.object(dingtalk, "_configured", return_value=True), \
                    mock.patch.object(dingtalk, "send_markdown",
                                      side_effect=fake_send_markdown), \
                    mock.patch.object(dingtalk, "send_file", side_effect=fake_send_file), \
                    mock.patch.object(dingtalk.shutil, "make_archive",
                                      side_effect=fake_archive):
                ok, err = dingtalk.send_zip_dir(ddir, "[GWAS] 标题", "正文")
            self.assertTrue(ok, err)
            self.assertEqual([s[0] for s in sends], ["md", "file"])
            self.assertEqual(sends[1][1], "260422_20260914.zip")
            for zp in made:     # 临时包已随 tmp 目录清理
                self.assertFalse(os.path.exists(zp))

    def test_send_zip_dir_oversize_splits_volumes(self):
        """★ >20MB → 分卷压缩发送（DEC-27，RUN-41）：每卷独立合法 zip、体积 ≤ 上限、
        partNNofMM 命名；说明消息点名卷数与"解压到同一目录"合并方法；
        单卷装不下的文件点名跳过（提示到服务器取）"""
        sent_md, sent_files, vol_files = [], [], set()

        def fake_send_markdown(title, text, logger=None):
            sent_md.append(text)
            return True, None

        def fake_send_file(path, logger=None):
            self.assertLessEqual(os.path.getsize(path), 200 * 1024)  # 每卷 ≤ 上限
            with zipfile.ZipFile(path) as zf:                         # 每卷是合法 zip
                self.assertIsNone(zf.testzip())
                vol_files.update(zf.namelist())
            sent_files.append(os.path.basename(path))
            return True, None

        with tempfile.TemporaryDirectory() as td:
            ddir = os.path.join(td, "B_20260101")
            os.makedirs(ddir)
            blob = os.urandom(64 * 1024)          # 随机不可压 → zip 体积 ≈ 文件和
            for i in range(4):                    # 4×64KB，预算 190KB → 每卷 2 文件
                open(os.path.join(ddir, f"f{i}.bin"), "wb").write(blob)
            open(os.path.join(ddir, "huge.bin"), "wb").write(os.urandom(200 * 1024))
            log = _Log()
            with mock.patch.object(dingtalk, "_configured", return_value=True), \
                    mock.patch.object(dingtalk, "send_markdown",
                                      side_effect=fake_send_markdown), \
                    mock.patch.object(dingtalk, "send_file",
                                      side_effect=fake_send_file), \
                    mock.patch.object(dingtalk, "FILE_SIZE_LIMIT", 200 * 1024):
                ok, _ = dingtalk.send_zip_dir(ddir, "t", "正文", logger=log)
            self.assertTrue(ok)
            self.assertEqual(sent_files, ["B_20260101_part01of02.zip",
                                          "B_20260101_part02of02.zip"])
            self.assertEqual(vol_files,
                             {f"B_20260101/f{i}.bin" for i in range(4)})  # 无丢失无重复
            self.assertIn("分 2 卷", sent_md[0])
            self.assertIn("解压到同一目录", sent_md[0])
            self.assertIn("huge.bin", sent_md[0])     # 单卷装不下 → 消息点名
            self.assertTrue(any("分卷" in m for lv, m in log.lines if lv == "warn"))

    def test_send_zip_dir_volume_cap_falls_back(self):
        """★ 卷数超 MAX_VOLUMES → 回落纯说明消息（钉钉 20 条/分钟限流防刷屏），
        一个文件卡片都不发"""
        sent_md = []

        def fake_send_markdown(title, text, logger=None):
            sent_md.append(text)
            return True, None

        with tempfile.TemporaryDirectory() as td:
            ddir = os.path.join(td, "B_20260101")
            os.makedirs(ddir)
            blob = os.urandom(64 * 1024)
            for i in range(4):
                open(os.path.join(ddir, f"f{i}.bin"), "wb").write(blob)
            log = _Log()
            with mock.patch.object(dingtalk, "_configured", return_value=True), \
                    mock.patch.object(dingtalk, "send_markdown",
                                      side_effect=fake_send_markdown), \
                    mock.patch.object(dingtalk, "send_file") as sf, \
                    mock.patch.object(dingtalk, "FILE_SIZE_LIMIT", 200 * 1024), \
                    mock.patch.object(dingtalk, "MAX_VOLUMES", 1):
                ok, _ = dingtalk.send_zip_dir(ddir, "t", "正文", logger=log)
            self.assertTrue(ok)
            sf.assert_not_called()
            self.assertIn("超限", sent_md[0])
            self.assertTrue(any("卷数超上限" in m for lv, m in log.lines
                                if lv == "warn"))


if __name__ == "__main__":
    unittest.main()
