#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""config .env 加载器单测（环境参数三源：进程环境变量 > pipeline/.env > 内置默认；
钉钉 webhook 等环境参数曾硬编码在 config.py 默认值里，RUN-26 迁入 .env——
本用例集防止回潮并锁定解析/优先级语义）。纯标准库，无网络。"""

import importlib
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config    # noqa: E402
import dingtalk  # noqa: E402


def _write_env(content):
    fd, path = tempfile.mkstemp(suffix=".env")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(content)
    return path


class TestParseEnvFile(unittest.TestCase):

    def test_missing_file_returns_empty(self):
        self.assertEqual(config.parse_env_file("/nonexistent/.env"), {})

    def test_comments_blank_export_and_bad_lines(self):
        p = _write_env("# 整行注释\n"
                       "\n"
                       "   \n"
                       "export A=1\n"
                       "B = 2\n"
                       "没有等号的行\n"
                       "1BAD=x\n")
        try:
            parsed = config.parse_env_file(p)
        finally:
            os.unlink(p)
        self.assertEqual(parsed, {"A": "1", "B": "2"})   # 非法键/无= 行全部忽略

    def test_webhook_style_value_kept_verbatim(self):
        """含 =?& 的 URL 整串保留，不被行内注释规则误伤"""
        url = "https://oapi.dingtalk.com/robot/send?access_token=tok123&x=1"
        p = _write_env(f"DINGTALK_WEBHOOK={url}\n")
        try:
            self.assertEqual(config.parse_env_file(p)["DINGTALK_WEBHOOK"], url)
        finally:
            os.unlink(p)

    def test_inline_comment_stripped_only_when_unquoted(self):
        p = _write_env('A=1 # 标注\n'
                       'B="va # lue"\n'
                       "C='x # y'\n"
                       'D=值#紧贴\n')
        try:
            parsed = config.parse_env_file(p)
        finally:
            os.unlink(p)
        self.assertEqual(parsed["A"], "1")            # 裸值：剥 " #" 注释
        self.assertEqual(parsed["B"], "va # lue")     # 双引号：原样保留
        self.assertEqual(parsed["C"], "x # y")        # 单引号：原样保留
        self.assertEqual(parsed["D"], "值#紧贴")      # 无空格前导的 # 不算注释

    def test_empty_and_equals_in_value(self):
        p = _write_env("EMPTY=\nKEY=a=b=c\n")
        try:
            parsed = config.parse_env_file(p)
        finally:
            os.unlink(p)
        self.assertEqual(parsed["EMPTY"], "")
        self.assertEqual(parsed["KEY"], "a=b=c")      # 只按首个 = 切分


class TestEnvPrecedence(unittest.TestCase):
    """三源优先级：进程环境变量 > .env > 内置默认"""

    def test_apply_env_file_setdefault(self):
        p = _write_env("UNITTEST_FROM_FILE=yes\n"
                       "UNITTEST_EXISTING=file_value\n")
        env = dict(os.environ)
        env.pop("UNITTEST_FROM_FILE", None)
        env["UNITTEST_EXISTING"] = "process_value"
        try:
            with mock.patch.dict(os.environ, env, clear=True):
                config.parse_env_file(p)              # 纯解析不碰环境
                self.assertNotIn("UNITTEST_FROM_FILE", os.environ)
                for k, v in config.parse_env_file(p).items():
                    os.environ.setdefault(k, v)       # 与 config 导入时同一语义
                self.assertEqual(os.environ["UNITTEST_FROM_FILE"], "yes")
                self.assertEqual(os.environ["UNITTEST_EXISTING"], "process_value")
        finally:
            os.unlink(p)

    def test_config_webhook_traceable_to_source(self):
        """配置值必须可溯源：非空时等于进程环境中同名值（来自真实环境或 .env 并入）"""
        if config.DINGTALK_WEBHOOK:
            self.assertEqual(config.DINGTALK_WEBHOOK,
                             os.environ.get("DINGTALK_WEBHOOK"))

    def test_reload_webhook_default_empty_without_env(self):
        """无环境变量且 .env 不含 webhook → 默认空串（代码不再自带密钥默认值）"""
        p = _write_env("DINGTALK_KEYWORD=unittest\n")
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("DINGTALK_", "GWAS_"))}
        env["GWAS_ENV_FILE"] = p
        try:
            with mock.patch.dict(os.environ, env, clear=True):
                importlib.reload(config)
                self.assertEqual(config.DINGTALK_WEBHOOK, "")
                self.assertEqual(config.DINGTALK_KEYWORD, "unittest")  # 来自 .env
                self.assertEqual(config.ENV_FILE, p)
        finally:
            os.unlink(p)
            importlib.reload(config)   # 恢复真实配置，避免污染其余用例


class TestTargetsBedEnvOverride(unittest.TestCase):
    """★ 靶区 bed 改址锚（RUN-38）：GWAS_TARGETS_BED 只指 bed 本体，
    派生 sorted.bed/interval_list 随其同目录；未配置走内置默认 $WORK/reference/targets.bed"""

    def test_default_without_key(self):
        p = _write_env("DINGTALK_KEYWORD=unittest\n")   # 不含 GWAS_TARGETS_BED
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("GWAS_", "DINGTALK_"))}
        env["GWAS_ENV_FILE"] = p
        try:
            with mock.patch.dict(os.environ, env, clear=True):
                importlib.reload(config)
                self.assertEqual(config.TARGETS_BED,
                                 os.path.join(config.REF_DIR, "targets.bed"))
                self.assertEqual(config.TARGETS_SORTED_BED,
                                 os.path.join(config.REF_DIR, "targets.sorted.bed"))
        finally:
            os.unlink(p)
            importlib.reload(config)   # 恢复真实配置，避免污染其余用例

    def test_env_file_overrides_bed_location(self):
        p = _write_env("GWAS_TARGETS_BED=/private/panel/targets.bed\n")
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("GWAS_", "DINGTALK_"))}
        env["GWAS_ENV_FILE"] = p
        try:
            with mock.patch.dict(os.environ, env, clear=True):
                importlib.reload(config)
                self.assertEqual(config.TARGETS_BED, "/private/panel/targets.bed")
                # 派生文件与 bed 同目录（只指定 bed 本体即可）
                self.assertEqual(config.TARGETS_SORTED_BED,
                                 "/private/panel/targets.sorted.bed")
                self.assertEqual(config.TARGETS_INTERVAL_LIST,
                                 "/private/panel/targets.sorted.interval_list")
                # 参考文件目录不受影响（genome/dbsnp 等仍在 REF_DIR）
                self.assertTrue(config.GENOME_FA.startswith(config.REF_DIR))
        finally:
            os.unlink(p)
            importlib.reload(config)


class TestNoHardcodedEnvParams(unittest.TestCase):
    """防回潮锚：环境参数（密钥类）不得再硬编码进源码"""

    def test_dingtalk_webhook_not_in_source(self):
        with open(config.__file__, encoding="utf-8") as f:
            src = f.read()
        self.assertNotIn("access_token=", src)
        self.assertNotIn("oapi.dingtalk.com", src)

    def test_send_markdown_skips_network_without_credentials(self):
        """未配置企业凭证时 (False, 指引) 且不构造网络请求（v2.15.0 语义；
        防回潮：测试绝不触网——曾因 mock token 缓存泄漏真实请求钉钉 API）"""
        with mock.patch.object(config, "DINGTALK_CLIENT_ID", ""), \
                mock.patch.object(config, "DINGTALK_CLIENT_SECRET", ""), \
                mock.patch.object(config, "DINGTALK_CONVERSATION_ID", ""), \
                mock.patch.object(dingtalk, "_request") as req:
            ok, err = dingtalk.send_markdown("t", "text")
        self.assertFalse(ok)
        self.assertIn("DINGTALK_CLIENT_ID", err)
        self.assertIn(".env", err)
        req.assert_not_called()

    def test_notify_without_credentials_no_exception(self):
        with mock.patch.object(config, "DINGTALK_CLIENT_ID", ""), \
                mock.patch.object(config, "DINGTALK_CLIENT_SECRET", ""), \
                mock.patch.object(config, "DINGTALK_CONVERSATION_ID", ""):
            dingtalk._no_config_warned = False
            dingtalk.notify("t", "text")   # 不抛异常即通过（降级不中断流程）


if __name__ == "__main__":
    unittest.main()
