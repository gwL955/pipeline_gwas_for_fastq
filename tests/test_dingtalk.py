#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dingtalk.py 单测（消息结构规范化：表格降级/换行/长度——官方 markdown 子集
不支持表格且换行需 \\n\\n；纯函数，不发网络请求）"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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


if __name__ == "__main__":
    unittest.main()
