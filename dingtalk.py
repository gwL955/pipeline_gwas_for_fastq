#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""钉钉机器人 markdown 通知（urllib.request 直连，不引入 requests）。

依据官方文档《自定义机器人发送消息的消息类型》与《企业内部机器人发送 Markdown 消息》：
  - 消息类型：text / link / markdown / actionCard / feedCard，本流程用 markdown
    （title=会话列表预览，text=消息体）
  - markdown 官方支持子集：标题(#~######)、引用(>)、文字效果(**加粗** 等)、
    链接、图片、无序列表(-)、有序列表(1.)；**不支持表格**（手机端渲染为竖线串）
  - 换行需 "\\n\\n"（单个 \\n 不生效）
  - 发送频率限制 20 条/分钟；发送失败只降级写日志，绝不中断分析流程

本模块在发送前做三层结构规范化：
  1) 段落规范：单换行 → 双换行（确保客户端真正换行）
  2) 表格降级：| 表格 | 块自动转为 "- **列名**: 值" 列表（防模板回退成表格）
  3) 长度保护：超出上限截断并附提示
"""

import json
import re
import urllib.request

import config

# markdown text 有效内容上限（字节，UTF-8）；客户端对超长消息折叠，过长的截断保尾部完整
MAX_TEXT_BYTES = 18000

# webhook 未配置时只告警一次（notify 会被高频调用，避免逐条刷屏）
_no_webhook_warned = False


def _no_webhook_err():
    return (f"未配置 DINGTALK_WEBHOOK（通知已跳过）。配置方法：编辑 "
            f"{config.ENV_FILE} 写入 DINGTALK_WEBHOOK=…（模板见 .env.example），"
            "或 export DINGTALK_WEBHOOK 环境变量")


def md_table_to_list(text):
    """把 markdown 表格块降级为钉钉可渲染的列表块。
    | a | b |\\n| --- | --- |\\n| 1 | 2 |  →  - **a**: 1（b: 2）样式逐行"""
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i].strip()
        if not (line.startswith("|") and line.endswith("|") and line.count("|") >= 2):
            out.append(lines[i])
            i += 1
            continue
        block = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            block.append(lines[i].strip())
            i += 1
        rows = []
        for b in block:
            cells = [c.strip() for c in b.strip("|").split("|")]
            if all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
                continue   # 分隔行
            rows.append(cells)
        if not rows:
            continue
        header, data = rows[0], rows[1:]
        for cells in data:
            parts = []
            for h, v in zip(header, cells):
                parts.append(f"**{h}**: {v}" if len(header) > 1 else v)
            out.append("- " + "｜".join(parts))
        out.append("")   # 列表块后补空行，确保下一段落正常
    return "\n".join(out)


def normalize_markdown(text):
    """发送前三步规范化：表格降级 → 段落换行 → 长度保护（顺序不可换：表块识别依赖连续行）"""
    # ① 表格 → 列表
    text = md_table_to_list(text)
    # ② 已是 \n\n 的段落不动；单换行提升为双换行（列表/引用行间同样需要）
    text = re.sub(r"(?<!\n)\n(?!\n)", "\n\n", text)
    # ③ 长度保护（UTF-8 字节）
    raw = text.encode("utf-8")
    if len(raw) > MAX_TEXT_BYTES:
        keep = raw[:MAX_TEXT_BYTES].decode("utf-8", errors="ignore")
        text = keep.rsplit("\n", 1)[0] + \
            f"\n\n> ⚠️ 消息超长（{len(raw)}B）已截断，完整内容见日志/报告"
    return text


def send_markdown(title, text, logger=None):
    """POST markdown 消息。返回 (ok, err)。
    未配置 DINGTALK_WEBHOOK 时不发网络请求，直接返回 (False, 配置指引)。
    关键词校验（实测）：作用于正文 text 且大小写敏感——标题可保持 [GWAS][P1]
    大写模板样式，只需正文含有小写关键词，缺失时自动补一行引用兜底。"""
    if not config.DINGTALK_WEBHOOK:
        return False, _no_webhook_err()
    text = normalize_markdown(text)
    if config.DINGTALK_KEYWORD and config.DINGTALK_KEYWORD not in text:
        text += f"\n\n> {config.DINGTALK_KEYWORD}"
    body = json.dumps({"msgtype": "markdown",
                       "markdown": {"title": title, "text": text}},
                      ensure_ascii=False).encode("utf-8")
    req = urllib.request.Request(
        config.DINGTALK_WEBHOOK, data=body,
        headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            resp_data = json.loads(resp.read().decode("utf-8"))
        ok = resp_data.get("errcode") == 0
        err = None if ok else str(resp_data)
        if not ok and logger:
            logger.warn(f"钉钉通知被拒绝: {resp_data}")
        return ok, err
    except Exception as e:   # noqa: BLE001  网络/限流等，降级
        if logger:
            logger.warn(f"钉钉通知发送失败（不中断流程）: {e}")
        return False, str(e)


def notify(title, text, logger=None, enabled=True):
    """流程内通知入口：webhook 未配置时静默跳过（仅首次写一条 WARN 指引配置）。"""
    global _no_webhook_warned
    if not enabled:
        return
    if not config.DINGTALK_WEBHOOK:
        if not _no_webhook_warned:
            _no_webhook_warned = True
            if logger:
                logger.warn(f"钉钉通知未启用：{_no_webhook_err()}")
        return
    send_markdown(title, text, logger=logger)
