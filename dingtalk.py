#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""钉钉企业内部机器人通知（urllib.request 直连，纯标准库）：markdown + 文件。

v2.15.0（DEC-24/RUN-39）由自定义 webhook 机器人切换为企业内部应用机器人——
webhook 机器人不能发文件，交付产物推送需要「媒体上传 + sampleFile 文件卡片」。

配置（.env 三源，DEC-18；客户端实现参考工作区 dingtalk_test/dingtalk_bot.py）：
  DINGTALK_CLIENT_ID / DINGTALK_CLIENT_SECRET   企业应用凭证（appKey/appSecret）
  DINGTALK_ROBOT_CODE        机器人编码（通常=Client ID，缺省复用）
  DINGTALK_CONVERSATION_ID   目标群 openConversationId（形如 cidXXXX==）
  （DINGTALK_WEBHOOK / DINGTALK_KEYWORD 为 webhook 时代遗留键，已不参与发送）

官方硬限制：
  - msgParam 为 JSON 字符串且 ≤15000B → 正文截断 MAX_TEXT_BYTES 预留转义余量
  - 文件消息 ≤20MB 且后缀限 xlsx/pdf/zip/rar/doc/docx → 交付目录打包 zip 即为此；
    整包超 20MB 时分卷多发（每卷独立合法 zip，DEC-27——真分卷后缀不在白名单）
  - accessToken 有效期 7200s：进程内缓存、提前 5 分钟过期（频繁取 token 会被限流）
  - markdown 官方子集：不支持表格、单换行不生效 → normalize_markdown 三层规范化
  - 发送频率 20 条/分钟；发送失败只降级写日志，绝不中断分析流程
"""

import json
import mimetypes
import os
import re
import shutil
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
import zipfile

import config

API_BASE = "https://api.dingtalk.com"
OAPI_BASE = "https://oapi.dingtalk.com"

# markdown text 有效内容上限（字节，UTF-8）：msgParam 限 15000B，预留 JSON 转义
# （换行/引号转义会膨胀）与 title/键名开销后取 14000
MAX_TEXT_BYTES = 14000

# 钉钉 type=file 媒体上传的硬限制（清单外格式一律失败，先包成 zip）
ALLOWED_FILE_EXT = {".xlsx", ".pdf", ".zip", ".rar", ".doc", ".docx"}
FILE_SIZE_LIMIT = 20 * 1024 * 1024   # 20MB

# 交付超限分卷（DEC-27）：白名单只认 zip 等五种后缀，.zip.001/.z01 这类真分卷
# 后缀上传必被拒 → 每卷打成独立合法 zip（partNNofMM 命名，收方全下后解压到
# 同一目录即还原交付结构）。卷预算按上限 95%（压缩膨胀+zip 结构余量）；
# 卷数超上限回落纯说明消息（钉钉 20 条/分钟限流，防文件卡片刷屏）
MAX_VOLUMES = 25

# webhook 未配置时只告警一次（notify 会被高频调用，避免逐条刷屏）
_no_config_warned = False


def _no_config_err():
    return (f"未配置企业机器人凭证（通知已跳过）。配置方法：编辑 "
            f"{config.ENV_FILE} 写入 DINGTALK_CLIENT_ID / DINGTALK_CLIENT_SECRET / "
            "DINGTALK_CONVERSATION_ID（模板见 .env.example），或 export 同名环境变量")


def _configured():
    return bool(config.DINGTALK_CLIENT_ID and config.DINGTALK_CLIENT_SECRET
                and config.DINGTALK_CONVERSATION_ID)


# ── HTTP 底座 ────────────────────────────────────────────────────────────
def _request(url, method="GET", data=None, headers=None, timeout=60):
    """→ (status, body_text)；HTTPError 也归一返回（调用方按响应体判成败）"""
    req = urllib.request.Request(url, data=data, method=method)
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode("utf-8", "replace")


# ── token（进程内缓存，提前 5 分钟过期）──────────────────────────────────
_token_lock = threading.Lock()
_token, _token_exp = "", 0.0
_oapi_token, _oapi_token_exp = "", 0.0


def _reset_token_cache():
    """测试用：清空进程内 token 缓存"""
    global _token, _token_exp, _oapi_token, _oapi_token_exp
    with _token_lock:
        _token, _token_exp = "", 0.0
        _oapi_token, _oapi_token_exp = "", 0.0


def access_token():
    """v1.0 accessToken（groupMessages/send 用）"""
    global _token, _token_exp
    with _token_lock:
        if _token and time.time() < _token_exp:
            return _token
        status, body = _request(
            f"{API_BASE}/v1.0/oauth2/accessToken", method="POST",
            data=json.dumps({"appKey": config.DINGTALK_CLIENT_ID,
                             "appSecret": config.DINGTALK_CLIENT_SECRET}).encode(),
            headers={"Content-Type": "application/json"})
        payload = json.loads(body or "{}")
        token = payload.get("accessToken")
        if not token:
            raise RuntimeError(f"获取 accessToken 失败（HTTP {status}）：{body[:300]}")
        _token = token
        _token_exp = time.time() + int(payload.get("expireIn", 7200)) - 300
        return _token


def _oapi_access_token():
    """旧版 oapi access_token（官方文档明确 media/upload 用 gettoken 的 token）"""
    global _oapi_token, _oapi_token_exp
    with _token_lock:
        if _oapi_token and time.time() < _oapi_token_exp:
            return _oapi_token
        qs = urllib.parse.urlencode({"appkey": config.DINGTALK_CLIENT_ID,
                                     "appsecret": config.DINGTALK_CLIENT_SECRET})
        status, body = _request(f"{OAPI_BASE}/gettoken?{qs}")
        payload = json.loads(body or "{}")
        token = payload.get("access_token")
        if not token:
            raise RuntimeError(f"获取 oapi access_token 失败（HTTP {status}）：{body[:300]}")
        _oapi_token = token
        _oapi_token_exp = time.time() + int(payload.get("expires_in", 7200)) - 300
        return _oapi_token


# ── markdown 结构规范化（沿用：表格降级 → 段落换行 → 长度保护）─────────
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


# ── 发送 ────────────────────────────────────────────────────────────────
def _group_send(msg_key, msg_param, logger):
    """v1.0 群消息统一入口（msgParam 必须是 JSON 字符串，传对象钉钉报
    invalidParameter）；→ (ok, err)"""
    try:
        msg_param_str = json.dumps(msg_param, ensure_ascii=False)
        if len(msg_param_str.encode()) > 15000:
            return False, f"msgParam 超过 15000B 限制（{len(msg_param_str.encode())}）"
        body_obj = {
            "robotCode": config.DINGTALK_ROBOT_CODE or config.DINGTALK_CLIENT_ID,
            "openConversationId": config.DINGTALK_CONVERSATION_ID,
            "msgKey": msg_key,
            "msgParam": msg_param_str,
        }
        status, body = _request(
            f"{API_BASE}/v1.0/robot/groupMessages/send", method="POST",
            data=json.dumps(body_obj, ensure_ascii=False).encode(),
            headers={"Content-Type": "application/json",
                     "x-acs-dingtalk-access-token": access_token()})
        payload = json.loads(body or "{}")
        # 成功返回 processQueryKey；失败返回 code/message
        if status != 200 or payload.get("code"):
            err = f"发送失败（HTTP {status}）：{body[:300]}"
            if logger:
                logger.warn(f"钉钉通知被拒绝（不中断流程）: {err}")
            return False, err
        return True, None
    except Exception as e:   # noqa: BLE001  网络/限流/token 失败等，降级
        if logger:
            logger.warn(f"钉钉通知发送失败（不中断流程）: {e}")
        return False, str(e)


def send_markdown(title, text, logger=None):
    """发 markdown 群消息。返回 (ok, err)；未配置凭证时不发网络请求。
    成功落一行 INFO（RUN-29：实跑 notify=on 需可事后确认钉钉真的发出）。"""
    if not _configured():
        return False, _no_config_err()
    ok, err = _group_send("sampleMarkdown",
                          {"title": title, "text": normalize_markdown(text)},
                          logger)
    if ok and logger:
        logger.info(f"钉钉已发送: {title}")
    return ok, err


def notify(title, text, logger=None, enabled=True):
    """流程内通知入口：凭证未配置时静默跳过（仅首次写一条 WARN 指引配置）。"""
    global _no_config_warned
    if not enabled:
        return
    if not _configured():
        if not _no_config_warned:
            _no_config_warned = True
            if logger:
                logger.warn(f"钉钉通知未启用：{_no_config_err()}")
        return
    send_markdown(title, text, logger=logger)


# ── 文件（媒体上传 + sampleFile 卡片）───────────────────────────────────
def file_reject_reason(path):
    """发文件前置校验 → None 可发 / str 拒绝原因（体积、格式）"""
    if not os.path.isfile(path):
        return f"文件不存在：{path}"
    size = os.path.getsize(path)
    if size > FILE_SIZE_LIMIT:
        return (f"文件 {size / 1048576:.1f}MB 超过钉钉 20MB 上限"
                "（拆包或改外链下载方案）")
    ext = os.path.splitext(path)[1].lower()
    if ext not in ALLOWED_FILE_EXT:
        return (f"钉钉不支持 {ext or '无后缀'} 格式（支持："
                f"{', '.join(sorted(ALLOWED_FILE_EXT))}；包成 zip 最省事）")
    return None


def upload_media(path):
    """上传文件到钉钉媒体库 → mediaId（有效期约 30 天，勿长期复用）。
    官方要求 media/upload 用 oapi gettoken 的 token，v1.0 token 兜底重试。"""
    ctype = mimetypes.guess_type(os.path.basename(path))[0] or "application/octet-stream"
    boundary = "----" + uuid.uuid4().hex
    with open(path, "rb") as f:
        content = f.read()
    body = b"".join([
        f'--{boundary}\r\nContent-Disposition: form-data; name="type"\r\n\r\nfile\r\n'.encode(),
        f'--{boundary}\r\nContent-Disposition: form-data; name="media"; '
        f'filename="{os.path.basename(path)}"\r\nContent-Type: {ctype}\r\n\r\n'.encode(),
        content,
        b"\r\n--" + boundary.encode() + b"--\r\n",
    ])

    def _do(token):
        return _request(
            f"{OAPI_BASE}/media/upload?access_token={token}", method="POST",
            data=body,
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
            timeout=300)

    status, raw = _do(_oapi_access_token())
    payload = json.loads(raw or "{}")
    if payload.get("errcode") in (40014, 40001, 42001, 41001):   # token 失效兜底
        status, raw = _do(access_token())
        payload = json.loads(raw or "{}")
    media_id = payload.get("media_id")
    if not media_id:
        raise RuntimeError(f"媒体上传失败（HTTP {status}）：{raw[:300]}")
    return media_id


def send_file(path, logger=None):
    """发原生文件卡片（上传拿 mediaId → sampleFile）。返回 (ok, err)。"""
    reason = file_reject_reason(path)
    if reason:
        if logger:
            logger.warn(f"钉钉文件发送跳过：{reason}")
        return False, reason
    try:
        media_id = upload_media(path)
    except Exception as e:   # noqa: BLE001
        if logger:
            logger.warn(f"钉钉媒体上传失败（不中断流程）: {e}")
        return False, str(e)
    ok, err = _group_send(
        "sampleFile",
        {"mediaId": media_id, "fileName": os.path.basename(path),
         "fileType": os.path.splitext(path)[1].lstrip(".").lower()},
        logger)
    if ok and logger:
        logger.info(f"钉钉文件已发送: {os.path.basename(path)}")
    return ok, err


# ── 交付分卷（>20MB 时按卷拆多发，DEC-27）────────────────────────────────
def _volume_budget():
    """单卷文件体积预算（未压缩口径）：上限 95%，留压缩膨胀与 zip 结构余量"""
    return max(FILE_SIZE_LIMIT * 19 // 20, 1)


def _iter_delivery_files(directory):
    """稳定顺序（目录名/文件名排序）枚举交付文件，arcname 保留顶层目录结构"""
    root = os.path.abspath(directory)
    base = os.path.basename(os.path.normpath(directory))
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            full = os.path.join(dirpath, name)
            yield full, os.path.join(base, os.path.relpath(full, root))


def _make_volumes(directory, tmp):
    """按卷预算把交付文件分组打成多卷独立 zip → (volumes, skipped)。
    单文件超卷预算无法入卷 → skipped（说明消息点名，到服务器取）。
    DEFLATE 压缩后体积 ≤ 存储和，卷预算按未压缩口径即保证每卷 ≤ 上限。"""
    budget = _volume_budget()
    groups, cur, cur_sz, skipped = [], [], 0, []
    for full, arc in _iter_delivery_files(directory):
        sz = os.path.getsize(full)
        if sz > budget:
            skipped.append(arc)
            continue
        if cur and cur_sz + sz > budget:
            groups.append(cur)
            cur, cur_sz = [], 0
        cur.append((full, arc))
        cur_sz += sz
    if cur:
        groups.append(cur)
    stem = os.path.basename(os.path.normpath(directory))
    volumes = []
    for i, group in enumerate(groups, 1):
        path = os.path.join(tmp, f"{stem}_part{i:02d}of{len(groups):02d}.zip")
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            for full, arc in group:
                zf.write(full, arc)
        volumes.append(path)
    return volumes, skipped


def send_zip_dir(directory, title, text, logger=None):
    """目录打包 zip → 先发说明消息再发文件（两条独立消息）→ 清理临时包。
    交付推送入口（DEC-24）：zip 后缀/20MB 限制由 file_reject_reason 把关，
    整包超限 → 分卷多发（DEC-27：每卷独立合法 zip ≤ 上限，说明消息点名卷数
    与合并方法；单文件超卷预算点名跳过；卷数超 MAX_VOLUMES 回落纯说明消息）；
    任何失败不中断流程。"""
    if not _configured():
        return False, _no_config_err()
    tmp = tempfile.mkdtemp(prefix="gwas_delivery_")
    try:
        base = os.path.join(tmp, os.path.basename(directory.rstrip("/")))
        zip_path = shutil.make_archive(base, "zip", os.path.dirname(
            os.path.abspath(directory)), os.path.basename(directory.rstrip("/")))
        size_mb = os.path.getsize(zip_path) / 1048576
        reason = file_reject_reason(zip_path)
        if reason:
            volumes, skipped = _make_volumes(directory, tmp)
            if not volumes or len(volumes) > MAX_VOLUMES:
                if logger:
                    logger.warn(f"交付 zip 不发送文件卡片：{reason}（分卷不可用："
                                f"{'卷数超上限 ' + str(MAX_VOLUMES) if volumes else '无有效文件'}）")
                return send_markdown(title, text + f"\n\n> 交付 zip {size_mb:.1f}MB 超限，"
                                     "文件请到服务器 Output/ 目录获取", logger=logger)
            if logger:
                logger.warn(f"交付 zip {size_mb:.1f}MB 超限，改分卷发送（{len(volumes)} 卷）")
            note = (f"\n\n> 交付 zip {size_mb:.1f}MB 超钉钉单文件 20MB 上限，"
                    f"已分 {len(volumes)} 卷发送——**全部下载后解压到同一目录**即还原交付结构")
            if skipped:
                note += ("\n\n> ⚠️ 以下文件单卷装不下未发送，请到服务器 Output/ 获取："
                         + "、".join(skipped))
            ok1, err1 = send_markdown(title, text + note, logger=logger)
            results = [send_file(v, logger=logger) for v in volumes]
            if logger:
                logger.info(f"交付分卷推送完成: {len(volumes)} 卷 "
                            f"{sum(1 for o, _ in results if o)}/{len(volumes)} 成功"
                            f"（markdown={'OK' if ok1 else err1}）")
            return all(o for o, _ in results) and ok1, \
                next((e for e in [err1] + [err for _, err in results] if e), None)
        ok1, err1 = send_markdown(title, text, logger=logger)
        ok2, err2 = send_file(zip_path, logger=logger)
        if logger:
            logger.info(f"交付 zip 推送完成: {os.path.basename(zip_path)}"
                        f"（{size_mb:.1f}MB，markdown={'OK' if ok1 else err1}"
                        f"，file={'OK' if ok2 else err2}）")
        return ok1 and ok2, err1 or err2
    except Exception as e:   # noqa: BLE001
        if logger:
            logger.warn(f"交付 zip 推送失败（不中断流程）: {e}")
        return False, str(e)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
