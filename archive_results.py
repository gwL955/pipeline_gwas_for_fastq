#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""results/ 批次目录归档（DEC-25，v2.15.0）：执行日期超过 --days 天的批次目录
打包为 7z 并删除源目录（-sdel：压缩成功后才删，失败自动保留）。

压缩口径（用户指定）：
  7z a -t7z -mx=9 -mfb=192 -ms=on -md=256m -snl -mmt -sdel <dest.7z> <src_dir>

用法：
  python3 archive_results.py --dry-run     # 只列出将归档的目录（不动任何文件）
  python3 archive_results.py               # 归档 results/ 下超 30 天批次目录
  python3 archive_results.py --days 60     # 自定义保留天数
默认 results=config.RESULTS_ROOT（GWAS_RESULTS 可改）、输出=config.ARCHIVE_DIR
（默认 $WORK/archive/，GWAS_ARCHIVE_DIR 可改）。

判定与安全：
  - 只认 `<批次名>_<YYYYMMDD>` 目录名后缀的执行日期（DEC-01 口径），其余跳过
  - 幂等：输出目录已存在同名 .7z → 跳过该目录（不重压不覆盖）
  - -sdel 由 7z 保证"成功才删"；rc≠0 时源目录保留并报错，逐目录独立继续
  - 需要 7z 在 PATH（本机 /usr/bin/7z）
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config                                       # noqa: E402

# 用户指定的压缩口径（顺序即命令行顺序）；-sdel=压缩成功后删源
SEVEN_ZIP_ARGS = ["-t7z", "-mx=9", "-mfb=192", "-ms=on", "-md=256m", "-snl", "-mmt", "-sdel"]
BATCH_DIR_RE = re.compile(r"^(?P<name>.+)_(?P<date>\d{8})$")


def select_batch_dirs(results_root, today, days):
    """→ [(dir_path, 执行日期)] 目录名日期距今超过 days 天的批次目录（按名排序）"""
    out = []
    for name in sorted(os.listdir(results_root)):
        p = os.path.join(results_root, name)
        m = BATCH_DIR_RE.match(name)
        if not m or not os.path.isdir(p):
            continue
        s = m.group("date")
        try:
            d = date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        except ValueError:
            continue   # 非法日期（如 20991301）跳过
        if (today - d).days > days:
            out.append((p, d))
    return out


def archive_one(seven_zip, src_dir, archive_dir, dry_run=False):
    """打包单个批次目录 → (status, msg)；status ∈ ok/skip/dry/fail"""
    dest = os.path.join(archive_dir, os.path.basename(src_dir) + ".7z")
    if os.path.exists(dest):
        return "skip", f"已存在 {os.path.basename(dest)}（不重压不覆盖）"
    if dry_run:
        return "dry", f"将归档 {src_dir} → {dest}"
    os.makedirs(archive_dir, exist_ok=True)
    r = subprocess.run([seven_zip, "a", *SEVEN_ZIP_ARGS, dest, src_dir],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return "fail", (f"7z 失败 rc={r.returncode}（源目录保留）: "
                        f"{(r.stderr or r.stdout)[-300:]}")
    return "ok", f"{src_dir} → {dest}（-sdel 已删源目录）"


def main():
    ap = argparse.ArgumentParser(
        description="results/ 超期批次目录归档为 7z（-sdel 压缩成功后删源）")
    ap.add_argument("--days", type=int, default=30,
                    help="保留天数，执行日期距今超过该值才归档（默认 30）")
    ap.add_argument("--results", default=config.RESULTS_ROOT,
                    help=f"批次结果根目录（默认 {config.RESULTS_ROOT}）")
    ap.add_argument("--out", default=config.ARCHIVE_DIR,
                    help=f"7z 输出目录（默认 {config.ARCHIVE_DIR}）")
    ap.add_argument("--dry-run", "-n", action="store_true",
                    help="只列出将归档的目录，不执行 7z、不删除任何文件")
    args = ap.parse_args()

    if not os.path.isdir(args.results):
        print(f"[FAIL] results 目录不存在: {args.results}")
        return 2
    seven_zip = shutil.which("7z")
    if not seven_zip and not args.dry_run:
        print("[FAIL] 未找到 7z（需要 p7zip；本机预期 /usr/bin/7z）")
        return 2

    today = date.today()
    targets = select_batch_dirs(args.results, today, args.days)
    print(f"归档扫描: {args.results}｜保留 >{args.days} 天｜输出 {args.out}"
          f"｜命中 {len(targets)} 个目录")
    counts = {"ok": 0, "skip": 0, "dry": 0, "fail": 0}
    for src, d in targets:
        status, msg = archive_one(seven_zip or "7z", src, args.out, args.dry_run)
        counts[status] += 1
        print(f"[{status.upper():4s}] {os.path.basename(src)}（执行日期 {d}，"
              f"距今 {(today - d).days} 天）: {msg}")
    print(f"完成: 归档 {counts['ok']}｜跳过(已存在) {counts['skip']}｜"
          f"dry {counts['dry']}｜失败 {counts['fail']}")
    return 1 if counts["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
