#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fastqc：运行 + zip 内 fastqc_data.txt 解析（各模块状态汇总、Adapter 判定）。
对应笔记《1-原始数据 QC + 修剪》代码块 0/1/12。"""

import io
import zipfile
from collections import OrderedDict

import config


def run_fastqc(runner, fastqs_host, outdir_host, threads, logger):
    """fastqc -t N -o outdir <fastq...>（对原始/合并后与 clean fastq 各跑一次）。
    fastqc 同为 JVM 但其启动器不收 --java-options——SIGHUP 防护（DEC-33/RUN-48）
    经 _JAVA_OPTIONS=-Xrs 环境前缀注入（singularity 默认透传宿主环境进容器；
    代价：每条命令 stderr 多一行 "Picked up _JAVA_OPTIONS: -Xrs" 提示）"""
    outs = [f"{outdir_host}/{fq.rsplit('/', 1)[-1].replace('.fastq.gz', '')}_fastqc.zip"
            for fq in fastqs_host]
    rc = runner.run(
        "_JAVA_OPTIONS=-Xrs "
        + runner.tool("fastqc",
                      f"fastqc -t {threads} -o {runner.cpath(outdir_host)} "
                      + " ".join(runner.cpath(f) for f in fastqs_host)),
        logger=logger, outputs=outs)
    return rc == 0


def parse_fastqc_zip(zip_path):
    """解析 zip 内 fastqc_data.txt → (OrderedDict[模块→状态], basic dict)"""
    modules, basic = OrderedDict(), {}
    try:
        with zipfile.ZipFile(zip_path) as zf:
            inner = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")]
            if not inner:
                return modules, basic
            with zf.open(inner[0]) as fh:
                for raw in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
                    line = raw.rstrip("\n")
                    if line.startswith(">>"):
                        if not line.startswith(">>END_MODULE"):
                            parts = line[2:].split("\t")
                            if len(parts) >= 2:
                                modules[parts[0]] = parts[1]
                        continue
                    if line.startswith("#") or line.startswith(">"):
                        continue
                    if "\t" in line:   # Basic Statistics 键值表
                        k, _, v = line.partition("\t")
                        if k in ("Filename", "Total Sequences", "Sequence length",
                                 "%GC", "Encoding"):
                            basic[k] = v
    except (OSError, zipfile.BadZipFile):
        pass
    return modules, basic


def adapter_judgement(status):
    """Adapter Content：PASS<5%、WARN 5-20%、FAIL>20%（fastqc 已按此打标）"""
    return {"PASS": "<5%", "WARN": "5-20%", "FAIL": ">20%"}.get(status, status)


def summarize(zip_paths, logger, tag):
    """汇总多个 zip 的模块状态；列出全部非 PASS 项 + 重点模块"""
    focus = ["Basic Statistics", "Per base sequence quality", "Adapter Content",
             "Overrepresented sequences"]
    any_fail = False
    for zp in zip_paths:
        modules, basic = parse_fastqc_zip(zp)
        if not modules:
            logger.warn(f"[{tag}] 无法解析 {zp}")
            continue
        bad = {m: s for m, s in modules.items() if s != "PASS"}
        name = zp.rsplit("/", 1)[-1].replace("_fastqc.zip", "")
        logger.result(
            f"[{tag}] {name}: " + (", ".join(f"{m}={s}" for m, s in bad.items()) if bad
                                   else "全模块 PASS")
            + f" | reads={basic.get('Total Sequences', '?')}"
            + (f" | Adapter({adapter_judgement(modules.get('Adapter Content', '?'))})"
               if "Adapter Content" in modules else ""))
        any_fail = any_fail or any(s == "FAIL" for s in modules.values())
    return any_fail


def check_adapter_cleared(raw_zip, trim_zip, logger):
    """FastQC(trim) 复检：确认 Adapter 转 PASS（fastqc 状态值为小写）"""
    raw_status = parse_fastqc_zip(raw_zip)[0].get("Adapter Content") if raw_zip else None
    trim_status = parse_fastqc_zip(trim_zip)[0].get("Adapter Content") if trim_zip else None
    if raw_status and trim_status:
        if trim_status.upper() == "PASS":
            logger.result(f"Adapter 复检转 PASS（raw={raw_status} → trim={trim_status}）")
        else:
            logger.warn(f"Adapter 复检未转 PASS（raw={raw_status} → trim={trim_status}）")
