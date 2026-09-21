#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWAS 靶向测序 FASTQ→VCF 通用流程主控。

用法示例：
  python3 run_pipeline.py --dry-run --batch 260422     # 打印全流程命令与资源计划
  python3 run_pipeline.py --batch 260422               # 冒烟批次（Step 0-6）
  nohup python3 run_pipeline.py --input 0_raw_data \
      > results/full_run.log 2>&1 &                    # 全量四批次（幂等断点续跑）
  python3 run_pipeline.py --resource-profile low --dry-run
  python3 run_pipeline.py --notify-test
"""

import os
import re
import sys
import glob
import json
import time
import shutil
import signal
import zipfile
import argparse
import tempfile
import threading
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import config
import resource
import scanner
import dingtalk
import alerts
import report as report_mod
from logger import (Logger, SampleLoggerFactory,
                    capture_stdio, tee_set_log_path, tee_set_echo)
from runner import Runner, nonempty
from modules import fastqc as mfastqc
from modules import fastp as mfastp
from modules import bwa_mem2 as mbwa
from modules import samtools as msam
from modules import gatk as mgatk
from modules import bcftools as mbc
from modules import mosdepth as mmos
from modules import multiqc as mmultiqc

# 名称白名单（批次名/样本名共用，DEC-19/RUN-38）：两者进入 shell 命令拼接与
# 结果路径（样本名另进 bwa @RG 头），非常规字符属注入面；中文路径在 singularity
# 容器内因 locale 报错——含中文即 P0 阻断
NAME_RE = re.compile(r"[A-Za-z0-9_.\-]+")


def parse_args():
    p = argparse.ArgumentParser(
        description="GWAS 靶向测序 FASTQ→VCF 通用流程（批次独立、幂等续跑、钉钉通知）")
    p.add_argument("--input", default=config.RAW_DATA_DIR,
                   help="多批次输入根目录（默认 0_raw_data，遍历其下全部批次子目录）")
    p.add_argument("--batch", help="指定单个批次名（与 --input 同时给出时优先）")
    p.add_argument("--step", type=int, default=6, choices=range(0, 7),
                   help="执行到第几步（0=清点合并 … 6=质量汇总，默认 6）")
    p.add_argument("--samples", help="仅处理指定样本（逗号分隔）")
    p.add_argument("--exclude-samples", default="NTC",
                   help="从联合变异检测排除的对照样本（逗号分隔，默认 NTC；空串关闭）")
    p.add_argument("--workers", type=int, help="并行样本数（覆盖资源规划，最高优先）")
    p.add_argument("--threads", type=int, help="每样本线程（覆盖资源规划）")
    p.add_argument("--max-memory", type=int, help="覆盖可用内存探测值（NG）")
    p.add_argument("--resource-profile", default="auto", choices=["auto", "low", "high"],
                   help="资源档位：auto=实测（默认）/low=16线程20G/high=100线程900G")
    p.add_argument("--serial", action="store_true", help="强制单样本串行")
    p.add_argument("--dry-run", action="store_true", help="打印命令与资源计划，不执行")
    p.add_argument("--notify", default="on", choices=["on", "off"], help="钉钉通知开关")
    p.add_argument("--notify-test", action="store_true", help="发送测试通知后退出")
    p.add_argument("--version", action="version",
                   version=f"%(prog)s {config.PIPELINE_VERSION}",
                   help="显示流程版本后退出")
    p.add_argument("--out", help="覆盖批次结果目录（默认 results/<批次>_<执行日期>）")
    return p.parse_args()


def _parallel(jobs, workers):
    """jobs: {sm: fn}；返回 {sm: fn 返回值}，异常捕获为 ('EXC', traceback)"""
    results = {}
    if not jobs:
        return results
    if workers <= 1 or len(jobs) == 1:
        for sm, fn in jobs.items():
            try:
                results[sm] = fn()
            except Exception:   # noqa: BLE001
                results[sm] = ("EXC", traceback.format_exc())
        return results
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(fn): sm for sm, fn in jobs.items()}
        for fut in as_completed(futs):
            sm = futs[fut]
            try:
                results[sm] = fut.result()
            except Exception:   # noqa: BLE001
                results[sm] = ("EXC", traceback.format_exc())
    return results


def _fmt_gb(n):
    return f"{n / 1e9:.2f}GB"


class BatchCtx:
    """单批次上下文：目录、runner、日志、计划、指标收集"""

    def __init__(self, batch, batch_dir, work, args, plan, main_logger):
        self.batch, self.batch_dir, self.work, self.args, self.plan = \
            batch, batch_dir, work, args, plan
        self.dry_run = args.dry_run
        # dry-run 完全不落盘：批次/样本日志仅控制台（由全局 tee 进 run 日志）
        self.log = Logger(None if self.dry_run else
                          os.path.join(work, "logs",
                                       f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"))
        self.factory = SampleLoggerFactory()
        self.runner = Runner(dry_run=args.dry_run)
        self.failed = {}           # sm → 失败步骤
        self.metrics = {}
        self.started = time.time()
        # 跨步骤状态（v2.7.0 拆分后由各 step 方法读写）
        self.bdata = None            # 批次结果 dict（process_batch 挂载）
        self.notify_on = False
        self.run_date = ""
        self.merged = {}             # Step0 合并成功样本 → (r1, r2)
        self.excluded = set()        # 联合检测排除的对照样本
        self.cohort_stats = {}       # Step5 cohort 级统计
        self._t0 = 0.0               # 当前步骤起始时间

    def slog(self, sm):
        if self.dry_run:
            return Logger(None, prefix=sm)
        return self.factory.get(sm, os.path.join(self.work, "logs"))

    def disk_guard(self, after_step):
        """P0 运行中磁盘复查（RUN-34）：Step 间检查剩余磁盘，低于阈值立即终止
        批次——长批次开跑前只查一次，中途磁盘满仅以"产物缺失"WARN 收场为时已晚；
        dry-run 零落盘不检查"""
        if self.runner.dry_run:
            return
        p = self.work if os.path.isdir(self.work) else config.RESULTS_ROOT
        if not os.path.isdir(p):   # 沿祖先取存在路径（与开跑前检查同策略，
            p = os.path.dirname(os.path.abspath(p))   # 新机器首轮 results/ 未建也不崩）
            while p and not os.path.isdir(p):
                p = os.path.dirname(p)
        free_gb = shutil.disk_usage(p or "/").free / 1e9
        if free_gb < config.DISK_MIN_FREE_GB:
            raise RuntimeError(
                f"运行中磁盘复查 P0（{after_step} 后）：剩余 {free_gb:.0f}GB < "
                f"{config.DISK_MIN_FREE_GB}GB，终止批次（防写入中途磁盘满）")

    def _ntc_reads_anoms(self):
        """NTC reads 占批次中位样本比例异常 → P2（污染维度之二，RUN-34）"""
        reads_m = self.metrics.get("reads") or {}
        calling_reads = sorted(v for sm, v in reads_m.items()
                               if sm not in self.excluded and isinstance(v, (int, float)))
        med = calling_reads[len(calling_reads) // 2] if calling_reads else None
        anoms = []
        for sm in sorted(self.excluded & set(reads_m)):
            anoms += alerts.check_ntc_reads(reads_m[sm], med)
        return anoms

    def step_time(self, name):
        """记录步骤耗时（self._t0 由各 step 方法起始处设置）"""
        self.bdata["steps"][name] = {"duration_s": round(time.time() - self._t0, 1)}

    def makedirs(self):
        w = self.work
        if self.runner.dry_run:
            return   # dry-run 不落任何目录/文件（日志文件由 Logger 自建）
        for d in ["fastq_merged", "fastq_clean",
                  "bam", "gvcf", "cohort", "matrix", "per_sample_vcf",
                  "logs",
                  "qc/fastqc_raw", "qc/fastqc_trim", "qc/fastp", "qc/flagstat",
                  "qc/stats", "qc/hsmetrics", "qc/mosdepth", "qc/bcftools_stats",
                  "qc/multiqc"]:
            os.makedirs(os.path.join(w, d), exist_ok=True)


    # ── step0_scan ──
    def step0_scan(self):
        """Step 0：清点/md5/磁盘依赖与样本名预检（DEC-19：P0 阻断）/Lane 合并/samples.tsv。
        启动通知先于 md5 校验发出——第一条消息必须是批次启动信息（DEC-38）。
        无有效样本时返回 bdata（批次 skipped），正常返回 None。"""
        self._t0 = time.time()
        scan = scanner.scan_batch(self.batch_dir)
        valid, invalid = scan
        # Undetermined 等非样本条目（DEC-31）：剔除不进分析，记 INFO 日志并写入
        # run_summary 的 samples.ignored 供追溯（不算 invalid、不触发告警）
        self.bdata["samples"]["ignored"] = list(scan.ignored)
        for sm in scan.ignored:
            self.log.info(f"忽略非样本条目: {sm}（Illumina 下机自带未匹配 index reads，"
                          f"不进分析，DEC-31）")
        if self.args.samples:
            keep = {s.strip() for s in self.args.samples.split(",") if s.strip()}
            for sm in list(valid):
                if sm not in keep:
                    invalid[sm] = "未在 --samples 白名单"
                    del valid[sm]
        self.bdata["samples"]["valid"] = sorted(valid)
        self.bdata["samples"]["invalid"] = invalid
        self.bdata["samples"]["unmatched"] = list(scan.unmatched)   # DEC-37 追溯
        if scan.unmatched:
            self.log.info(f"未识别为样本的文件 {len(scan.unmatched)} 个"
                          f"（不构成样本，忽略；扩展名应为 .fastq.gz/.fq.gz，DEC-37）")
        if not valid:
            um = scan.unmatched
            um_note = ""
            if um:
                shown = "、".join(um[:5]) + (f" 等 {len(um)} 个" if len(um) > 5 else "")
                um_note = (f"\n\n- **未识别文件（不构成样本，扩展名应为"
                           f" .fastq.gz/.fq.gz）**：{shown}")
                self.log.warn(f"未识别文件 {len(um)} 个: {shown}")
            self.log.warn(f"批次 {self.batch} 无任何有效样本，跳过该批次（退出码仍为 0）")
            if self.notify_on:
                reasons = "".join(f"\n\n- **{k}**：{v}" for k, v in invalid.items())
                dingtalk.notify(f"批次 {self.batch} 跳过",
                                f"#### 批次 {self.batch} 跳过（无有效样本）"
                                f"\n\n> 无效样本清单如下，退出码保持 0{reasons}{um_note}",
                                logger=self.log)
            self.step_time("step0")
            return self.bdata
        for sm, reason in invalid.items():
            self.log.warn(f"无效样本（跳过）: {sm} —— {reason}")

        input_bytes = sum(os.path.getsize(f) for si in valid.values()
                          for f in si.r1 + si.r2)
        lanes_expected = sum(len(si.r1) for si in valid.values())
        n_initial = len(valid)
        self.log.info(f"有效样本 {n_initial} 个，输入体量 {_fmt_gb(input_bytes)}")

        # ── 开跑前检查（磁盘 / 依赖文件 / 样本名规范）──
        pre_anoms = []
        disk_path = self.work if os.path.isdir(self.work) else config.RESULTS_ROOT
        if not os.path.isdir(disk_path):   # dry-run 零落盘不建目录：沿祖先取存在的路径
            disk_path = os.path.dirname(os.path.abspath(disk_path))
            while disk_path and not os.path.isdir(disk_path):
                disk_path = os.path.dirname(disk_path)
        free_gb = shutil.disk_usage(disk_path or "/").free / 1e9
        need_gb = round(n_initial * config.DISK_PER_SAMPLE_GB)
        if free_gb < config.DISK_MIN_FREE_GB:
            pre_anoms.append(("P0", f"磁盘剩余 {free_gb:.0f}GB < {config.DISK_MIN_FREE_GB}GB"
                                    f"（{n_initial} 样本估算约需 {need_gb}GB）"))
        dep_files = ([config.GENOME_FA, config.GENOME_FA + ".fai", config.GENOME_DICT,
                      config.TARGETS_BED, config.DBSNP_VCF, config.MILLS_VCF]
                     + list(config.SIF.values()))
        missing_deps = [p for p in dep_files if not nonempty(p)]
        if missing_deps:
            pre_anoms.append(("P0", "依赖文件缺失或 0 字节: "
                              + ", ".join(os.path.basename(p) for p in missing_deps)))
        # 批次名/样本名白名单（DEC-19，v2.8.0 起 P0 阻断；RUN-38 补批次名）：
        # 两者都进入 shell 命令拼接与结果路径，样本名另进 bwa @RG 头——非常规
        # 字符属注入面；中文路径在 singularity 容器内因 locale 直接报错（用户实测）
        bad_batch = None if NAME_RE.fullmatch(self.batch) else self.batch
        bad_names = [sm for sm in valid if not NAME_RE.fullmatch(sm)]
        if bad_batch:
            pre_anoms.append(("P0", f"批次名 {bad_batch} 含中文或非常规字符"
                              f"（须 A-Za-z0-9_.-，容器环境不支持中文路径）"))
        if bad_names:
            pre_anoms.append(("P0", "样本名含中文或非常规字符"
                              "（须 A-Za-z0-9_.-，容器环境不支持中文路径）: "
                              + ", ".join(sorted(bad_names))))
        for lv, msg in pre_anoms:
            (self.log.error if lv == "P0" else self.log.warn)(f"[开跑前-{lv}] {msg}")

        # 启动通知（样本数/输入体量/资源计划 + 开跑前检查结论）——**第一条消息，
        # 先于 md5 校验发出（DEC-38）**：md5 哈希 GB 级输入可达数分钟，期间群里
        # 应已收到启动信息；md5 结论改由 Step0 里程碑三态（DEC-36）与
        # P0 失败通知兜底，不阻塞启动播报
        if self.notify_on:
            p = self.plan
            lvl = alerts.worst_level(pre_anoms)
            body = (f"#### [GWAS][{lvl}] {self.batch}批次 · 启动"
                    f"\n\n样本: {len(valid)}/{n_initial} 有效｜输入 {_fmt_gb(input_bytes)}｜"
                    f"{datetime.now().strftime('%m-%d %H:%M')}"
                    f"\n\n指标: 探测 {p.cpu_detected} 线程 / {p.mem_detected} GB｜"
                    f"profile {p.profile}（预留 {p.reserve_cores} 核 + {p.reserve_mem_gb} GB）"
                    f"\n\n计划: 比对 {p.workers}×{p.threads} 线程｜GATK 类 {p.workers_gatk} 路"
                    f"｜IO 类 {p.workers_io} 路｜sort -m {p.sort_mem}｜GATK -Xmx {p.gatk_mem}"
                    f"｜cohort {p.cohort_mem}｜单样本峰值 {p.peak_per_sample_gb} GB")
            for lv, msg in pre_anoms:
                body += f"\n\n异常: [{lv}] {msg}" + {"P0": " ← 阻断级（中断分析）", "P1": " ← 需确认"}.get(lv, " ← 提示")
            if not pre_anoms:
                body += "\n\n异常: 无"
            body += f"\n\n产物: {self.work}"
            dingtalk.notify(f"[GWAS][{lvl}] {self.batch}批次 · 启动", body, logger=self.log)

        # md5 完整性校验：有清单且失败=输入数据损坏，属"严重影响分析"——P0 整批
        # 阻断（v2.9.0 DEC-21，原为剔除该样本继续）；清单识别与三态
        # （OK/FAIL/SKIPPED 无清单跳过）见 DEC-36；启动通知已先行发出（DEC-38）
        md5_failed, md5_state = set(), "SKIPPED"
        if not self.args.dry_run:
            md5_failed, _, md5_state = scanner.verify_md5(self.batch_dir, self.log,
                                                          workers=self.plan.workers_io)
            for sm in sorted(md5_failed & set(valid)):
                self.log.error(f"md5 校验失败，样本 {sm} 输入数据损坏")

        # P0=阻断级（严重影响分析→直接中断批次，DEC-21）：依赖/样本名/磁盘/md5
        p0_reasons = []
        if missing_deps:
            p0_reasons.append(f"依赖文件缺失 → {'; '.join(missing_deps)}")
        if bad_batch:
            p0_reasons.append(f"批次名含中文或非常规字符（容器不支持）→ {bad_batch}")
        if bad_names:
            p0_reasons.append(f"样本名中文或非常规字符（注入面/容器不支持）→ {'; '.join(sorted(bad_names))}")
        if free_gb < config.DISK_MIN_FREE_GB:
            p0_reasons.append(f"磁盘剩余 {free_gb:.0f}GB < {config.DISK_MIN_FREE_GB}GB"
                              f"（{n_initial} 样本估算约需 {need_gb}GB）")
        if md5_failed:
            p0_reasons.append(f"md5 校验失败（输入损坏）→ {'; '.join(sorted(md5_failed))}")
        if p0_reasons:
            raise RuntimeError("开跑前检查 P0（严重影响分析，中断批次）："
                               + "；".join(p0_reasons))

        self.merged, merge_failed = scanner.merge_all(
            valid, os.path.join(self.work, "fastq_merged"), self.runner,
            self.slog, self.plan.workers_io)   # IO 类（DEC-35）
        for sm, reason in merge_failed.items():
            invalid[sm] = reason
            valid.pop(sm, None)
            self.bdata["samples"]["invalid"][sm] = reason
        if not self.merged:
            raise RuntimeError(f"批次 {self.batch} 全部样本 Lane 合并失败")
        if not self.args.dry_run:
            scanner.write_samples_tsv(os.path.join(self.work, "samples.tsv"), self.merged, valid)
            self.bdata["artifacts"]["samples_tsv"] = os.path.join(self.work, "samples.tsv")
        self.bdata["samples"]["valid"] = sorted(self.merged)
        self.log.result(f"Step 0 完成: {len(self.merged)} 样本合并就绪，samples.tsv 已生成")
        lanes_merged = sum(len(si.r1) for sm, si in valid.items() if sm in self.merged)
        self.step_time("step0")
        md5_txt = (f"FAIL {len(md5_failed)}" if md5_state == "FAIL"
                   else {"OK": "OK", "SKIPPED": "SKIPPED（无清单）"}[md5_state])
        _step_notify(self.notify_on, self.batch, 0, "清点与Lane合并", self.log,
                     total_steps=self.args.step,
                     samples=f"{len(self.merged)}/{n_initial} 成功 | Lane 合并 {lanes_merged}/{lanes_expected}",
                     metrics=f"输入 {_fmt_gb(input_bytes)} | md5 {md5_txt}"
                             f" | 无效样本 {len(invalid)}",
                     anomalies=[("P1", f"{sm} {rs}（样本终止，其余照常）")
                                for sm, rs in merge_failed.items()],
                     artifacts=f"fastq_merged/*_R*.fastq.gz "
                               + alerts.artifact_summary(os.path.join(self.work, "fastq_merged", "*_R*.fastq.gz")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")

        self.excluded = {s.strip() for s in (self.args.exclude_samples or "").split(",") if s.strip()}

    # ── step1_qc_trim ──
    def step1_qc_trim(self):
        """Step 1：FastQC(raw) → fastp → FastQC(trim)。"""
        self._t0 = time.time()
        self.log.step("Step 1: 原始 QC + 修剪")

        def s1(sm):
            slog = self.slog(sm)
            r1, r2 = self.merged[sm]
            qr1 = os.path.join(self.work, "qc", "fastqc_raw", f"{sm}_R1_fastqc.zip")
            qr2 = os.path.join(self.work, "qc", "fastqc_raw", f"{sm}_R2_fastqc.zip")
            if not mfastqc.run_fastqc(self.runner, [r1, r2],
                                      os.path.join(self.work, "qc", "fastqc_raw"),
                                      self.plan.fastqc_threads, slog):
                return sm, False, "FastQC(raw) 失败"
            mfastqc.summarize([qr1, qr2], slog, "raw")
            c1 = os.path.join(self.work, "fastq_clean", f"{sm}_R1.fastq.gz")
            c2 = os.path.join(self.work, "fastq_clean", f"{sm}_R2.fastq.gz")
            fh = os.path.join(self.work, "qc", "fastp", f"{sm}.html")
            fj = os.path.join(self.work, "qc", "fastp", f"{sm}.json")
            if not mfastp.run_fastp(self.runner, sm, r1, r2, c1, c2, fh, fj,
                                    self.plan.fastp_threads, slog):
                return sm, False, "fastp 失败"
            met = mfastp.parse_json(fj)
            slog.result(f"fastp: 保留率 {met.get('retention_pct')}% "
                        f"reads {met.get('before_reads')}→{met.get('after_reads')}")
            mfastp.check_retention(met, slog)
            t1 = os.path.join(self.work, "qc", "fastqc_trim", f"{sm}_R1_fastqc.zip")
            t2 = os.path.join(self.work, "qc", "fastqc_trim", f"{sm}_R2_fastqc.zip")
            if not mfastqc.run_fastqc(self.runner, [c1, c2],
                                      os.path.join(self.work, "qc", "fastqc_trim"),
                                      self.plan.fastqc_threads, slog):
                slog.warn("FastQC(trim) 失败（不阻断）")
            else:
                mfastqc.summarize([t1, t2], slog, "trim")
                mfastqc.check_adapter_cleared(qr1, t1, slog)
            return sm, True, met

        res = _parallel({sm: (lambda sm=sm: s1(sm)) for sm in self.merged},
                        self.plan.workers_io)   # IO 类（DEC-35）
        for sm, r in res.items():
            if isinstance(r, tuple) and len(r) == 3:
                _, ok, note = r
                if ok:
                    met = note or {}
                    self.metrics.setdefault("fastp_retention", {})[sm] = met.get("retention_pct")
                    self.metrics.setdefault("reads", {})[sm] = met.get("after_reads")
                    # Q30 百分数口径（parse_json 已换算）；整体优先，回退 R1/R2 均值
                    q30 = met.get("q30_pct")
                    if q30 is None and met.get("q30_pct_r1") is not None:
                        q1 = met["q30_pct_r1"]
                        q2 = met.get("q30_pct_r2", q1)
                        q30 = round((q1 + q2) / 2, 2)
                    if q30 is not None:
                        self.metrics.setdefault("q30_pct", {})[sm] = q30
                else:
                    self.failed[sm] = f"step1: {note}"
                    self.log.error(f"{sm} Step 1 失败: {note}")
            else:
                self.failed[sm] = f"step1 异常: {r}"
                self.log.error(f"{sm} Step 1 异常\n{r[1] if isinstance(r, tuple) else ''}")
        self.step_time("step1")
        self.log.result(f"Step 1 完成: 成功 {len(self.merged) - len(self.failed)}/{len(self.merged)}")
        _step_notify(self.notify_on, self.batch, 1, "QC+修剪", self.log,
                     total_steps=self.args.step,
                     samples=f"{len(self.merged) - len(self.failed)}/{len(self.merged)} 成功",
                     # fastp 保留率/Q30 均值保留含对照原口径（DEC-32：修剪口径对
                     # 对照同样成立；reads 低已由 DEC-31 OK 提示行专属播报）
                     metrics=f"fastp 保留率 {_avg(self.metrics.get('fastp_retention'))}% | "
                             f"Q30 {_avg(self.metrics.get('q30_pct'))}%",
                     anomalies=_step_anoms(self, "step1")
                               + alerts.check_fastp(self.metrics.get("fastp_retention"),
                                                    self.metrics.get("q30_pct"),
                                                    excluded=self.excluded)
                               + alerts.check_reads_low(self.metrics.get("reads"),
                                                        excluded=self.excluded)
                               + self._ntc_reads_anoms(),
                     artifacts=f"fastq_clean/*_R*.fastq.gz "
                               + alerts.artifact_summary(os.path.join(self.work, "fastq_clean", "*_R*.fastq.gz")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")


    # ── step2_align ──
    def step2_align(self):
        """Step 2：bwa-mem2 比对 + sort + flagstat/stats 质检。"""
        self._t0 = time.time()
        self.log.step(f"Step 2: bwa-mem2 比对（比对类 {self.plan.workers} workers × {self.plan.threads} 线程，"
                 f"sort -m {self.plan.sort_mem}）")

        def s2(sm):
            slog = self.slog(sm)
            c1 = os.path.join(self.work, "fastq_clean", f"{sm}_R1.fastq.gz")
            c2 = os.path.join(self.work, "fastq_clean", f"{sm}_R2.fastq.gz")
            bam_dir = os.path.join(self.work, "bam", sm)
            if not self.runner.dry_run:
                os.makedirs(bam_dir, exist_ok=True)
            sort_bam = os.path.join(bam_dir, f"{sm}.sort.bam")
            cmd = mbwa.build_align_pipe(self.runner, sm, c1, c2, sort_bam,
                                        self.plan.threads, self.plan.threads, self.plan.sort_mem)
            rc = self.runner.run(cmd, logger=slog, outputs=[sort_bam])
            if rc != 0:
                return sm, False, "比对/排序失败"
            if not msam.run_index(self.runner, sort_bam, slog):
                return sm, False, "index 失败"
            fl = os.path.join(self.work, "qc", "flagstat", f"{sm}.sorted.flagstat")
            st = os.path.join(self.work, "qc", "stats", f"{sm}.sorted.samtools.stats")
            msam.run_flagstat(self.runner, sort_bam, fl, slog)
            msam.run_stats(self.runner, sort_bam, st, slog)
            fs = msam.parse_flagstat(fl)
            slog.result(f"mapped {fs.get('mapped_pct')}% properly_paired "
                        f"{fs.get('pp_pct')}% singletons {fs.get('sgl_pct')}%")
            # 质量口径仅记录（v2.9.0：<90 报 P1 通知，不再判样本失败，DEC-21）
            msam.qc_judgement(fs, slog)
            return sm, True, fs

        res = _parallel({sm: (lambda sm=sm: s2(sm)) for sm in self.merged
                         if sm not in self.failed}, self.plan.workers)
        for sm, r in res.items():
            if isinstance(r, tuple) and len(r) == 3:
                _, ok, fs = r
                if isinstance(fs, dict):
                    self.metrics.setdefault("mapped_pct", {})[sm] = fs.get("mapped_pct")
                    self.metrics.setdefault("pp_pct", {})[sm] = fs.get("pp_pct")
            else:
                self.failed[sm] = f"step2 异常: {r}"
                self.log.error(f"{sm} Step 2 异常\n{r[1] if isinstance(r, tuple) else ''}")
        self.step_time("step2")
        self.log.result(f"Step 2 完成: 成功 {len(self.merged) - len(self.failed)}/{len(self.merged)}")
        _step_notify(self.notify_on, self.batch, 2, "比对", self.log,
                     total_steps=self.args.step,
                     samples=f"{len(self.merged) - len(self.failed)}/{len(self.merged)} 成功",
                     metrics=f"mapped {_avg(self.metrics.get('mapped_pct'), self.excluded)}% | "
                             f"proper pair {_avg(self.metrics.get('pp_pct'), self.excluded)}%",
                     anomalies=_step_anoms(self, "step2")
                               + alerts.check_flagstat(self.metrics.get("mapped_pct"),
                                                       self.metrics.get("pp_pct"),
                                                       excluded=self.excluded)
                               + [("P1", f"{sm} mapped {v}% < {config.MAPPED_MIN_PCT}%"
                                         f"（QC 口径，报错不中断）")
                                  for sm, v in sorted((self.metrics.get("mapped_pct") or {}).items())
                                  if sm not in self.excluded
                                  and v is not None and v < config.MAPPED_MIN_PCT],
                     artifacts=f"bam/*/*.sort.bam "
                               + alerts.artifact_summary(os.path.join(self.work, "bam", "*", "*.sort.bam")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")


    # ── step3_markdup ──
    def step3_markdup(self):
        """Step 3：MarkDuplicates 去重 + metrics/flagstat/stats。"""
        self._t0 = time.time()
        self.log.step("Step 3: MarkDuplicates 去重")

        def s3(sm):
            slog = self.slog(sm)
            bam_dir = os.path.join(self.work, "bam", sm)
            sort_bam = os.path.join(bam_dir, f"{sm}.sort.bam")
            md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
            metrics_f = os.path.join(bam_dir, f"{sm}.markdup.metrics")
            if not mgatk.markdup(self.runner, sort_bam, md_bam, metrics_f,
                                 self.plan.gatk_mem, slog):
                return sm, False, None
            met = mgatk.parse_markdup_metrics(metrics_f)
            slog.result(f"READ_PAIRS={met.get('READ_PAIRS_EXAMINED')} "
                        f"DUP={_pct(met.get('PERCENT_DUPLICATION'))} "
                        f"ELS={met.get('ESTIMATED_LIBRARY_SIZE')}")
            mgatk.qc_duplication(met, slog)
            fl = os.path.join(self.work, "qc", "flagstat", f"{sm}.markdup.flagstat")
            st = os.path.join(self.work, "qc", "stats", f"{sm}.markdup.samtools.stats")
            msam.run_flagstat(self.runner, md_bam, fl, slog)
            msam.run_stats(self.runner, md_bam, st, slog)
            # 去重前 flagstat 的 duplicates 行无意义，不采集（硬性要求）
            return sm, True, met

        res = _parallel({sm: (lambda sm=sm: s3(sm)) for sm in self.merged
                         if sm not in self.failed},
                        self.plan.workers_gatk)   # GATK 单线程类（DEC-35）
        for sm, r in res.items():
            if isinstance(r, tuple) and len(r) == 3:
                _, ok, met = r
                self.metrics.setdefault("dup_pct", {})[sm] = \
                    round(met["PERCENT_DUPLICATION"] * 100, 2) \
                    if met and met.get("PERCENT_DUPLICATION") is not None else None
                if met:
                    self.metrics.setdefault("els", {})[sm] = \
                        met.get("ESTIMATED_LIBRARY_SIZE")
                if not ok:
                    self.failed[sm] = "step3: MarkDuplicates 失败"
            else:
                self.failed[sm] = f"step3 异常: {r}"
        self.step_time("step3")
        self.log.result(f"Step 3 完成: 成功 {len(self.merged) - len(self.failed)}/{len(self.merged)}")
        # ELS 播报（DEC-32）：排除对照后的 均值/方差/最低值+样本——NTC reads 近 0，
        # 其 ELS 无统计意义且必然占据最低值（误读为文库复杂度不足）
        els_str = alerts.els_summary(self.metrics.get("els"), excluded=self.excluded)
        _step_notify(self.notify_on, self.batch, 3, "去重", self.log,
                     total_steps=self.args.step,
                     samples=f"{len(self.merged) - len(self.failed)}/{len(self.merged)} 成功",
                     metrics=f"重复率 {_avg(self.metrics.get('dup_pct'), self.excluded)}% | "
                             f"ELS {els_str if els_str else 'n/a'}",
                     anomalies=_step_anoms(self, "step3")
                               + alerts.check_dup(self.metrics.get("dup_pct"),
                                                  excluded=self.excluded),
                     artifacts=f"bam/*/*.markdup.bam "
                               + alerts.artifact_summary(os.path.join(self.work, "bam", "*", "*.markdup.bam")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")


    # ── step4_bqsr ──
    def step4_bqsr(self):
        """Step 4：BQSR 校准 + 前后 flagstat 逐行一致断言。"""
        self._t0 = time.time()
        self.log.step("Step 4: BQSR 校准")

        def s4(sm):
            slog = self.slog(sm)
            bam_dir = os.path.join(self.work, "bam", sm)
            md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
            table = os.path.join(bam_dir, f"{sm}.recal.table")
            bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
            if not mgatk.base_recalibrator(self.runner, md_bam, table,
                                           self.plan.gatk_mem, slog):
                return sm, False, "BaseRecalibrator 失败"
            if not mgatk.apply_bqsr(self.runner, md_bam, table, bqsr_bam,
                                    self.plan.gatk_mem, slog):
                return sm, False, "ApplyBQSR 失败"
            fl_md = os.path.join(self.work, "qc", "flagstat", f"{sm}.markdup.flagstat")
            fl_rc = os.path.join(self.work, "qc", "flagstat", f"{sm}.recal.flagstat")
            st_rc = os.path.join(self.work, "qc", "stats", f"{sm}.recal.samtools.stats")
            msam.run_flagstat(self.runner, bqsr_bam, fl_rc, slog)
            msam.run_stats(self.runner, bqsr_bam, st_rc, slog)
            if not self.runner.dry_run and not msam.flagstat_identical(fl_md, fl_rc):
                return sm, False, "markdup 与 BQSR flagstat 不一致（判 FAIL）"
            slog.result("BQSR 前后 flagstat 逐行一致 ✓" if not self.runner.dry_run
                        else "BQSR 前后 flagstat 断言（dry-run 跳过实际比对）")
            return sm, True, None

        res = _parallel({sm: (lambda sm=sm: s4(sm)) for sm in self.merged
                         if sm not in self.failed},
                        self.plan.workers_gatk)   # GATK 单线程类（DEC-35）
        for sm, r in res.items():
            if isinstance(r, tuple) and len(r) == 3:
                _, ok, note = r
                if not ok:
                    self.failed[sm] = f"step4: {note}"
                    self.log.error(f"{sm}: {note}")
            else:
                self.failed[sm] = f"step4 异常: {r}"
        # BQSR 校准可信度（RUN-34）：RecalTable1 M 事件观测数——known-sites
        # 覆盖崩坏时骤降 → P2（报错不中断）
        if not self.runner.dry_run:
            for sm in self.merged:
                if sm in self.failed:
                    continue
                obs = mgatk.parse_recal_observations(
                    os.path.join(self.work, "bam", sm, f"{sm}.recal.table"))
                if obs is not None:
                    self.metrics.setdefault("recal_obs", {})[sm] = int(obs)
        self.step_time("step4")
        n_bqsr_ok = len([sm for sm in self.merged if sm not in self.failed])
        self.log.result(f"Step 4 完成: 成功 {n_bqsr_ok}/{len(self.merged)}")
        _step_notify(self.notify_on, self.batch, 4, "BQSR校准", self.log,
                     total_steps=self.args.step,
                     samples=f"{n_bqsr_ok}/{len(self.merged)} 成功",
                     metrics="markdup↔BQSR flagstat 逐行一致断言 "
                             f"{n_bqsr_ok}/{n_bqsr_ok} 通过",
                     anomalies=_step_anoms(self, "step4")
                               + alerts.check_recal_low(self.metrics.get("recal_obs"),
                                                        excluded=self.excluded),
                     artifacts=f"bam/*/*.markdup.BQSR.bam "
                               + alerts.artifact_summary(os.path.join(self.work, "bam", "*", "*.markdup.BQSR.bam")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")


    # ── step5_variant_calling ──
    def step5_variant_calling(self):
        """Step 5：HaplotypeCaller gVCF → 批次内联合分型 → norm → 硬过滤 → PASS
        → 双口径矩阵导出 → 对账（数量 -R 同口径/新鲜度）。"""
        self.cohort_stats = {}
        self._t0 = time.time()
        self.log.step("Step 5: 变异检测")
        # 排序保证与 VCF 样本列序（gvcf.list=sorted）一致——矩阵列映射/裁决/重建都依赖此顺序
        calling = sorted(sm for sm in self.merged
                         if sm not in self.failed and sm not in self.excluded)
        for sm in sorted(self.excluded & set(self.merged)):
            self.log.info(f"对照样本 {sm} 按配置排除出联合变异检测（纳入 QC）")
        self.bdata["samples"]["calling"] = calling
        if not calling:
            raise RuntimeError("无可用于联合变异检测的样本（全部失败或被排除）")
        if not mgatk.prep_interval_list(self.runner, self.log):
            raise RuntimeError("interval_list 准备失败")

        def s5(sm):
            slog = self.slog(sm)
            bam_dir = os.path.join(self.work, "bam", sm)
            bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
            gvcf = os.path.join(self.work, "gvcf", f"{sm}.g.vcf.gz")
            if not mgatk.haplotypecaller(self.runner, bqsr_bam, gvcf,
                                         self.plan.gatk_mem, self.plan.hc_hmm_threads, slog):
                return sm, False
            return sm, True

        res = _parallel({sm: (lambda sm=sm: s5(sm)) for sm in calling},
                        self.plan.workers_gatk)   # GATK 单线程类（DEC-35）
        hc_ok = []
        for sm, r in res.items():
            if r == (sm, True):
                hc_ok.append(sm)
            else:
                self.failed[sm] = "step5: HaplotypeCaller 失败"
                self.log.error(f"{sm} HC 失败\n{r[1] if isinstance(r, tuple) else ''}")
        if not hc_ok:
            raise RuntimeError("全部 HaplotypeCaller 失败")
        hc_ok = sorted(hc_ok)
        self.bdata["samples"]["calling"] = hc_ok   # 实际进入联合分型的样本（排序）
        self.log.result(f"HaplotypeCaller 完成 {len(hc_ok)}/{len(calling)}，进入联合分型")

        # cohort 级（串行，gvcf.list 每次重新生成；dry-run 零落盘只打印命令）
        coh = os.path.join(self.work, "cohort")
        gvcf_list = os.path.join(coh, "gvcf.list")
        if not self.runner.dry_run:
            with open(gvcf_list, "w", encoding="utf-8") as f:
                for sm in sorted(hc_ok):
                    f.write(self.runner.cpath(os.path.join(self.work, "gvcf", f"{sm}.g.vcf.gz")) + "\n")
            self.log.info(f"gvcf.list 重新生成: {len(hc_ok)} 个 gVCF（批次内联合，严禁跨批次）")
        combined = os.path.join(coh, "cohort.g.vcf.gz")
        raw_vcf = os.path.join(coh, "cohort.raw.vcf.gz")
        pass_v = os.path.join(coh, "cohort.PASS.vcf.gz")   # 前置定义：复跑守卫要查
        # cohort 复跑守卫（DEC-34/RUN-48）：同日断点续跑补回 HC 失败样本后
        # gvcf.list 样本集已变，但旧 cohort 链产物非空仍会被幂等 SKIP——旧口径
        # （缺样本）的矩阵/拆分/裁决一路沿用，补回样本静默丢失还报 success
        cohort_rerun_guard(self.runner, self.work, hc_ok, self.log,
                           checks=(("cohort.g.vcf.gz", combined),
                                   ("cohort.PASS.vcf.gz", pass_v)))
        if not mgatk.combine_gvcfs(self.runner, gvcf_list, combined,
                                   self.plan.cohort_mem, self.log):
            raise RuntimeError("CombineGVCFs 失败")
        if not mgatk.genotype_gvcfs(self.runner, combined, raw_vcf,
                                    self.plan.cohort_mem, self.log):
            raise RuntimeError("GenotypeGVCFs 失败")
        mbc.run_stats(self.runner, raw_vcf,
                      os.path.join(self.work, "qc", "bcftools_stats", "cohort.raw.stats"), self.log)
        self.cohort_stats["raw"] = mbc.parse_stats(
            os.path.join(self.work, "qc", "bcftools_stats", "cohort.raw.stats"))

        split_vcf = os.path.join(coh, "cohort.raw.split.vcf.gz")
        norm_stats_json = os.path.join(coh, "norm_stats.json")
        ok, norm_stats = mbc.norm_split(self.runner, raw_vcf, split_vcf, self.log)
        if not ok:
            raise RuntimeError("bcftools norm 失败")
        if norm_stats:   # 统计落盘，SKIP 续跑时仍可报告
            with open(norm_stats_json, "w", encoding="utf-8") as f:
                json.dump(norm_stats, f)
        else:
            try:
                with open(norm_stats_json, encoding="utf-8") as f:
                    norm_stats = json.load(f)
            except (OSError, ValueError):
                norm_stats = {}
        self.cohort_stats["norm"] = norm_stats
        self.log.result(f"norm 摊平: {norm_stats}（后续一切靶区提取加 ±100bp padding）")

        snp_v = os.path.join(coh, "cohort.snp.vcf.gz")
        indel_v = os.path.join(coh, "cohort.indel.vcf.gz")
        snp_f = os.path.join(coh, "cohort.snp.hardfiltered.vcf.gz")
        indel_f = os.path.join(coh, "cohort.indel.hardfiltered.vcf.gz")
        hard_v = os.path.join(coh, "cohort.hardfiltered.vcf.gz")
        if not mgatk.select_variants(self.runner, split_vcf, "SNP", snp_v,
                                     self.plan.gatk_mem, self.log) \
                or not mgatk.select_variants(self.runner, split_vcf, "INDEL", indel_v,
                                             self.plan.gatk_mem, self.log):
            raise RuntimeError("SelectVariants 失败")
        if not mgatk.variant_filtration(self.runner, snp_v, snp_f,
                                        config.SNP_HARD_FILTERS, self.plan.gatk_mem, self.log) \
                or not mgatk.variant_filtration(self.runner, indel_v, indel_f,
                                                config.INDEL_HARD_FILTERS,
                                                self.plan.gatk_mem, self.log):
            raise RuntimeError("VariantFiltration 失败")
        if not mbc.concat(self.runner, [snp_f, indel_f], hard_v, self.log):
            raise RuntimeError("bcftools concat 失败")
        mbc.index_tbi(self.runner, hard_v, self.log)
        filt_ok, filt_dist = mbc.filter_column_check(self.runner, hard_v, self.log)
        if not filt_ok:
            raise RuntimeError("FILTER 列出现 '.'——过滤漏跑，判 FAIL")
        self.log.result(f"FILTER 标签分布: {filt_dist}")
        if not mbc.view_pass(self.runner, hard_v, pass_v, self.log):
            raise RuntimeError("PASS 提取失败")
        mbc.index_tbi(self.runner, pass_v, self.log)

        # 双口径导出 + stats
        mbc.export_genotype_matrix(
            self.runner, hard_v, os.path.join(self.work, "matrix", "genotype_matrix.tsv"), self.log)
        mbc.export_detail_pass(
            self.runner, pass_v, os.path.join(self.work, "matrix",
                                             "genotype_detail_PASS.tsv"), self.log)
        mbc.run_stats(self.runner, hard_v, os.path.join(
            self.work, "qc", "bcftools_stats", "cohort.hardfiltered.stats"), self.log)
        mbc.run_stats(self.runner, pass_v, os.path.join(
            self.work, "qc", "bcftools_stats", "cohort.PASS.stats"), self.log)
        self.cohort_stats["hardfiltered"] = mbc.parse_stats(os.path.join(
            self.work, "qc", "bcftools_stats", "cohort.hardfiltered.stats"))
        self.cohort_stats["PASS"] = mbc.parse_stats(os.path.join(
            self.work, "qc", "bcftools_stats", "cohort.PASS.stats"))
        # 每样本 hardfiltered / PASS（未裁决版，供追溯）——名单必须是实际进入
        # cohort 的 hc_ok（DEC-34/RUN-48）：曾误用 HC 失败前的原始 calling，
        # 失败样本不在 cohort VCF 头里，bcftools view -s 连锁报"样本不存在"
        # （exit=255，无产物白跑 2 条命令）
        for sm in hc_ok:
            hf_sm = os.path.join(self.work, "per_sample_vcf", f"{sm}.hardfiltered.vcf.gz")
            ps_sm = os.path.join(self.work, "per_sample_vcf", f"{sm}.PASS.vcf.gz")
            mbc.split_sample(self.runner, hard_v, sm, hf_sm, self.slog(sm))
            mbc.index_tbi(self.runner, hf_sm, self.slog(sm))
            mbc.split_sample(self.runner, pass_v, sm, ps_sm, self.slog(sm))
            mbc.index_tbi(self.runner, ps_sm, self.slog(sm))

        tt_raw = self.cohort_stats["raw"].get("titv")
        tt_pass = self.cohort_stats["PASS"].get("titv")
        self.log.result(f"Ti/Tv: raw {tt_raw} → PASS {tt_pass}（应上升；panel 参考区间 2.5-3.5，"
                   f"看趋势不看绝对值）")
        self.step_time("step5")
        self.bdata["artifacts"].update({
            "cohort_PASS_vcf": pass_v,
            "genotype_matrix": os.path.join(self.work, "matrix", "genotype_matrix.tsv"),
            "genotype_detail_PASS": os.path.join(self.work, "matrix",
                                                 "genotype_detail_PASS.tsv")})

        # 对账·数量：矩阵行数 vs hardfiltered 限定 targets 的记录数（应一致）
        recon_anoms = _step_anoms(self, "step5")
        if not self.runner.dry_run:
            matrix_n = sum(1 for _ in open(
                os.path.join(self.work, "matrix", "genotype_matrix.tsv"), encoding="utf-8"))
            # 同口径核对（矩阵由 query -R 生成，核对也用 -R：
            # -R/-T 在区间边界对跨界 indel 的取舍不同，混用会有固有差额）
            view_n_out = self.runner.out(
                f"{self.runner.tool('bcftools', 'bcftools view -R ' + self.runner.cpath(config.TARGETS_SORTED_BED) + ' -H ' + self.runner.cpath(hard_v))} | wc -l",
                timeout=config.STATS_TIMEOUT_S)
            try:
                view_n = int(view_n_out.strip().split()[-1])
            except (ValueError, IndexError):
                view_n = None
            self.log.result(f"对账·数量: 矩阵 {matrix_n} 行 vs view -R 计数 {view_n}")
            if view_n is not None and matrix_n != view_n:
                recon_anoms.append(("P1", f"对账·数量：矩阵行数 {matrix_n} ≠ 靶区记录数 "
                                          f"{view_n}（导出可能不完整）"))
            # 对账·新鲜度：关键 VCF mtime < 本次启动（断点续跑复用旧产物）
            stale = []
            for key_vcf in [pass_v, hard_v]:
                if nonempty(key_vcf) and os.path.getmtime(key_vcf) < self.started:
                    stale.append(os.path.basename(key_vcf))
            adj_all = sorted(glob.glob(os.path.join(
                self.work, "per_sample_vcf", "*.PASS.adjudicated.vcf.gz")))
            if adj_all and all(os.path.getmtime(p) < self.started for p in adj_all):
                stale.append("每样本裁决VCF×" + str(len(adj_all)))
            if stale:
                recon_anoms.append(("P1", "对账·新鲜度：" + ", ".join(stale)
                                    + " 为历史运行产物（断点续跑复用）"))

        _step_notify(self.notify_on, self.batch, 5, "变异检测", self.log,
                     total_steps=self.args.step,
                     samples=f"gVCF {len(hc_ok)}/{len(calling)} | 联合分型样本 {len(hc_ok)}",
                     metrics=f"raw {self.cohort_stats['raw'].get('records')} → PASS "
                             f"{self.cohort_stats['PASS'].get('records')} | "
                             f"SNP/INDEL PASS {self.cohort_stats['PASS'].get('snps')}/"
                             f"{self.cohort_stats['PASS'].get('indels')} | "
                             f"Ti/Tv {tt_raw}→{tt_pass} | "
                             f"norm split/realigned {norm_stats.get('split')}/"
                             f"{norm_stats.get('realigned')}",
                     anomalies=recon_anoms,
                     artifacts=f"cohort/cohort.PASS.vcf.gz "
                               + alerts.artifact_summary(os.path.join(self.work, "cohort", "*.vcf.gz")),
                     log_hint=f"tail -f {self.work}/logs/pipeline_*.self.log")


    # ── step6_summary_delivery ──
    def step6_summary_delivery(self):
        """Step 6：mosdepth×2 + HsMetrics → 矩阵 ./. 裁决 + 每样本 PASS 重建
        → MultiQC 最终报告（全部流程结束后）→ 交付导出 Output/（DEC-02）。"""
        self._t0 = time.time()
        self.log.step("Step 6: mosdepth ×2 + HsMetrics + 矩阵裁决 + MultiQC")

        def s6(sm):
            slog = self.slog(sm)
            bam_dir = os.path.join(self.work, "bam", sm)
            md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
            bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
            ms = {}
            for tag, bam in (("md", md_bam), ("bqsr", bqsr_bam)):
                prefix = os.path.join(self.work, "qc", "mosdepth", f"{sm}.{tag}")
                if not mmos.run_mosdepth(self.runner, bam, prefix,
                                         self.plan.mosdepth_threads, slog):
                    slog.warn(f"mosdepth({tag}) 失败（不阻断）")
                    continue
                summ = mmos.parse_summary(prefix)
                ms[tag] = summ["mean"]
                slog.result(f"mosdepth[{tag}] 靶区均值={ms[tag]}×")
            hs_txt = os.path.join(self.work, "qc", "hsmetrics", f"{sm}.hs_metrics.txt")
            if not mgatk.collect_hsmetrics(self.runner, bqsr_bam, hs_txt,
                                           self.plan.gatk_mem, slog):
                return sm, False, (ms, None)
            hs = mgatk.parse_hsmetrics(hs_txt)
            mgatk.qc_hsmetrics(hs, slog, control=sm in self.excluded)
            return sm, True, (ms, hs)

        res = _parallel({sm: (lambda sm=sm: s6(sm)) for sm in self.merged
                         if sm not in self.failed},
                        self.plan.workers_gatk)   # GATK 单线程类（DEC-35）
        for sm, r in res.items():
            if isinstance(r, tuple) and len(r) == 3:
                _, ok, (ms, hs) = r
                self.metrics.setdefault("mosdepth_mean", {})[sm] = ms.get("bqsr")
                if hs:
                    self.metrics.setdefault("mean_target_coverage", {})[sm] = \
                        hs.get("MEAN_TARGET_COVERAGE")
                    self.metrics.setdefault("pct_20x", {})[sm] = \
                        round((hs.get("PCT_TARGET_BASES_20X") or 0) * 100, 2)
                    # 捕获效率告警口径 = PCT_SELECTED_BASES（DEC-29）；
                    # on_target_pct 保留为信息指标（1bp SNP panel 下 ≈0.6%，不再告警）
                    self.metrics.setdefault("pct_selected", {})[sm] = \
                        round((hs.get("PCT_SELECTED_BASES") or 0) * 100, 2)
                    self.metrics.setdefault("on_target_pct", {})[sm] = \
                        hs.get("ON_TARGET_PCT")
                if not ok:
                    self.failed[sm] = "step6: HsMetrics 失败"
            else:
                self.failed[sm] = f"step6 异常: {r}"

        # 矩阵 ./. 裁决（mosdepth bqsr regions，DP≥20 改判 0/0）+ 每样本 PASS 重建
        calling = sorted(self.bdata["samples"].get("calling") or
                         [sm for sm in self.merged if sm not in self.failed])
        adj_tsv = os.path.join(self.work, "matrix", "genotype_matrix.adjudicated.tsv")
        adj_stats = None
        if calling and self.cohort_stats and not self.runner.dry_run:
            regions_cache = {}
            for sm in calling:
                prefix = os.path.join(self.work, "qc", "mosdepth", f"{sm}.bqsr")
                regions_cache[sm] = mmos.load_regions(prefix)

            def region_of(sm, chrom, pos):
                return mmos.region_depth(regions_cache.get(sm, []), chrom, pos)

            matrix = os.path.join(self.work, "matrix", "genotype_matrix.tsv")
            out_lines, adj_stats = mbc.adjudicate_matrix(matrix, region_of, calling,
                                                     config.DP_MIN)
            with open(adj_tsv, "w", encoding="utf-8") as f:
                f.write("\n".join(out_lines) + "\n")
            self.log.result(f"矩阵 ./. 裁决: 总 ./. {adj_stats['dotdot_total']} → "
                       f"改判 0/0 {adj_stats['filled_00']} · 保留 ./. "
                       f"{adj_stats['kept_dotdot']}（DP_MIN={config.DP_MIN}）")
            self.bdata["artifacts"]["genotype_matrix_adjudicated"] = adj_tsv

            # 每样本 PASS VCF 重建（bcftools view -s 拆分 + GT 替换，其余字段原样保留）
            pass_v = os.path.join(self.work, "cohort", "cohort.PASS.vcf.gz")
            adj_vcfs = {}
            import tempfile
            import shlex as _shlex
            for sm in calling:
                slog = self.slog(sm)
                out_vcf = os.path.join(self.work, "per_sample_vcf",
                                       f"{sm}.PASS.adjudicated.vcf.gz")
                if nonempty(out_vcf) and nonempty(out_vcf + ".tbi"):
                    adj_vcfs[sm] = out_vcf
                    continue
                with tempfile.TemporaryDirectory() as td:
                    raw_plain = os.path.join(td, f"{sm}.raw.vcf")
                    gt_tsv = os.path.join(td, f"{sm}.gt.tsv")
                    adj_plain = os.path.join(td, f"{sm}.adj.vcf")
                    rc = self.runner.run(
                        f"{self.runner.tool('bcftools', 'bcftools view -s ' + _shlex.quote(sm) + ' --min-ac 0 ' + self.runner.cpath(pass_v))} > {raw_plain}",
                        logger=slog)
                    if rc != 0:
                        continue
                    with open(adj_tsv, encoding="utf-8") as fa, \
                            open(gt_tsv, "w", encoding="utf-8") as fo:
                        idx = calling.index(sm) + 4
                        for line in fa:
                            p = line.rstrip("\n").split("\t")
                            if len(p) > idx:
                                fo.write("\t".join([p[0], p[1], p[idx]]) + "\n")
                    st = mbc.rebuild_sample_vcf(raw_plain, gt_tsv, adj_plain, sm)
                    slog.result(f"重建 {sm}: 记录 {st['records']} GT已改 {st['gt_changed']}")
                    rc = self.runner.run(
                        f"{self.runner.tool('bcftools', 'bcftools view -Oz -o ' + self.runner.cpath(out_vcf) + ' ' + adj_plain)} "
                        f"&& {self.runner.tool('bcftools', 'bcftools index -t ' + self.runner.cpath(out_vcf))}",
                        logger=slog, outputs=[out_vcf])
                    if rc == 0:
                        adj_vcfs[sm] = out_vcf
            self.bdata["artifacts"]["per_sample_adjudicated_dir"] = \
                os.path.join(self.work, "per_sample_vcf")
            # per-sample PASS 裁决 VCF 空结果检查（RUN-34）：交付为空 → P1
            for sm, v in adj_vcfs.items():
                nrec = mbc.count_records(self.runner, v)
                self.metrics.setdefault("pass_records", {})[sm] = nrec
                if nrec == 0:
                    self.log.error(f"{sm} PASS 裁决 VCF 0 条记录（交付结果为空）")
            self.bdata["_adj_vcfs"] = adj_vcfs
        elif self.runner.dry_run:
            self.log.info("[DRY-RUN] 跳过裁决与每样本 VCF 重建的实际计算")

        # MultiQC（串行）：其他全部流程结束、qc/ 产物生成完全后，
        # 才生成最终汇总报告（fastqc/fastp/比对/去重/BQSR/覆盖度全量输入）
        if not self.runner.dry_run:
            mmultiqc.run_multiqc(self.runner, os.path.join(self.work, "qc"),
                                 os.path.join(self.work, "qc", "multiqc"), self.log)
            mq = sorted(glob.glob(os.path.join(self.work, "qc", "multiqc",
                                               "*multiqc_report.html")))
            if mq:
                self.bdata["artifacts"]["multiqc"] = mq[-1]

            # 交付导出：独立 Output/<批次>_<日期>/
            # （VCF+tbi+MultiQC 报告+md5sum+MANIFEST，幂等；裁决未产出则不导出）
            if self.bdata.get("_adj_vcfs"):
                delivery_dir = export_delivery(
                    self.runner, f"{self.batch}_{self.run_date}", self.bdata["_adj_vcfs"], self.log,
                    batch_note=f"｜批次 {self.batch}",
                    extra_files=(self.bdata["artifacts"]["multiqc"],)
                    if self.bdata["artifacts"].get("multiqc") else ())
                if delivery_dir:
                    self.bdata["artifacts"]["delivery_dir"] = delivery_dir
        self.step_time("step6")

        # Step6 里程碑：捕获效率 / 覆盖达标 / 结论口径 / NTC 污染
        cells = adj_stats.get("cells_total") if adj_stats else None
        call_rate = None
        if adj_stats and cells:
            call_rate = round((cells - adj_stats["kept_dotdot"]) * 100.0 / cells, 2)
        ntc_depth = (self.metrics.get("mosdepth_mean") or {}).get("NTC") \
            if self.excluded else None
        _step_notify(self.notify_on, self.batch, 6, "质量汇总", self.log,
                     total_steps=self.args.step,
                     samples=f"{len(self.merged) - len(self.failed)}/{len(self.merged)} 成功 | "
                             f"裁决 {len(self.bdata.get('_adj_vcfs') or {})} 样本",
                     # 均值排除对照（DEC-32）：NTC 深度近 0 会显著拉低批均值；
                     # 其自身深度由行尾"NTC 深度"单独播报、污染由 TH-21 负责判定
                     metrics=f"mean depth {_avg(self.metrics.get('mean_target_coverage'), self.excluded)}× | "
                             f"20X {_avg(self.metrics.get('pct_20x'), self.excluded)}% | "
                             f"捕获效率(selected) {_avg(self.metrics.get('pct_selected'), self.excluded)}% | "
                             f"call rate {call_rate}% | "
                             f"Ti/Tv PASS {self.cohort_stats.get('PASS', {}).get('titv')}"
                             + (f" | NTC 深度 {ntc_depth}×" if ntc_depth is not None else ""),
                     anomalies=_step_anoms(self, "step6")
                               + alerts.check_capture(
                                   self.metrics.get("mean_target_coverage"),
                                   self.metrics.get("pct_20x"),
                                   self.metrics.get("pct_selected"),
                                   excluded=self.excluded)
                               + alerts.check_variantqc(
                                   self.cohort_stats.get("PASS", {}).get("titv"), call_rate)
                               + alerts.check_ntc(ntc_depth)
                               + [(f"P1", f"{sm} PASS 裁决 VCF 0 条记录（交付结果为空）")
                                  for sm, v in sorted((self.metrics.get("pass_records") or {}).items())
                                  if v == 0]
                               + alerts.check_depth_cv(
                                   self.metrics.get("mean_target_coverage"),
                                   excluded=self.excluded),
                     artifacts=f"MultiQC + 裁决VCF "
                               + alerts.artifact_summary(os.path.join(
                                   self.work, "per_sample_vcf", "*.PASS.adjudicated.vcf.gz")),
                     log_hint=f"tail -f {self.work}/logs/sample_<样本>.self.log")



def process_batch(batch, batch_dir, args, plan, main_logger, run_date=""):
    """单批次全流程编排。返回批次结果 dict（status: success/failed/skipped）。
    输出目录 results/<批次>_<执行日期>/（多次运行隔离；同日重跑同目录幂等续跑）。
    各步骤实现拆分为 BatchCtx.step0~step6 方法（v2.7.0，原为单函数 ~700 行）"""
    work = args.out if args.out else config.batch_result_dir(batch, run_date)
    ctx = BatchCtx(batch, batch_dir, work, args, plan, main_logger)
    log = ctx.log
    ctx.bdata = {"status": "skipped", "batch_dir": batch_dir, "work_dir": work,
                 "samples": {"valid": [], "invalid": {}, "calling": []},
                 "steps": {}, "metrics": {}, "artifacts": {},
                 "error": None}
    ctx.notify_on = args.notify == "on" and not args.dry_run
    ctx.run_date = run_date
    bdata = ctx.bdata

    try:
        log.step(f"═══ 批次 {batch} 启动（{batch_dir} → {work}）═══")
        ctx.makedirs()

        early = ctx.step0_scan()          # 无有效样本 → 返回 bdata（批次 skipped）
        if early is not None:
            return early
        for upto, step_fn in ((1, ctx.step1_qc_trim), (2, ctx.step2_align),
                              (3, ctx.step3_markdup), (4, ctx.step4_bqsr),
                              (5, ctx.step5_variant_calling), (6, ctx.step6_summary_delivery)):
            if args.step >= upto:
                step_fn()
                ctx.disk_guard(f"step{upto}")   # P0 Step 间磁盘复查（RUN-34）

        # ── 汇总 ──
        merged, cohort_stats = ctx.merged, ctx.cohort_stats
        bdata["status"] = "failed" if ctx.failed and len(ctx.failed) >= len(merged) \
            else ("success" if not ctx.failed else "partial")
        bdata["failed_samples"] = dict(ctx.failed)
        bdata["samples"]["excluded"] = sorted(ctx.excluded)   # 对照清单入 run_summary（DEC-32）
        m = dict(ctx.metrics)
        if cohort_stats:
            m.update({
                "snp_raw": cohort_stats.get("raw", {}).get("snps"),
                "indel_raw": cohort_stats.get("raw", {}).get("indels"),
                "records_raw": cohort_stats.get("raw", {}).get("records"),
                "snp_pass": cohort_stats.get("PASS", {}).get("snps"),
                "indel_pass": cohort_stats.get("PASS", {}).get("indels"),
                "records_pass": cohort_stats.get("PASS", {}).get("records"),
                "titv_raw": cohort_stats.get("raw", {}).get("titv"),
                "titv_pass": cohort_stats.get("PASS", {}).get("titv"),
                "norm_split": cohort_stats.get("norm"),
            })
        bdata["metrics"] = m
        bdata["duration_s"] = round(time.time() - ctx.started, 1)
        if not ctx.runner.dry_run:
            report_mod.render_run_report(batch, bdata,
                                         os.path.join(work, "run_report.md"))

        # 成功通知
        if ctx.notify_on:
            _notify_result(batch, bdata, plan, log)
        log.step(f"═══ 批次 {batch} 结束: {bdata['status']}"
                 f"（{bdata['duration_s']}s）═══")
    except Exception as e:   # noqa: BLE001
        bdata["status"] = "failed"
        bdata["error"] = f"{e}\n{traceback.format_exc()}"
        log.error(f"批次 {batch} 失败: {e}")
        log.error(traceback.format_exc())
        bdata["duration_s"] = round(time.time() - ctx.started, 1)
        # ctx.notify_on 在进入 try 之前即已赋值——曾用局部 notify_on（仅全步骤
        # 成功后的汇总段才赋值），step 中途抛异常即 UnboundLocalError：失败通知
        # 发不出且异常炸穿 main()（RUN-42）
        if ctx.notify_on:
            tail = log.tail(20).replace("`", "'")
            dingtalk.notify(
                f"批次 {batch} 失败",
                f"#### 批次 {batch} 失败"
                f"\n\n> {e}"
                f"\n\n- **主机**：{os.uname().nodename}"
                f"\n\n- **时间**：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                f"\n\n**日志尾部 20 行**"
                f"\n\n> {tail}",
                logger=log)
    finally:
        ctx.factory.close_all()
        ctx.log.close()
    return bdata


def _pct(v):
    return f"{v * 100:.1f}%" if isinstance(v, (int, float)) else "?"


def _avg(d, excluded=()):
    """指标播报批均值：剔除对照样本（v2.22.0/DEC-32——DEC-31 只豁免了告警点名，
    "指标:"行均值此前仍含 NTC，其 reads 近 0 使 dup/深度/mapped 等均值失真）"""
    vals = [v for sm, v in (d or {}).items()
            if isinstance(v, (int, float)) and sm not in excluded]
    return round(sum(vals) / len(vals), 2) if vals else None


# ── 最终交付导出：独立 Output/<批次>/ 目录 ────────────────────────────
def export_delivery(runner, batch, adj_vcfs, log, batch_note="", extra_files=()):
    """把最终交付文件 *.PASS.adjudicated.vcf.gz(+.tbi) 与附加文件（MultiQC
    报告等）复制到独立交付目录，生成标准 md5sum.txt（可 md5sum -c 校验）、
    MANIFEST.tsv（含记录数/md5/来源）与交付说明 README.md。
    幂等：源未更新则不复制，manifest 每次重生成。"""
    import hashlib
    dbatch = config.batch_delivery_dir(batch)   # batch 参数已是"<批次>_<日期>"标签
    os.makedirs(dbatch, exist_ok=True)
    rows, md5_lines = [], []

    def _copy_idempotent(src, logname="复制"):
        dst = os.path.join(dbatch, os.path.basename(src))
        if not nonempty(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy2(src, dst)
            log.info(f"[交付] {logname} {os.path.basename(src)}")
        return dst

    def _md5(path):
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        return h.hexdigest()

    for sm in sorted(adj_vcfs):
        src = str(adj_vcfs[sm])
        if not nonempty(src):
            continue
        dst = _copy_idempotent(src)
        if nonempty(src + ".tbi"):
            _copy_idempotent(src + ".tbi", logname="复制索引")
        md5 = _md5(dst)
        records = mbc.count_records(runner, dst) if not runner.dry_run else None
        rows.append((sm, os.path.basename(dst), os.path.getsize(dst),
                     md5, records, src))
        md5_lines.append(f"{md5}  {os.path.basename(dst)}")
    n_vcf = len(rows)
    # 附加交付文件（MultiQC 报告等）：平铺进交付目录，同样幂等并纳入校验清单
    extra_names = []
    for src in extra_files:
        src = str(src)
        if not nonempty(src):
            continue
        dst = _copy_idempotent(src, logname="复制附件")
        kind = "multiqc" if "multiqc" in os.path.basename(src).lower() else "extra"
        md5 = _md5(dst)
        rows.append((kind, os.path.basename(dst), os.path.getsize(dst),
                     md5, "-", src))
        md5_lines.append(f"{md5}  {os.path.basename(dst)}")
        extra_names.append(os.path.basename(dst))
    if not rows:
        return None
    with open(os.path.join(dbatch, "md5sum.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(md5_lines) + "\n")
    with open(os.path.join(dbatch, "MANIFEST.tsv"), "w", encoding="utf-8") as f:
        f.write("sample\tfile\tsize_bytes\tmd5\trecords\tsource\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    readme = os.path.join(dbatch, "README.md")
    have = nonempty(readme)
    with open(readme, "w", encoding="utf-8") as f:
        f.write(f"# {batch} 批次最终交付\n\n"
                f"- 导出时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{batch_note}\n"
                f"- 交付文件：`<样本>.PASS.adjudicated.vcf.gz`（+ `.tbi` 索引）× {n_vcf}\n")
        if extra_names:
            f.write(f"- MultiQC 汇总报告：`{extra_names[0] if len(extra_names) == 1 else '…'}"
                    f"`（fastqc/fastp/比对/去重/BQSR/覆盖度等全流程 QC 汇总）\n")
        f.write(f"- 口径：GATK HaplotypeCaller（gVCF 联合分型）→ SNP/INDEL 硬过滤 → PASS →"
                f" 矩阵 `./.` 按 mosdepth 靶区深度 DP≥{config.DP_MIN} 裁决后重建的"
                f"每样本 VCF（仅替换 GT，其余字段原样保留）\n"
                f"- 校验：`md5sum -c md5sum.txt`；明细见 MANIFEST.tsv"
                f"（样本/文件/大小/md5/记录数/来源）\n"
                f"- 生成流程：GWAS pipeline v{config.PIPELINE_VERSION}"
                f"（{config.WORK_DIR}/pipeline，README 含完整口径与差异记录）\n")
    log.result(f"[交付] {dbatch}: VCF ×{n_vcf} + md5sum.txt + MANIFEST.tsv"
               + (f" + 附件 ×{len(extra_names)}" if extra_names else "")
               + ("（README 更新）" if have else ""))
    return dbatch


def write_delivery_index(delivered_dirs):
    """Output/INDEX.md 交付总索引：**累积合并**——扫描 Output/ 下全部批次交付
    目录 ∪ 本次运行的交付目录（多次运行隔离下历史交付不丢索引；RUN-29 前只写
    本次运行批次，跨日运行会把历史交付从索引中挤掉）。无任何目录时返回 None。"""
    dirs = set(delivered_dirs)
    if os.path.isdir(config.DELIVERY_DIR):
        dirs |= {os.path.join(config.DELIVERY_DIR, d)
                 for d in os.listdir(config.DELIVERY_DIR)
                 if os.path.isdir(os.path.join(config.DELIVERY_DIR, d))}
    if not dirs:
        return None
    idx = os.path.join(config.DELIVERY_DIR, "INDEX.md")
    os.makedirs(config.DELIVERY_DIR, exist_ok=True)
    with open(idx, "w", encoding="utf-8") as f:
        f.write(f"# 交付索引\n\n- 生成时间："
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
                f"| 批次_日期 | 交付目录 | VCF 数 |\n| --- | --- | --- |\n")
        for p in sorted(dirs, key=os.path.basename):
            n = len([x for x in os.listdir(p) if x.endswith(".vcf.gz")])
            f.write(f"| {os.path.basename(p)} | `{p}` | {n} |\n")
    return idx


def _step_anoms(ctx, step_prefix):
    """本步失败样本 → P0 异常行"""
    return [(f"P1", f"{sm} 执行失败：{reason.split(':', 1)[-1].strip()}")
            for sm, reason in ctx.failed.items() if reason.startswith(step_prefix)]


# ── cohort 复跑守卫（DEC-34/RUN-48；RUN-49 补漏）───────────────────────
def cohort_rerun_guard(runner, work, hc_ok, log, checks=()):
    """断点续跑时 HC 失败样本补回 → gvcf.list 样本集变化，但旧 cohort 链产物
    非空会被幂等 SKIP——旧口径（缺样本）的 cohort/矩阵/每样本 VCF 一路沿用：
    矩阵按旧列数裁决、补回样本 view -s 静默失败，交付残缺还报 success。
    守卫：读现存关键 VCF 的 header 样本清单（bcftools query -l，秒级）与本次
    hc_ok（sorted，与 VCF 列序同口径，DEC-05）比对，不一致即作废
    cohort/matrix/per_sample_vcf 全部派生产物（均可在链上重算，REQ-04），
    **连同 qc/bcftools_stats/cohort.*.stats 与 qc/multiqc/ 一并作废（RUN-49）**
    ——run_stats 按产物非空幂等、run_multiqc 按报告存在幂等，不删则陈旧
    统计（Ti/Tv 等）与陈旧 MultiQC 报告被 SKIP 沿用进通知与交付。
    dry-run / header 不可读 → 不判不作废（零副作用）。返回 True=已作废重算。"""
    stale = []
    for tag, vcf in checks:
        if nonempty(vcf):
            have = mbc.list_samples(runner, vcf)
            if have and have != list(hc_ok):
                stale.append(f"{tag}（{len(have)} 样本 ≠ 本次 {len(hc_ok)}）")
    if not stale:
        return False
    log.warn("cohort 样本集与本次联合分型名单不一致——" + "；".join(stale)
             + "。旧 cohort/矩阵/每样本 VCF/cohort 统计/MultiQC 作废重算"
               "（DEC-34：HC 失败样本补回后的断点续跑，防旧口径产物被幂等 SKIP 沿用；"
               "RUN-49 补漏 stats 与 MultiQC）")
    for d, pats in ((os.path.join(work, "cohort"), ("*.vcf.gz", "*.vcf.gz.tbi")),
                    (os.path.join(work, "matrix"), ("*.tsv",)),
                    (os.path.join(work, "per_sample_vcf"), ("*.vcf.gz", "*.vcf.gz.tbi")),
                    (os.path.join(work, "qc", "bcftools_stats"), ("cohort.*.stats",)),
                    (os.path.join(work, "qc", "multiqc"), ("*",))):
        for pat in pats:
            for p in glob.glob(os.path.join(d, pat)):
                if os.path.isdir(p):
                    shutil.rmtree(p)   # multiqc_data/ 等子目录
                else:
                    os.remove(p)
    return True


def _step_notify(notify_on, batch, step_no, name, log, total_steps=None, **kw):
    """步骤里程碑通知。total_steps=--step（编号步数：Step 0 为清点预备步不计入，
    全流程共 6 步；标题"（共 X 步）"，DEC-32，RUN-47 修正首版 +1 口径）；
    DINGTALK_MILESTONES=0 时只发异常级（P0/P1），OK 级静默。"""
    if not notify_on:
        return
    title, text = alerts.step_milestone(batch, step_no, name,
                                        total_steps=total_steps, **kw)
    level = alerts.worst_level(kw.get("anomalies"))
    if not config.DINGTALK_MILESTONES and level == "OK":
        return
    dingtalk.notify(title, text, logger=log)


def _notify_result(batch, bdata, plan, log):
    m = bdata.get("metrics", {})
    # 质量均值行排除对照（DEC-32）：mapped/dup/20X 等样本级指标对 NTC 无统计
    # 意义（reads 近 0 拉低批均值）；保留率含对照影响可忽略但同口径一并排除
    excluded = set(bdata.get("samples", {}).get("excluded") or ())
    final_anoms = []
    if bdata.get("failed_samples"):
        final_anoms += [("P1", f"{len(bdata['failed_samples'])} 个样本失败: "
                          + ", ".join(list(bdata['failed_samples'])[:5]))]
    lvl = alerts.worst_level(final_anoms) if final_anoms else \
        ("OK" if bdata['status'] == 'success' else "P1")
    dur = bdata.get("duration_s", 0)
    hh, mm, ss = dur // 3600, dur % 3600 // 60, dur % 60
    title = f"[GWAS][{lvl}] {batch}批次 · 全流程完成"
    text = (
        f"#### {title}"
        f"\n\n样本: {len(bdata['samples']['valid'])} 有效｜"
        f"联合分型 {len(bdata['samples'].get('calling') or [])}"
        f"\n\n指标: 总耗时 {int(hh)}h{int(mm)}m{int(ss)}s｜"
        f"workers 比对 {plan.workers}×{plan.threads}/GATK {plan.workers_gatk}/IO {plan.workers_io}"
        f"\n\n质量: fastp 保留率 {_avg(m.get('fastp_retention'), excluded)}%｜"
        f"mapped {_avg(m.get('mapped_pct'), excluded)}%｜dup {_avg(m.get('dup_pct'), excluded)}%｜"
        f"20X {_avg(m.get('pct_20x'), excluded)}%"
        f"\n\n变异: raw→PASS SNP {m.get('snp_raw')}→{m.get('snp_pass')}｜"
        f"INDEL {m.get('indel_raw')}→{m.get('indel_pass')}｜"
        f"Ti/Tv {m.get('titv_raw')}→{m.get('titv_pass')}")
    if final_anoms:
        for lv, msg in final_anoms:
            text += f"\n\n异常: [{lv}] {msg}" + {"P0": " ← 阻断级（中断分析）", "P1": " ← 需确认"}.get(lv, " ← 提示")
    else:
        text += "\n\n异常: 无"
    delivery = (bdata.get("artifacts") or {}).get("delivery_dir")
    if delivery:
        dnames = os.listdir(delivery)
        text += f"\n\n交付: {delivery}（*.PASS.adjudicated.vcf.gz ×"
        text += f"{len([f for f in dnames if f.endswith('.vcf.gz')])}"
        if any("multiqc_report.html" in f for f in dnames):
            text += " + MultiQC"
        text += " + md5sum.txt + MANIFEST.tsv）"
    text += (f"\n\n产物: {bdata.get('work_dir')}（cohort.PASS / 矩阵 / 每样本裁决VCF / "
             f"MultiQC）"
             f"\n\n日志: tail -f {bdata.get('work_dir')}/logs/pipeline_*.log")
    dingtalk.notify(title, text, logger=log)


def ignore_sighup():
    """启动即忽略 SIGHUP——代码级 nohup（DEC-33/RUN-48）：关闭启动命令所在的
    终端/SSH 会话时，内核向会话与前台进程组发 SIGHUP；nohup 只让本进程忽略，
    而 SIG_IGN 经 fork/exec 继承即可覆盖 sh/singularity/bcftools 等全部子进程，
    唯独 JVM 启动时会安装自己的 SIGHUP 处理器覆盖继承位（除非 -Xrs）——曾在
    Step 5 同瞬杀死 4 个 HC JVM（"Hangup"，exit=129）。GATK 命令已全部加 -Xrs
    （modules/gatk），fastqc 经 _JAVA_OPTIONS 注入；此处兜底主进程——未套
    nohup 的后台启动同样免疫。控制台输出由 TeeStream 兜底（写失败只落盘）。"""
    try:
        signal.signal(signal.SIGHUP, signal.SIG_IGN)
    except (AttributeError, ValueError, OSError):
        pass   # 平台无 SIGHUP / 非主线程：保持默认


def guard_host_python():
    """当前 Python 位于 singularity 容器内 → 启动即失败（快速失败原则）。
    本机 `python` 是容器别名（~/.bashrc: singularity exec ... mamba.sif python），
    容器内看不到宿主的 singularity——嵌套运行时所有子进程命令都会 127
    （实测 RUN-22/23）；本流程纯标准库，必须用宿主系统 Python 运行。"""
    if os.path.exists("/.singularity.d") or "SINGULARITY_CONTAINER" in os.environ \
            or "SINGULARITY_NAME" in os.environ:
        sys.exit(
            "[ERROR] 检测到当前 Python 运行在 singularity 容器内（alias python = 容器内解释器）。\n"
            "  本流程需在子进程中调用宿主 singularity，容器内嵌套运行会全部\n"
            "  `singularity: not found`(exit=127) 且时区/PATH 均为容器环境。\n"
            "  请改用宿主系统 Python：/usr/bin/python3 run_pipeline.py ...")


def main():
    guard_host_python()
    ignore_sighup()
    args = parse_args()

    if args.notify_test:
        # 企业机器人连通性（DEC-24）：markdown 链路 + 文件链路（媒体上传权限）
        ok, err = dingtalk.send_markdown(
            "[GWAS] 测试通知",
            "#### [GWAS] 流程测试通知\n > 这是 GWAS pipeline 的钉钉通知测试\n"
            "| 项目 | 值 |\n| --- | --- |\n| 状态 | OK |")
        file_note = ""
        if ok:
            with tempfile.TemporaryDirectory() as td:
                zp = os.path.join(td, "notify_test.zip")
                with zipfile.ZipFile(zp, "w") as zf:
                    zf.writestr("readme.txt", "GWAS pipeline 钉钉文件链路测试\n")
                fok, ferr = dingtalk.send_file(zp)
            file_note = "" if fok else f"；文件链路失败: {ferr}"
            ok = fok
        print(f"notify-test: {'发送成功（markdown+文件）' if ok else f'发送失败: {err}{file_note}'}")
        sys.exit(0 if ok else 1)

    # 先接管 stdio（缓冲模式）：资源计划表、快速失败报错、任何未捕获输出
    # 都被缓存；目标文件等批次发现后确定——results/ 根不允许散文件，
    # 运行日志归属批次结果目录（单批次=该批次；多批次=首个批次目录）
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    capture_stdio()

    # 资源规划（--serial 强制 workers=1；--workers 优先级最高）
    workers_override = 1 if args.serial else args.workers
    plan = resource.plan(profile=args.resource_profile,
                         workers_override=workers_override,
                         threads_override=args.threads,
                         max_memory_gb=args.max_memory)
    print(plan.table(), flush=True)

    # 主日志仅控制台输出：tee 已把控制台镜像写入 run_<ts>.log
    main_logger = Logger(None)
    main_logger.step(f"GWAS pipeline 启动: {' '.join(sys.argv)}")
    main_logger.info(f"dry-run={args.dry_run} notify={args.notify} "
                     f"profile={args.resource_profile} step≤{args.step}")
    main_logger.info(
        f"环境配置: {config.ENV_FILE}"
        + (f"（已加载 {len(config.ENV_FILE_KEYS)} 项）" if config.ENV_FILE_KEYS
           else "（不存在或为空——仅用进程环境变量+内置默认值，模板见 .env.example）"))
    if args.notify == "on" and not config.DINGTALK_WEBHOOK:
        main_logger.warn("钉钉通知开启但未配置 DINGTALK_WEBHOOK——通知将全部跳过"
                         f"（在 {config.ENV_FILE} 中配置，模板 .env.example）")

    # 批次发现：--batch 优先；否则遍历 --input 下全部批次子目录逐批独立执行
    # 相对路径按 $WORK 解析（如 --input 0_raw_data，无论从哪个 cwd 启动）
    # 相对 --input 按 RAW_DATA_DIR 语境解析（尊重 GWAS_RAW_DATA 覆盖——RUN-32
    # 事故：曾回落 $WORK/0_raw_data，验证脚本带着覆盖变量却扫到真实批次实跑）：
    # 同名即原始数据根本身，其余视为其子目录；绝不回落 $WORK
    input_root = args.input
    if not os.path.isdir(input_root):
        input_root = config.RAW_DATA_DIR \
            if os.path.basename(config.RAW_DATA_DIR) == args.input \
            else os.path.join(config.RAW_DATA_DIR, args.input)
    if args.batch:
        batches = [(args.batch, os.path.join(input_root, args.batch))]
        if not os.path.isdir(batches[0][1]):
            main_logger.error(f"批次目录不存在: {batches[0][1]}")
            sys.exit(2)
    else:
        if not os.path.isdir(input_root):
            main_logger.error(f"输入目录不存在: {input_root}")
            sys.exit(2)
        batches = sorted((d, os.path.join(input_root, d))
                         for d in os.listdir(input_root)
                         if os.path.isdir(os.path.join(input_root, d)))
    main_logger.info(f"待处理批次: {[b for b, _ in batches]}")

    # 运行日志落盘（dry-run / 无批次保持缓冲丢弃 → 零散文件不产生）
    run_date = ts.split("_")[0]   # 本次运行的执行日期（各批次共用，跨午夜不切换）
    run_log_path = None
    if batches and not args.dry_run:
        first_out = args.out if args.out \
            else config.batch_result_dir(batches[0][0], run_date)
        run_log_path = os.path.join(first_out, "logs", f"run_{ts}.log")
        tee_set_log_path(run_log_path)
    main_logger.info("运行日志（自动落盘，无需 shell 重定向）: "
                     + (run_log_path if run_log_path else "（dry-run 不落盘）"))
    if run_log_path:
        # 实跑取消控制台外显（DEC-09 v2.13.0）：上面的路径提示是控制台最后一行，
        # 此后全量日志只进 run_<ts>.log——nohup 后台不再向 nohup.out 倾倒；
        # dry-run/启动即退（无 run_log_path）不关，控制台照常可见
        tee_set_echo(False)

    summary = {
        "run": {"timestamp": ts, "argv": " ".join(sys.argv),
                "pipeline_version": config.PIPELINE_VERSION,
                "host": os.uname().nodename,
                "user": os.environ.get("USER", ""),
                "work_dir": config.WORK_DIR, "dry_run": args.dry_run,
                "step": args.step},
        "resource_plan": plan.to_dict(),
        "batches": {},
    }
    # run_summary.json 逐批写入该批次结果目录（快照含截至该批的运行信息，
    # 最后完成的批次文件即全貌）——results/ 根不产生任何散文件；
    # dry-run 不写（避免覆盖实跑的 run_summary）
    any_failed = False
    last_summary_path = None

    for batch, bdir in batches:
        main_logger.step(f"▶▶▶ 批次 {batch} 开始（输出 results/{batch}_{run_date}/）")
        bdata = process_batch(batch, bdir, args, plan, main_logger, run_date=run_date)
        bdata.pop("_adj_vcfs", None)
        summary["batches"][batch] = bdata
        out_dir = args.out if args.out else config.batch_result_dir(batch, run_date)
        last_summary_path = os.path.join(out_dir, "run_summary.json")
        if not args.dry_run:
            report_mod.write_run_summary(summary, last_summary_path)
        main_logger.result(f"批次 {batch}: {bdata['status']}")
        if bdata["status"] in ("failed", "partial"):
            any_failed = True

    # 交付总索引：Output/INDEX.md（累积合并全部历史交付目录；仅实跑交付后触发，
    # dry-run 零落盘不重写）
    delivered = {d.get("artifacts", {}).get("delivery_dir")
                 for d in summary["batches"].values()}
    delivered = {p for p in delivered if p}
    if delivered:
        idx = write_delivery_index(delivered)
        summary["delivery_index"] = idx
        if last_summary_path and not args.dry_run:
            report_mod.write_run_summary(summary, last_summary_path)

    # 交付文件打包推送钉钉（DEC-24，v2.15.0；分卷 DEC-27，v2.17.0）：
    # **批次全部结束后统一发送**——多批次逐批推会被通知淹没；每个 success 批次的
    # Output 交付目录 → 临时 zip（钉钉文件卡片仅支持 zip 等白名单后缀且 ≤20MB；
    # 超限自动分卷：每卷独立合法 zip ≤20MB，逐卷发卡片）→ 说明消息 + 文件卡片。
    # 发送失败只降级 WARN（通知链路问题不影响分析结果与退出码）
    if args.notify == "on" and not args.dry_run:
        for batch, d in sorted(summary["batches"].items()):
            ddir = d.get("artifacts", {}).get("delivery_dir")
            if d.get("status") != "success" or not ddir or not os.path.isdir(ddir):
                continue
            n_sm = len(d.get("samples", {}).get("calling")
                       or d.get("samples", {}).get("valid") or [])
            main_logger.step(f"▶ 交付 zip 推送钉钉: Output/{os.path.basename(ddir)}")
            dingtalk.send_zip_dir(
                ddir, f"[GWAS] {batch} 交付文件",
                f"#### 批次 {batch} 交付文件"
                f"\n\n> {len(d['samples']['valid'])} 样本｜{n_sm} 进入联合检测"
                f"\n\n- **目录**：`Output/{os.path.basename(ddir)}/`"
                f"\n\n- **内容**：每样本 `*.PASS.adjudicated.vcf.gz(+.tbi)`、"
                f"MultiQC 汇总报告、md5sum.txt（可 `md5sum -c` 校验）、MANIFEST.tsv"
                f"\n\n- **说明**：zip 为目录原样打包，解压后校验 md5 再取用",
                logger=main_logger)

    # 总通知
    if args.notify == "on" and not args.dry_run and len(batches) > 1:
        rows = "".join(
            f"\n\n{idx}. **{b}**：{d['status']}｜{d.get('duration_s', 0)}s｜"
            f"{len(d['samples']['valid'])} 样本"
            for idx, (b, d) in enumerate(sorted(summary["batches"].items()), 1))
        dingtalk.notify(
            "全部批次结束",
            f"#### 全部批次结束（{len(batches)} 批）"
            f"\n\n> 运行 {ts}｜汇总文件见下方"
            f"{rows}"
            f"\n\n- **run_summary**：`{last_summary_path}`",
            logger=main_logger)

    main_logger.step("运行结束: run_summary → "
                     + (last_summary_path if last_summary_path else "无（dry-run）"))
    main_logger.close()
    sys.exit(1 if any_failed else 0)


if __name__ == "__main__":
    main()
