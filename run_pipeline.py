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
import argparse
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
from logger import Logger, SampleLoggerFactory, capture_stdio, tee_set_log_path
from runner import Runner, nonempty
from modules import fastqc as mfastqc
from modules import fastp as mfastp
from modules import bwa_mem2 as mbwa
from modules import samtools as msam
from modules import gatk as mgatk
from modules import bcftools as mbc
from modules import mosdepth as mmos
from modules import multiqc as mmultiqc


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

    def slog(self, sm):
        if self.dry_run:
            return Logger(None, prefix=sm)
        return self.factory.get(sm, os.path.join(self.work, "logs"))

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


def process_batch(batch, batch_dir, args, plan, main_logger, run_date=""):
    """单批次全流程。返回批次结果 dict（status: success/failed/skipped）。
    输出目录 results/<批次>_<执行日期>/（多次运行隔离；同日重跑同目录幂等续跑）"""
    work = args.out if args.out else config.batch_result_dir(batch, run_date)
    ctx = BatchCtx(batch, batch_dir, work, args, plan, main_logger)
    log = ctx.log
    bdata = {"status": "skipped", "batch_dir": batch_dir, "work_dir": work,
             "samples": {"valid": [], "invalid": {}, "calling": []},
             "steps": {}, "metrics": {}, "artifacts": {},
             "error": None}
    notify_on = args.notify == "on" and not args.dry_run

    def step_time(name):
        bdata["steps"][name] = {"duration_s": round(time.time() - t0, 1)}

    try:
        log.step(f"═══ 批次 {batch} 启动（{batch_dir} → {work}）═══")
        ctx.makedirs()

        # ── Step 0：样本清点 + md5 + Lane 合并 + samples.tsv ──
        t0 = time.time()
        valid, invalid = scanner.scan_batch(batch_dir)
        if args.samples:
            keep = {s.strip() for s in args.samples.split(",") if s.strip()}
            for sm in list(valid):
                if sm not in keep:
                    invalid[sm] = "未在 --samples 白名单"
                    del valid[sm]
        bdata["samples"]["valid"] = sorted(valid)
        bdata["samples"]["invalid"] = invalid
        if not valid:
            log.warn(f"批次 {batch} 无任何有效样本，跳过该批次（退出码仍为 0）")
            if notify_on:
                reasons = "".join(f"\n\n- **{k}**：{v}" for k, v in invalid.items())
                dingtalk.notify(f"批次 {batch} 跳过",
                                f"#### 批次 {batch} 跳过（无有效样本）"
                                f"\n\n> 无效样本清单如下，退出码保持 0{reasons}",
                                logger=log)
            step_time("step0")
            return bdata
        for sm, reason in invalid.items():
            log.warn(f"无效样本（跳过）: {sm} —— {reason}")

        input_bytes = sum(os.path.getsize(f) for si in valid.values()
                          for f in si.r1 + si.r2)
        lanes_expected = sum(len(si.r1) for si in valid.values())
        n_initial = len(valid)
        log.info(f"有效样本 {n_initial} 个，输入体量 {_fmt_gb(input_bytes)}")

        # ── 开跑前检查（磁盘 / 依赖文件 / 样本名规范）──
        pre_anoms = []
        disk_path = work if os.path.isdir(work) else config.RESULTS_ROOT
        free_gb = shutil.disk_usage(disk_path).free / 1e9
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
        bad_names = [sm for sm in valid if not re.fullmatch(r"[A-Za-z0-9_.\-]+", sm)]
        if bad_names:
            pre_anoms.append(("P1", f"样本名含非常规字符（建议 A-Za-z0-9_.-）: {bad_names}"))
        for lv, msg in pre_anoms:
            (log.error if lv == "P0" else log.warn)(f"[开跑前-{lv}] {msg}")

        # md5 完整性校验（并行）；失败样本终止分析并进入通知
        md5_failed = set()
        if not args.dry_run:
            md5_failed, _ = scanner.verify_md5(batch_dir, log, workers=plan.workers)
            for sm in md5_failed & set(valid):
                log.error(f"md5 校验失败，样本 {sm} 终止分析")
                invalid[sm] = "md5 校验失败"
        for sm in md5_failed:
            valid.pop(sm, None)
        bdata["samples"]["valid"] = sorted(valid)

        # 启动通知（样本数/输入体量/资源计划 + 开跑前检查结论）
        if notify_on:
            p = plan
            lvl = alerts.worst_level(pre_anoms)
            body = (f"#### [GWAS][{lvl}] {batch}批次 · 启动"
                    f"\n\n样本: {len(valid)}/{n_initial} 有效｜输入 {_fmt_gb(input_bytes)}｜"
                    f"{datetime.now().strftime('%m-%d %H:%M')}"
                    f"\n\n指标: 探测 {p.cpu_detected} 线程 / {p.mem_detected} GB｜"
                    f"profile {p.profile}（预留 {p.reserve_cores} 核 + {p.reserve_mem_gb} GB）"
                    f"\n\n计划: workers {p.workers} × 每样本 {p.threads} 线程｜"
                    f"sort -m {p.sort_mem}｜GATK -Xmx {p.gatk_mem}｜cohort {p.cohort_mem}"
                    f"｜单样本峰值 {p.peak_per_sample_gb} GB")
            for lv, msg in pre_anoms:
                body += f"\n\n异常: [{lv}] {msg}" + (" ← 需确认" if lv == "P1" else " ← 阻断级")
            if not pre_anoms:
                body += "\n\n异常: 无"
            body += f"\n\n产物: {work}"
            dingtalk.notify(f"[GWAS][{lvl}] {batch}批次 · 启动", body, logger=log)
        if missing_deps:
            raise RuntimeError(f"开跑前检查 P0：依赖文件缺失 → {'; '.join(missing_deps)}")

        merged, merge_failed = scanner.merge_all(
            valid, os.path.join(work, "fastq_merged"), ctx.runner,
            ctx.slog, plan.workers)
        for sm, reason in merge_failed.items():
            invalid[sm] = reason
            valid.pop(sm, None)
            bdata["samples"]["invalid"][sm] = reason
        if not merged:
            raise RuntimeError(f"批次 {batch} 全部样本 Lane 合并失败")
        if not args.dry_run:
            scanner.write_samples_tsv(os.path.join(work, "samples.tsv"), merged, valid)
            bdata["artifacts"]["samples_tsv"] = os.path.join(work, "samples.tsv")
        bdata["samples"]["valid"] = sorted(merged)
        log.result(f"Step 0 完成: {len(merged)} 样本合并就绪，samples.tsv 已生成")
        lanes_merged = sum(len(si.r1) for sm, si in valid.items() if sm in merged)
        step_time("step0")
        _step_notify(notify_on, batch, 0, "清点与Lane合并", log,
                     samples=f"{len(merged)}/{n_initial} 成功 | Lane 合并 {lanes_merged}/{lanes_expected}",
                     metrics=f"输入 {_fmt_gb(input_bytes)} | md5 "
                             f"{'FAIL ' + str(len(md5_failed)) if md5_failed else 'OK'}"
                             f" | 无效样本 {len(invalid)}",
                     anomalies=[("P0", f"{sm} md5 校验失败（终止分析）") for sm in sorted(md5_failed)]
                               + [("P0", f"{sm} {rs}") for sm, rs in merge_failed.items()],
                     artifacts=f"fastq_merged/*_R*.fastq.gz "
                               + alerts.artifact_summary(os.path.join(work, "fastq_merged", "*_R*.fastq.gz")),
                     log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        excluded = {s.strip() for s in (args.exclude_samples or "").split(",") if s.strip()}

        # ── Step 1：FastQC(raw) → fastp → FastQC(trim) ──
        if args.step >= 1:
            t0 = time.time()
            log.step("Step 1: 原始 QC + 修剪")

            def s1(sm):
                slog = ctx.slog(sm)
                r1, r2 = merged[sm]
                qr1 = os.path.join(work, "qc", "fastqc_raw", f"{sm}_R1_fastqc.zip")
                qr2 = os.path.join(work, "qc", "fastqc_raw", f"{sm}_R2_fastqc.zip")
                if not mfastqc.run_fastqc(ctx.runner, [r1, r2],
                                          os.path.join(work, "qc", "fastqc_raw"),
                                          plan.fastqc_threads, slog):
                    return sm, False, "FastQC(raw) 失败"
                mfastqc.summarize([qr1, qr2], slog, "raw")
                c1 = os.path.join(work, "fastq_clean", f"{sm}_R1.fastq.gz")
                c2 = os.path.join(work, "fastq_clean", f"{sm}_R2.fastq.gz")
                fh = os.path.join(work, "qc", "fastp", f"{sm}.html")
                fj = os.path.join(work, "qc", "fastp", f"{sm}.json")
                if not mfastp.run_fastp(ctx.runner, sm, r1, r2, c1, c2, fh, fj,
                                        plan.fastp_threads, slog):
                    return sm, False, "fastp 失败"
                met = mfastp.parse_json(fj)
                slog.result(f"fastp: 保留率 {met.get('retention_pct')}% "
                            f"reads {met.get('before_reads')}→{met.get('after_reads')}")
                mfastp.check_retention(met, slog)
                t1 = os.path.join(work, "qc", "fastqc_trim", f"{sm}_R1_fastqc.zip")
                t2 = os.path.join(work, "qc", "fastqc_trim", f"{sm}_R2_fastqc.zip")
                if not mfastqc.run_fastqc(ctx.runner, [c1, c2],
                                          os.path.join(work, "qc", "fastqc_trim"),
                                          plan.fastqc_threads, slog):
                    slog.warn("FastQC(trim) 失败（不阻断）")
                else:
                    mfastqc.summarize([t1, t2], slog, "trim")
                    mfastqc.check_adapter_cleared(qr1, t1, slog)
                return sm, True, met

            res = _parallel({sm: (lambda sm=sm: s1(sm)) for sm in merged}, plan.workers)
            for sm, r in res.items():
                if isinstance(r, tuple) and len(r) == 3:
                    _, ok, note = r
                    if ok:
                        met = note or {}
                        ctx.metrics.setdefault("fastp_retention", {})[sm] = met.get("retention_pct")
                        # Q30 百分数口径（parse_json 已换算）；整体优先，回退 R1/R2 均值
                        q30 = met.get("q30_pct")
                        if q30 is None and met.get("q30_pct_r1") is not None:
                            q1 = met["q30_pct_r1"]
                            q2 = met.get("q30_pct_r2", q1)
                            q30 = round((q1 + q2) / 2, 2)
                        if q30 is not None:
                            ctx.metrics.setdefault("q30_pct", {})[sm] = q30
                    else:
                        ctx.failed[sm] = f"step1: {note}"
                        log.error(f"{sm} Step 1 失败: {note}")
                else:
                    ctx.failed[sm] = f"step1 异常: {r}"
                    log.error(f"{sm} Step 1 异常\n{r[1] if isinstance(r, tuple) else ''}")
            step_time("step1")
            log.result(f"Step 1 完成: 成功 {len(merged) - len(ctx.failed)}/{len(merged)}")
            _step_notify(notify_on, batch, 1, "QC+修剪", log,
                         samples=f"{len(merged) - len(ctx.failed)}/{len(merged)} 成功",
                         metrics=f"fastp 保留率 {_avg(ctx.metrics.get('fastp_retention'))}% | "
                                 f"Q30 {_avg(ctx.metrics.get('q30_pct'))}%",
                         anomalies=_step_anoms(ctx, "step1")
                                   + alerts.check_fastp(ctx.metrics.get("fastp_retention"),
                                                        ctx.metrics.get("q30_pct")),
                         artifacts=f"fastq_clean/*_R*.fastq.gz "
                                   + alerts.artifact_summary(os.path.join(work, "fastq_clean", "*_R*.fastq.gz")),
                         log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        # ── Step 2：比对 ──
        if args.step >= 2:
            t0 = time.time()
            log.step(f"Step 2: bwa-mem2 比对（{plan.workers} workers × {plan.threads} 线程，"
                     f"sort -m {plan.sort_mem}）")

            def s2(sm):
                slog = ctx.slog(sm)
                c1 = os.path.join(work, "fastq_clean", f"{sm}_R1.fastq.gz")
                c2 = os.path.join(work, "fastq_clean", f"{sm}_R2.fastq.gz")
                bam_dir = os.path.join(work, "bam", sm)
                if not ctx.runner.dry_run:
                    os.makedirs(bam_dir, exist_ok=True)
                sort_bam = os.path.join(bam_dir, f"{sm}.sort.bam")
                cmd = mbwa.build_align_pipe(ctx.runner, sm, c1, c2, sort_bam,
                                            plan.threads, plan.threads, plan.sort_mem)
                rc = ctx.runner.run(cmd, logger=slog, outputs=[sort_bam])
                if rc != 0:
                    return sm, False, "比对/排序失败"
                if not msam.run_index(ctx.runner, sort_bam, slog):
                    return sm, False, "index 失败"
                fl = os.path.join(work, "qc", "flagstat", f"{sm}.sorted.flagstat")
                st = os.path.join(work, "qc", "stats", f"{sm}.sorted.samtools.stats")
                msam.run_flagstat(ctx.runner, sort_bam, fl, slog)
                msam.run_stats(ctx.runner, sort_bam, st, slog)
                fs = msam.parse_flagstat(fl)
                slog.result(f"mapped {fs.get('mapped_pct')}% properly_paired "
                            f"{fs.get('pp_pct')}% singletons {fs.get('sgl_pct')}%")
                qc_ok = msam.qc_judgement(fs, slog)
                return sm, qc_ok, fs

            res = _parallel({sm: (lambda sm=sm: s2(sm)) for sm in merged
                             if sm not in ctx.failed}, plan.workers)
            for sm, r in res.items():
                if isinstance(r, tuple) and len(r) == 3:
                    _, ok, fs = r
                    if isinstance(fs, dict):
                        ctx.metrics.setdefault("mapped_pct", {})[sm] = fs.get("mapped_pct")
                        ctx.metrics.setdefault("pp_pct", {})[sm] = fs.get("pp_pct")
                    if not ok:
                        ctx.failed[sm] = "step2: 比对质检未达标"
                else:
                    ctx.failed[sm] = f"step2 异常: {r}"
                    log.error(f"{sm} Step 2 异常\n{r[1] if isinstance(r, tuple) else ''}")
            step_time("step2")
            log.result(f"Step 2 完成: 成功 {len(merged) - len(ctx.failed)}/{len(merged)}")
            _step_notify(notify_on, batch, 2, "比对", log,
                         samples=f"{len(merged) - len(ctx.failed)}/{len(merged)} 成功",
                         metrics=f"mapped {_avg(ctx.metrics.get('mapped_pct'))}% | "
                                 f"proper pair {_avg(ctx.metrics.get('pp_pct'))}%",
                         anomalies=_step_anoms(ctx, "step2")
                                   + alerts.check_flagstat(ctx.metrics.get("mapped_pct"),
                                                           ctx.metrics.get("pp_pct")),
                         artifacts=f"bam/*/*.sort.bam "
                                   + alerts.artifact_summary(os.path.join(work, "bam", "*", "*.sort.bam")),
                         log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        # ── Step 3：MarkDuplicates ──
        if args.step >= 3:
            t0 = time.time()
            log.step("Step 3: MarkDuplicates 去重")

            def s3(sm):
                slog = ctx.slog(sm)
                bam_dir = os.path.join(work, "bam", sm)
                sort_bam = os.path.join(bam_dir, f"{sm}.sort.bam")
                md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
                metrics_f = os.path.join(bam_dir, f"{sm}.markdup.metrics")
                if not mgatk.markdup(ctx.runner, sort_bam, md_bam, metrics_f,
                                     plan.gatk_mem, slog):
                    return sm, False, None
                met = mgatk.parse_markdup_metrics(metrics_f)
                slog.result(f"READ_PAIRS={met.get('READ_PAIRS_EXAMINED')} "
                            f"DUP={_pct(met.get('PERCENT_DUPLICATION'))} "
                            f"ELS={met.get('ESTIMATED_LIBRARY_SIZE')}")
                mgatk.qc_duplication(met, slog)
                fl = os.path.join(work, "qc", "flagstat", f"{sm}.markdup.flagstat")
                st = os.path.join(work, "qc", "stats", f"{sm}.markdup.samtools.stats")
                msam.run_flagstat(ctx.runner, md_bam, fl, slog)
                msam.run_stats(ctx.runner, md_bam, st, slog)
                # 去重前 flagstat 的 duplicates 行无意义，不采集（硬性要求）
                return sm, True, met

            res = _parallel({sm: (lambda sm=sm: s3(sm)) for sm in merged
                             if sm not in ctx.failed}, plan.workers)
            for sm, r in res.items():
                if isinstance(r, tuple) and len(r) == 3:
                    _, ok, met = r
                    ctx.metrics.setdefault("dup_pct", {})[sm] = \
                        round(met["PERCENT_DUPLICATION"] * 100, 2) \
                        if met and met.get("PERCENT_DUPLICATION") is not None else None
                    if met:
                        ctx.metrics.setdefault("els", {})[sm] = \
                            met.get("ESTIMATED_LIBRARY_SIZE")
                    if not ok:
                        ctx.failed[sm] = "step3: MarkDuplicates 失败"
                else:
                    ctx.failed[sm] = f"step3 异常: {r}"
            step_time("step3")
            log.result(f"Step 3 完成: 成功 {len(merged) - len(ctx.failed)}/{len(merged)}")
            _step_notify(notify_on, batch, 3, "去重", log,
                         samples=f"{len(merged) - len(ctx.failed)}/{len(merged)} 成功",
                         metrics=f"重复率 {_avg(ctx.metrics.get('dup_pct'))}% | "
                                 f"ELS 最小 {alerts._min_item(ctx.metrics.get('els'))[1]}",
                         anomalies=_step_anoms(ctx, "step3")
                                   + alerts.check_dup(ctx.metrics.get("dup_pct")),
                         artifacts=f"bam/*/*.markdup.bam "
                                   + alerts.artifact_summary(os.path.join(work, "bam", "*", "*.markdup.bam")),
                         log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        # ── Step 4：BQSR ──
        if args.step >= 4:
            t0 = time.time()
            log.step("Step 4: BQSR 校准")

            def s4(sm):
                slog = ctx.slog(sm)
                bam_dir = os.path.join(work, "bam", sm)
                md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
                table = os.path.join(bam_dir, f"{sm}.recal.table")
                bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
                if not mgatk.base_recalibrator(ctx.runner, md_bam, table,
                                               plan.gatk_mem, slog):
                    return sm, False, "BaseRecalibrator 失败"
                if not mgatk.apply_bqsr(ctx.runner, md_bam, table, bqsr_bam,
                                        plan.gatk_mem, slog):
                    return sm, False, "ApplyBQSR 失败"
                fl_md = os.path.join(work, "qc", "flagstat", f"{sm}.markdup.flagstat")
                fl_rc = os.path.join(work, "qc", "flagstat", f"{sm}.recal.flagstat")
                st_rc = os.path.join(work, "qc", "stats", f"{sm}.recal.samtools.stats")
                msam.run_flagstat(ctx.runner, bqsr_bam, fl_rc, slog)
                msam.run_stats(ctx.runner, bqsr_bam, st_rc, slog)
                if not ctx.runner.dry_run and not msam.flagstat_identical(fl_md, fl_rc):
                    return sm, False, "markdup 与 BQSR flagstat 不一致（判 FAIL）"
                slog.result("BQSR 前后 flagstat 逐行一致 ✓" if not ctx.runner.dry_run
                            else "BQSR 前后 flagstat 断言（dry-run 跳过实际比对）")
                return sm, True, None

            res = _parallel({sm: (lambda sm=sm: s4(sm)) for sm in merged
                             if sm not in ctx.failed}, plan.workers)
            for sm, r in res.items():
                if isinstance(r, tuple) and len(r) == 3:
                    _, ok, note = r
                    if not ok:
                        ctx.failed[sm] = f"step4: {note}"
                        log.error(f"{sm}: {note}")
                else:
                    ctx.failed[sm] = f"step4 异常: {r}"
            step_time("step4")
            n_bqsr_ok = len([sm for sm in merged if sm not in ctx.failed])
            log.result(f"Step 4 完成: 成功 {n_bqsr_ok}/{len(merged)}")
            _step_notify(notify_on, batch, 4, "BQSR校准", log,
                         samples=f"{n_bqsr_ok}/{len(merged)} 成功",
                         metrics="markdup↔BQSR flagstat 逐行一致断言 "
                                 f"{n_bqsr_ok}/{n_bqsr_ok} 通过",
                         anomalies=_step_anoms(ctx, "step4"),
                         artifacts=f"bam/*/*.markdup.BQSR.bam "
                                   + alerts.artifact_summary(os.path.join(work, "bam", "*", "*.markdup.BQSR.bam")),
                         log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        # ── Step 5：变异检测（gVCF → 联合分型 → 硬过滤 → VCF） ──
        cohort_stats = {}
        if args.step >= 5:
            t0 = time.time()
            log.step("Step 5: 变异检测")
            # 排序保证与 VCF 样本列序（gvcf.list=sorted）一致——矩阵列映射/裁决/重建都依赖此顺序
            calling = sorted(sm for sm in merged
                             if sm not in ctx.failed and sm not in excluded)
            for sm in sorted(excluded & set(merged)):
                log.info(f"对照样本 {sm} 按配置排除出联合变异检测（纳入 QC）")
            bdata["samples"]["calling"] = calling
            if not calling:
                raise RuntimeError("无可用于联合变异检测的样本（全部失败或被排除）")
            if not mgatk.prep_interval_list(ctx.runner, log):
                raise RuntimeError("interval_list 准备失败")

            def s5(sm):
                slog = ctx.slog(sm)
                bam_dir = os.path.join(work, "bam", sm)
                bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
                gvcf = os.path.join(work, "gvcf", f"{sm}.g.vcf.gz")
                if not mgatk.haplotypecaller(ctx.runner, bqsr_bam, gvcf,
                                             plan.gatk_mem, plan.hc_hmm_threads, slog):
                    return sm, False
                return sm, True

            res = _parallel({sm: (lambda sm=sm: s5(sm)) for sm in calling}, plan.workers)
            hc_ok = []
            for sm, r in res.items():
                if r == (sm, True):
                    hc_ok.append(sm)
                else:
                    ctx.failed[sm] = "step5: HaplotypeCaller 失败"
                    log.error(f"{sm} HC 失败\n{r[1] if isinstance(r, tuple) else ''}")
            if not hc_ok:
                raise RuntimeError("全部 HaplotypeCaller 失败")
            hc_ok = sorted(hc_ok)
            bdata["samples"]["calling"] = hc_ok   # 实际进入联合分型的样本（排序）
            log.result(f"HaplotypeCaller 完成 {len(hc_ok)}/{len(calling)}，进入联合分型")

            # cohort 级（串行，gvcf.list 每次重新生成）
            coh = os.path.join(work, "cohort")
            gvcf_list = os.path.join(coh, "gvcf.list")
            with open(gvcf_list, "w", encoding="utf-8") as f:
                for sm in sorted(hc_ok):
                    f.write(ctx.runner.cpath(os.path.join(work, "gvcf", f"{sm}.g.vcf.gz")) + "\n")
            log.info(f"gvcf.list 重新生成: {len(hc_ok)} 个 gVCF（批次内联合，严禁跨批次）")
            combined = os.path.join(coh, "cohort.g.vcf.gz")
            raw_vcf = os.path.join(coh, "cohort.raw.vcf.gz")
            if not mgatk.combine_gvcfs(ctx.runner, gvcf_list, combined,
                                       plan.cohort_mem, log):
                raise RuntimeError("CombineGVCFs 失败")
            if not mgatk.genotype_gvcfs(ctx.runner, combined, raw_vcf,
                                        plan.cohort_mem, log):
                raise RuntimeError("GenotypeGVCFs 失败")
            mbc.run_stats(ctx.runner, raw_vcf,
                          os.path.join(work, "qc", "bcftools_stats", "cohort.raw.stats"), log)
            cohort_stats["raw"] = mbc.parse_stats(
                os.path.join(work, "qc", "bcftools_stats", "cohort.raw.stats"))

            split_vcf = os.path.join(coh, "cohort.raw.split.vcf.gz")
            norm_stats_json = os.path.join(coh, "norm_stats.json")
            ok, norm_stats = mbc.norm_split(ctx.runner, raw_vcf, split_vcf, log)
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
            cohort_stats["norm"] = norm_stats
            log.result(f"norm 摊平: {norm_stats}（后续一切靶区提取加 ±100bp padding）")

            snp_v = os.path.join(coh, "cohort.snp.vcf.gz")
            indel_v = os.path.join(coh, "cohort.indel.vcf.gz")
            snp_f = os.path.join(coh, "cohort.snp.hardfiltered.vcf.gz")
            indel_f = os.path.join(coh, "cohort.indel.hardfiltered.vcf.gz")
            hard_v = os.path.join(coh, "cohort.hardfiltered.vcf.gz")
            pass_v = os.path.join(coh, "cohort.PASS.vcf.gz")
            if not mgatk.select_variants(ctx.runner, split_vcf, "SNP", snp_v,
                                         plan.gatk_mem, log) \
                    or not mgatk.select_variants(ctx.runner, split_vcf, "INDEL", indel_v,
                                                 plan.gatk_mem, log):
                raise RuntimeError("SelectVariants 失败")
            if not mgatk.variant_filtration(ctx.runner, snp_v, snp_f,
                                            config.SNP_HARD_FILTERS, plan.gatk_mem, log) \
                    or not mgatk.variant_filtration(ctx.runner, indel_v, indel_f,
                                                    config.INDEL_HARD_FILTERS,
                                                    plan.gatk_mem, log):
                raise RuntimeError("VariantFiltration 失败")
            if not mbc.concat(ctx.runner, [snp_f, indel_f], hard_v, log):
                raise RuntimeError("bcftools concat 失败")
            mbc.index_tbi(ctx.runner, hard_v, log)
            filt_ok, filt_dist = mbc.filter_column_check(ctx.runner, hard_v, log)
            if not filt_ok:
                raise RuntimeError("FILTER 列出现 '.'——过滤漏跑，判 FAIL")
            log.result(f"FILTER 标签分布: {filt_dist}")
            if not mbc.view_pass(ctx.runner, hard_v, pass_v, log):
                raise RuntimeError("PASS 提取失败")
            mbc.index_tbi(ctx.runner, pass_v, log)

            # 双口径导出 + stats
            mbc.export_genotype_matrix(
                ctx.runner, hard_v, os.path.join(work, "matrix", "genotype_matrix.tsv"), log)
            mbc.export_detail_pass(
                ctx.runner, pass_v, os.path.join(work, "matrix",
                                                 "genotype_detail_PASS.tsv"), log)
            mbc.run_stats(ctx.runner, hard_v, os.path.join(
                work, "qc", "bcftools_stats", "cohort.hardfiltered.stats"), log)
            mbc.run_stats(ctx.runner, pass_v, os.path.join(
                work, "qc", "bcftools_stats", "cohort.PASS.stats"), log)
            cohort_stats["hardfiltered"] = mbc.parse_stats(os.path.join(
                work, "qc", "bcftools_stats", "cohort.hardfiltered.stats"))
            cohort_stats["PASS"] = mbc.parse_stats(os.path.join(
                work, "qc", "bcftools_stats", "cohort.PASS.stats"))
            # 每样本 hardfiltered / PASS（未裁决版，供追溯）
            for sm in calling:
                hf_sm = os.path.join(work, "per_sample_vcf", f"{sm}.hardfiltered.vcf.gz")
                ps_sm = os.path.join(work, "per_sample_vcf", f"{sm}.PASS.vcf.gz")
                mbc.split_sample(ctx.runner, hard_v, sm, hf_sm, ctx.slog(sm))
                mbc.index_tbi(ctx.runner, hf_sm, ctx.slog(sm))
                mbc.split_sample(ctx.runner, pass_v, sm, ps_sm, ctx.slog(sm))
                mbc.index_tbi(ctx.runner, ps_sm, ctx.slog(sm))

            tt_raw = cohort_stats["raw"].get("titv")
            tt_pass = cohort_stats["PASS"].get("titv")
            log.result(f"Ti/Tv: raw {tt_raw} → PASS {tt_pass}（应上升；panel 参考区间 2.5-3.5，"
                       f"看趋势不看绝对值）")
            step_time("step5")
            bdata["artifacts"].update({
                "cohort_PASS_vcf": pass_v,
                "genotype_matrix": os.path.join(work, "matrix", "genotype_matrix.tsv"),
                "genotype_detail_PASS": os.path.join(work, "matrix",
                                                     "genotype_detail_PASS.tsv")})

            # 对账·数量：矩阵行数 vs hardfiltered 限定 targets 的记录数（应一致）
            recon_anoms = _step_anoms(ctx, "step5")
            if not ctx.runner.dry_run:
                matrix_n = sum(1 for _ in open(
                    os.path.join(work, "matrix", "genotype_matrix.tsv"), encoding="utf-8"))
                # 同口径核对（矩阵由 query -R 生成，核对也用 -R：
                # -R/-T 在区间边界对跨界 indel 的取舍不同，混用会有固有差额）
                view_n_out = ctx.runner.out(
                    f"{ctx.runner.tool('bcftools', 'bcftools view -R ' + ctx.runner.cpath(config.TARGETS_SORTED_BED) + ' -H ' + ctx.runner.cpath(hard_v))} | wc -l")
                try:
                    view_n = int(view_n_out.strip().split()[-1])
                except (ValueError, IndexError):
                    view_n = None
                log.result(f"对账·数量: 矩阵 {matrix_n} 行 vs view -T 计数 {view_n}")
                if view_n is not None and matrix_n != view_n:
                    recon_anoms.append(("P1", f"对账·数量：矩阵行数 {matrix_n} ≠ 靶区记录数 "
                                              f"{view_n}（导出可能不完整）"))
                # 对账·新鲜度：关键 VCF mtime < 本次启动（断点续跑复用旧产物）
                stale = []
                for key_vcf in [pass_v, hard_v]:
                    if nonempty(key_vcf) and os.path.getmtime(key_vcf) < ctx.started:
                        stale.append(os.path.basename(key_vcf))
                adj_all = sorted(glob.glob(os.path.join(
                    work, "per_sample_vcf", "*.PASS.adjudicated.vcf.gz")))
                if adj_all and all(os.path.getmtime(p) < ctx.started for p in adj_all):
                    stale.append("每样本裁决VCF×" + str(len(adj_all)))
                if stale:
                    recon_anoms.append(("P1", "对账·新鲜度：" + ", ".join(stale)
                                        + " 为历史运行产物（断点续跑复用）"))

            _step_notify(notify_on, batch, 5, "变异检测", log,
                         samples=f"gVCF {len(hc_ok)}/{len(calling)} | 联合分型样本 {len(hc_ok)}",
                         metrics=f"raw {cohort_stats['raw'].get('records')} → PASS "
                                 f"{cohort_stats['PASS'].get('records')} | "
                                 f"SNP/INDEL PASS {cohort_stats['PASS'].get('snps')}/"
                                 f"{cohort_stats['PASS'].get('indels')} | "
                                 f"Ti/Tv {tt_raw}→{tt_pass} | "
                                 f"norm split/realigned {norm_stats.get('split')}/"
                                 f"{norm_stats.get('realigned')}",
                         anomalies=recon_anoms,
                         artifacts=f"cohort/cohort.PASS.vcf.gz "
                                   + alerts.artifact_summary(os.path.join(work, "cohort", "*.vcf.gz")),
                         log_hint=f"tail -f {work}/logs/pipeline_*.log")

        # ── Step 6：测序质量与汇总 ──
        if args.step >= 6:
            t0 = time.time()
            log.step("Step 6: mosdepth ×2 + HsMetrics + 矩阵裁决 + MultiQC")

            def s6(sm):
                slog = ctx.slog(sm)
                bam_dir = os.path.join(work, "bam", sm)
                md_bam = os.path.join(bam_dir, f"{sm}.markdup.bam")
                bqsr_bam = os.path.join(bam_dir, f"{sm}.markdup.BQSR.bam")
                ms = {}
                for tag, bam in (("md", md_bam), ("bqsr", bqsr_bam)):
                    prefix = os.path.join(work, "qc", "mosdepth", f"{sm}.{tag}")
                    if not mmos.run_mosdepth(ctx.runner, bam, prefix,
                                             plan.mosdepth_threads, slog):
                        slog.warn(f"mosdepth({tag}) 失败（不阻断）")
                        continue
                    summ = mmos.parse_summary(prefix)
                    ms[tag] = summ["mean"]
                    slog.result(f"mosdepth[{tag}] 靶区均值={ms[tag]}×")
                hs_txt = os.path.join(work, "qc", "hsmetrics", f"{sm}.hs_metrics.txt")
                if not mgatk.collect_hsmetrics(ctx.runner, bqsr_bam, hs_txt,
                                               plan.gatk_mem, slog):
                    return sm, False, (ms, None)
                hs = mgatk.parse_hsmetrics(hs_txt)
                mgatk.qc_hsmetrics(hs, slog, control=sm in excluded)
                return sm, True, (ms, hs)

            res = _parallel({sm: (lambda sm=sm: s6(sm)) for sm in merged
                             if sm not in ctx.failed}, plan.workers)
            for sm, r in res.items():
                if isinstance(r, tuple) and len(r) == 3:
                    _, ok, (ms, hs) = r
                    ctx.metrics.setdefault("mosdepth_mean", {})[sm] = ms.get("bqsr")
                    if hs:
                        ctx.metrics.setdefault("mean_target_coverage", {})[sm] = \
                            hs.get("MEAN_TARGET_COVERAGE")
                        ctx.metrics.setdefault("pct_20x", {})[sm] = \
                            round((hs.get("PCT_TARGET_BASES_20X") or 0) * 100, 2)
                        ctx.metrics.setdefault("on_target_pct", {})[sm] = \
                            hs.get("ON_TARGET_PCT")
                    if not ok:
                        ctx.failed[sm] = "step6: HsMetrics 失败"
                else:
                    ctx.failed[sm] = f"step6 异常: {r}"

            # 矩阵 ./. 裁决（mosdepth bqsr regions，DP≥20 改判 0/0）+ 每样本 PASS 重建
            calling = sorted(bdata["samples"].get("calling") or
                             [sm for sm in merged if sm not in ctx.failed])
            adj_tsv = os.path.join(work, "matrix", "genotype_matrix.adjudicated.tsv")
            adj_stats = None
            if calling and cohort_stats and not ctx.runner.dry_run:
                regions_cache = {}
                for sm in calling:
                    prefix = os.path.join(work, "qc", "mosdepth", f"{sm}.bqsr")
                    regions_cache[sm] = mmos.load_regions(prefix)

                def region_of(sm, chrom, pos):
                    return mmos.region_depth(regions_cache.get(sm, []), chrom, pos)

                matrix = os.path.join(work, "matrix", "genotype_matrix.tsv")
                out_lines, adj_stats = mbc.adjudicate_matrix(matrix, region_of, calling,
                                                         config.DP_MIN)
                with open(adj_tsv, "w", encoding="utf-8") as f:
                    f.write("\n".join(out_lines) + "\n")
                log.result(f"矩阵 ./. 裁决: 总 ./. {adj_stats['dotdot_total']} → "
                           f"改判 0/0 {adj_stats['filled_00']} · 保留 ./. "
                           f"{adj_stats['kept_dotdot']}（DP_MIN={config.DP_MIN}）")
                bdata["artifacts"]["genotype_matrix_adjudicated"] = adj_tsv

                # 每样本 PASS VCF 重建（bcftools view -s 拆分 + GT 替换，其余字段原样保留）
                pass_v = os.path.join(work, "cohort", "cohort.PASS.vcf.gz")
                adj_vcfs = {}
                import tempfile
                import shlex as _shlex
                for sm in calling:
                    slog = ctx.slog(sm)
                    out_vcf = os.path.join(work, "per_sample_vcf",
                                           f"{sm}.PASS.adjudicated.vcf.gz")
                    if nonempty(out_vcf) and nonempty(out_vcf + ".tbi"):
                        adj_vcfs[sm] = out_vcf
                        continue
                    with tempfile.TemporaryDirectory() as td:
                        raw_plain = os.path.join(td, f"{sm}.raw.vcf")
                        gt_tsv = os.path.join(td, f"{sm}.gt.tsv")
                        adj_plain = os.path.join(td, f"{sm}.adj.vcf")
                        rc = ctx.runner.run(
                            f"{ctx.runner.tool('bcftools', 'bcftools view -s ' + _shlex.quote(sm) + ' --min-ac 0 ' + ctx.runner.cpath(pass_v))} > {raw_plain}",
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
                        rc = ctx.runner.run(
                            f"{ctx.runner.tool('bcftools', 'bcftools view -Oz -o ' + ctx.runner.cpath(out_vcf) + ' ' + adj_plain)} "
                            f"&& {ctx.runner.tool('bcftools', 'bcftools index -t ' + ctx.runner.cpath(out_vcf))}",
                            logger=slog, outputs=[out_vcf])
                        if rc == 0:
                            adj_vcfs[sm] = out_vcf
                bdata["artifacts"]["per_sample_adjudicated_dir"] = \
                    os.path.join(work, "per_sample_vcf")
                bdata["_adj_vcfs"] = adj_vcfs
                # 最终交付导出：独立 delivery/<批次>_<日期>/（VCF+md5sum+MANIFEST，幂等）
                delivery_dir = export_delivery(ctx.runner, f"{batch}_{run_date}",
                                               adj_vcfs, log,
                                               batch_note=f"｜批次 {batch}")
                if delivery_dir:
                    bdata["artifacts"]["delivery_dir"] = delivery_dir
            elif ctx.runner.dry_run:
                log.info("[DRY-RUN] 跳过裁决与每样本 VCF 重建的实际计算")

            # MultiQC（串行）
            if not ctx.runner.dry_run:
                mmultiqc.run_multiqc(ctx.runner, os.path.join(work, "qc"),
                                     os.path.join(work, "qc", "multiqc"), log)
                mq = sorted(glob.glob(os.path.join(work, "qc", "multiqc",
                                                   "*multiqc_report.html")))
                if mq:
                    bdata["artifacts"]["multiqc"] = mq[-1]
            step_time("step6")

            # Step6 里程碑：捕获效率 / 覆盖达标 / 结论口径 / NTC 污染
            cells = adj_stats.get("cells_total") if adj_stats else None
            call_rate = None
            if adj_stats and cells:
                call_rate = round((cells - adj_stats["kept_dotdot"]) * 100.0 / cells, 2)
            ntc_depth = (ctx.metrics.get("mosdepth_mean") or {}).get("NTC") \
                if excluded else None
            _step_notify(notify_on, batch, 6, "质量汇总", log,
                         samples=f"{len(merged) - len(ctx.failed)}/{len(merged)} 成功 | "
                                 f"裁决 {len(bdata.get('_adj_vcfs') or {})} 样本",
                         metrics=f"mean depth {_avg(ctx.metrics.get('mean_target_coverage'))}× | "
                                 f"20X {_avg(ctx.metrics.get('pct_20x'))}% | "
                                 f"on-target {_avg(ctx.metrics.get('on_target_pct'))}% | "
                                 f"call rate {call_rate}% | "
                                 f"Ti/Tv PASS {cohort_stats.get('PASS', {}).get('titv')}"
                                 + (f" | NTC 深度 {ntc_depth}×" if ntc_depth is not None else ""),
                         anomalies=_step_anoms(ctx, "step6")
                                   + alerts.check_capture(
                                       ctx.metrics.get("mean_target_coverage"),
                                       ctx.metrics.get("pct_20x"),
                                       ctx.metrics.get("on_target_pct"))
                                   + alerts.check_variantqc(
                                       cohort_stats.get("PASS", {}).get("titv"), call_rate)
                                   + alerts.check_ntc(ntc_depth),
                         artifacts=f"MultiQC + 裁决VCF "
                                   + alerts.artifact_summary(os.path.join(
                                       work, "per_sample_vcf", "*.PASS.adjudicated.vcf.gz")),
                         log_hint=f"tail -f {work}/logs/sample_<样本>.log")

        # ── 汇总 ──
        bdata["status"] = "failed" if ctx.failed and len(ctx.failed) >= len(merged) \
            else ("success" if not ctx.failed else "partial")
        bdata["failed_samples"] = dict(ctx.failed)
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
        if notify_on:
            _notify_result(batch, bdata, plan, log)
        log.step(f"═══ 批次 {batch} 结束: {bdata['status']}"
                 f"（{bdata['duration_s']}s）═══")
    except Exception as e:   # noqa: BLE001
        bdata["status"] = "failed"
        bdata["error"] = f"{e}\n{traceback.format_exc()}"
        log.error(f"批次 {batch} 失败: {e}")
        log.error(traceback.format_exc())
        bdata["duration_s"] = round(time.time() - ctx.started, 1)
        if notify_on:
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


def _avg(d):
    vals = [v for v in (d or {}).values() if isinstance(v, (int, float))]
    return round(sum(vals) / len(vals), 2) if vals else None


# ── 最终交付导出：独立 delivery/<批次>/ 目录 ────────────────────────────
def export_delivery(runner, batch, adj_vcfs, log, batch_note=""):
    """把最终交付文件 *.PASS.adjudicated.vcf.gz(+.tbi) 复制到独立交付目录，
    生成标准 md5sum.txt（可 md5sum -c 校验）、MANIFEST.tsv（含记录数/md5/来源）
    与交付说明 README.md。幂等：源未更新则不复制，manifest 每次重生成。"""
    import hashlib
    dbatch = config.batch_delivery_dir(batch)   # batch 参数已是"<批次>_<日期>"标签
    os.makedirs(dbatch, exist_ok=True)
    rows, md5_lines = [], []
    for sm in sorted(adj_vcfs):
        src = str(adj_vcfs[sm])
        if not nonempty(src):
            continue
        dst = os.path.join(dbatch, os.path.basename(src))
        if not nonempty(dst) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy2(src, dst)
            log.info(f"[交付] 复制 {os.path.basename(src)}")
        if nonempty(src + ".tbi"):
            dst_tbi = dst + ".tbi"
            if not nonempty(dst_tbi) or os.path.getmtime(src + ".tbi") > os.path.getmtime(dst_tbi):
                shutil.copy2(src + ".tbi", dst_tbi)
        h = hashlib.md5()
        with open(dst, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        md5 = h.hexdigest()
        records = mbc.count_records(runner, dst) if not runner.dry_run else None
        rows.append((sm, os.path.basename(dst), os.path.getsize(dst), md5, records, src))
        md5_lines.append(f"{md5}  {os.path.basename(dst)}")
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
                f"- 交付文件：`<样本>.PASS.adjudicated.vcf.gz`（+ `.tbi` 索引）× {len(rows)}\n"
                f"- 口径：GATK HaplotypeCaller（gVCF 联合分型）→ SNP/INDEL 硬过滤 → PASS →"
                f" 矩阵 `./.` 按 mosdepth 靶区深度 DP≥{config.DP_MIN} 裁决后重建的"
                f"每样本 VCF（仅替换 GT，其余字段原样保留）\n"
                f"- 校验：`md5sum -c md5sum.txt`；明细见 MANIFEST.tsv"
                f"（样本/文件/大小/md5/记录数/来源）\n"
                f"- 生成流程版本：{config.WORK_DIR}/pipeline（README 含完整口径与差异记录）\n")
    log.result(f"[交付] {dbatch}: VCF ×{len(rows)} + md5sum.txt + MANIFEST.tsv"
               + ("（README 更新）" if have else ""))
    return dbatch


def _step_anoms(ctx, step_prefix):
    """本步失败样本 → P0 异常行"""
    return [(f"P0", f"{sm} 执行失败：{reason.split(':', 1)[-1].strip()}")
            for sm, reason in ctx.failed.items() if reason.startswith(step_prefix)]


def _step_notify(notify_on, batch, step_no, name, log, **kw):
    """步骤里程碑通知。DINGTALK_MILESTONES=0 时只发异常级（P0/P1），OK 级静默。"""
    if not notify_on:
        return
    title, text = alerts.step_milestone(batch, step_no, name, **kw)
    level = alerts.worst_level(kw.get("anomalies"))
    if not config.DINGTALK_MILESTONES and level == "OK":
        return
    dingtalk.notify(title, text, logger=log)


def _notify_result(batch, bdata, plan, log):
    m = bdata.get("metrics", {})
    final_anoms = []
    if bdata.get("failed_samples"):
        final_anoms += [("P0", f"{len(bdata['failed_samples'])} 个样本失败: "
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
        f"workers {plan.workers}×{plan.threads}"
        f"\n\n质量: fastp 保留率 {_avg(m.get('fastp_retention'))}%｜"
        f"mapped {_avg(m.get('mapped_pct'))}%｜dup {_avg(m.get('dup_pct'))}%｜"
        f"20X {_avg(m.get('pct_20x'))}%"
        f"\n\n变异: raw→PASS SNP {m.get('snp_raw')}→{m.get('snp_pass')}｜"
        f"INDEL {m.get('indel_raw')}→{m.get('indel_pass')}｜"
        f"Ti/Tv {m.get('titv_raw')}→{m.get('titv_pass')}")
    if final_anoms:
        for lv, msg in final_anoms:
            text += f"\n\n异常: [{lv}] {msg}" + (" ← 需确认" if lv == "P1" else " ← 阻断/污染级")
    else:
        text += "\n\n异常: 无"
    delivery = (bdata.get("artifacts") or {}).get("delivery_dir")
    if delivery:
        text += f"\n\n交付: {delivery}（*.PASS.adjudicated.vcf.gz ×"
        text += f"{len([f for f in os.listdir(delivery) if f.endswith('.vcf.gz')])}"
        text += " + md5sum.txt + MANIFEST.tsv）"
    text += (f"\n\n产物: {bdata.get('work_dir')}（cohort.PASS / 矩阵 / 每样本裁决VCF / "
             f"MultiQC）"
             f"\n\n日志: tail -f {bdata.get('work_dir')}/logs/pipeline_*.log")
    dingtalk.notify(title, text, logger=log)


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
    args = parse_args()

    if args.notify_test:
        ok, err = dingtalk.send_markdown(
            "[GWAS] 测试通知",
            "#### [GWAS] 流程测试通知\n > 这是 GWAS pipeline 的钉钉通知测试\n"
            "| 项目 | 值 |\n| --- | --- |\n| 状态 | OK |")
        print(f"notify-test: {'发送成功' if ok else f'发送失败: {err}'}")
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
    input_root = args.input if os.path.isdir(args.input) \
        else os.path.join(config.WORK_DIR, args.input)
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

    summary = {
        "run": {"timestamp": ts, "argv": " ".join(sys.argv),
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

    # 交付总索引：delivery/INDEX.md（各批次交付目录一览）
    delivered = {b: d.get("artifacts", {}).get("delivery_dir")
                 for b, d in summary["batches"].items()}
    delivered = {b: p for b, p in delivered.items() if p}
    if delivered:
        idx = os.path.join(config.DELIVERY_DIR, "INDEX.md")
        os.makedirs(config.DELIVERY_DIR, exist_ok=True)
        with open(idx, "w", encoding="utf-8") as f:
            f.write(f"# 交付索引\n\n- 生成时间："
                    f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
                    f"| 批次 | 交付目录 | VCF 数 |\n| --- | --- | --- |\n")
            for b in sorted(delivered):
                n = len([x for x in os.listdir(delivered[b])
                         if x.endswith(".vcf.gz")])
                f.write(f"| {b} | `{delivered[b]}` | {n} |\n")
        summary["delivery_index"] = idx
        if last_summary_path and not args.dry_run:
            report_mod.write_run_summary(summary, last_summary_path)

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
