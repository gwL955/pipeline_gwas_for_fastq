#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gatk：MarkDuplicates / BQSR / HaplotypeCaller / CombineGVCFs / GenotypeGVCFs /
SelectVariants / VariantFiltration / CollectHsMetrics 分函数 + interval_list 准备。
命令与笔记《3-去重+校准》《4-变异检测》《5-测序质量》逐字对应；
笔记中写死的 -Xmx2g/-t 8 等资源参数由 resource.py 规划值代入（有意工程化，见 README）。"""

import os

import config
from runner import nonempty


# ── interval_list 准备（targets.sorted.bed / interval_list 缺失时自动生成，幂等跳过） ──
def prep_interval_list(runner, logger):
    if not nonempty(config.TARGETS_SORTED_BED):
        logger.info(f"生成 {config.TARGETS_SORTED_BED}（awk 前三列 + 版本排序去重）")
        rc = runner.run(
            f"awk 'NF>=3{{print $1\"\\t\"$2\"\\t\"$3}}' {config.TARGETS_BED} "
            f"| sort -k1,1V -k2,2n -u > {config.TARGETS_SORTED_BED}",
            logger=logger, outputs=[config.TARGETS_SORTED_BED])
        if rc != 0:
            return False
    if not nonempty(config.TARGETS_INTERVAL_LIST):
        logger.info(f"生成 {config.TARGETS_INTERVAL_LIST}（BedToIntervalList -SD genome.dict）")
        rc = runner.run(
            runner.tool("gatk",
                        f"gatk --java-options -Xmx2g BedToIntervalList "
                        f"-I {runner.cpath(config.TARGETS_SORTED_BED)} "
                        f"-O {runner.cpath(config.TARGETS_INTERVAL_LIST)} "
                        f"-SD {runner.cpath(config.GENOME_DICT)}"),
            logger=logger, outputs=[config.TARGETS_INTERVAL_LIST])
        if rc != 0:
            return False
    return True


# ── Step 3 MarkDuplicates ──
def markdup(runner, in_bam, out_bam, metrics, gatk_mem, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} MarkDuplicates "
                    f"-I {runner.cpath(in_bam)} -O {runner.cpath(out_bam)} "
                    f"-M {runner.cpath(metrics)} --CREATE_INDEX true"),
        logger=logger, outputs=[out_bam, metrics])
    return rc == 0


def parse_markdup_metrics(path):
    """按表头名解析：READ_PAIRS_EXAMINED / PERCENT_DUPLICATION / ESTIMATED_LIBRARY_SIZE"""
    d = {}
    try:
        with open(path, encoding="utf-8") as f:
            header = None
            for line in f:
                if line.startswith("LIBRARY"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header and line.strip() and not line.startswith("#"):
                    fields = line.rstrip("\n").split("\t")
                    row = dict(zip(header, fields))
                    for k, cast in (("READ_PAIRS_EXAMINED", int),
                                    ("READ_PAIR_DUPLICATES", int),
                                    ("PERCENT_DUPLICATION", float),
                                    ("READ_PAIR_OPTICAL_DUPLICATES", int)):
                        try:
                            d[k] = cast(row.get(k, "") or 0)
                        except ValueError:
                            d[k] = row.get(k)
                    els = row.get("ESTIMATED_LIBRARY_SIZE", "")
                    d["ESTIMATED_LIBRARY_SIZE"] = int(els) if els not in ("", "?") else None
                    break
    except OSError:
        pass
    return d


def qc_duplication(metrics, logger):
    """dup>30% 且 ELS 偏小 → 文库复杂度不足告警"""
    dup = metrics.get("PERCENT_DUPLICATION")
    pairs = metrics.get("READ_PAIRS_EXAMINED", 0) or 0
    els = metrics.get("ESTIMATED_LIBRARY_SIZE")
    if dup is not None and dup > 0.30 and els is not None and pairs \
            and els < pairs:
        logger.warn(f"dup {dup*100:.1f}% >30% 且 ELS({els}) < 测序对数({pairs})"
                    f" → 文库复杂度不足，建议重建文库")
    elif dup is not None and dup > 0.30:
        logger.warn(f"dup {dup*100:.1f}% >30%（ELS={els}，若 ELS 足够大仅为测太深）")


# ── Step 4 BQSR ──
def base_recalibrator(runner, markdup_bam, table, gatk_mem, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} BaseRecalibrator "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-I {runner.cpath(markdup_bam)} "
                    f"--known-sites {runner.cpath(config.DBSNP_VCF)} "
                    f"--known-sites {runner.cpath(config.MILLS_VCF)} "
                    f"-O {runner.cpath(table)}"),
        logger=logger, outputs=[table])
    return rc == 0


def apply_bqsr(runner, markdup_bam, table, out_bam, gatk_mem, logger):
    """ApplyBQSR 输入是 markdup.bam（不是任何中间产物）"""
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} ApplyBQSR "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-I {runner.cpath(markdup_bam)} "
                    f"--bqsr-recal-file {runner.cpath(table)} "
                    f"-O {runner.cpath(out_bam)}"),
        logger=logger, outputs=[out_bam])
    if rc == 0:   # HC 需要索引，MarkDuplicates 已建 .bai，BQSR 产物补 index
        from modules import samtools as msam
        if not nonempty(out_bam + ".bai"):
            msam.run_index(runner, out_bam, logger)
    return rc == 0


# ── Step 5 变异检测 ──
def haplotypecaller(runner, bqsr_bam, out_gvcf, gatk_mem, hmm_threads, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} HaplotypeCaller "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-I {runner.cpath(bqsr_bam)} "
                    f"-O {runner.cpath(out_gvcf)} "
                    f"-ERC GVCF "
                    f"-L {runner.cpath(config.TARGETS_SORTED_BED)} "
                    f"--interval-padding {config.HC_INTERVAL_PADDING} "
                    f"--native-pair-hmm-threads {hmm_threads} "
                    f"-ploidy {config.PLOIDY}"),
        logger=logger, outputs=[out_gvcf, str(out_gvcf) + ".tbi"])
    return rc == 0


def combine_gvcfs(runner, gvcf_list_host, out_host, cohort_mem, logger):
    """gvcf.list 每次重新生成；小 panel 用 CombineGVCFs，不用 GenomicsDB"""
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{cohort_mem} CombineGVCFs "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-V {runner.cpath(gvcf_list_host)} "
                    f"-O {runner.cpath(out_host)}"),
        logger=logger, outputs=[out_host, str(out_host) + ".tbi"])
    return rc == 0


def genotype_gvcfs(runner, combined_gvcf, out_host, cohort_mem, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{cohort_mem} GenotypeGVCFs "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-V {runner.cpath(combined_gvcf)} "
                    f"-O {runner.cpath(out_host)}"),
        logger=logger, outputs=[out_host, str(out_host) + ".tbi"])
    return rc == 0


def select_variants(runner, in_vcf, vtype, out_vcf, gatk_mem, logger):
    """SNP / INDEL 分拣（必须先 norm 摊平，防 MIXED 整条丢弃）"""
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} SelectVariants "
                    f"-V {runner.cpath(in_vcf)} --select-type {vtype} "
                    f"-O {runner.cpath(out_vcf)}"),
        logger=logger, outputs=[out_vcf])
    return rc == 0


def variant_filtration(runner, in_vcf, out_vcf, filters, gatk_mem, logger):
    """硬过滤（样本量约 30，不用 VQSR）。filters: [(expr, name), ...]"""
    fargs = " ".join(f'-filter "{expr}" --filter-name {name}' for expr, name in filters)
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} VariantFiltration "
                    f"-V {runner.cpath(in_vcf)} {fargs} "
                    f"-O {runner.cpath(out_vcf)}"),
        logger=logger, outputs=[out_vcf])
    return rc == 0


# ── Step 6 HsMetrics ──
def collect_hsmetrics(runner, bam, out_txt, gatk_mem, logger):
    """BAIT_INTERVALS = TARGET_INTERVALS = targets.sorted.interval_list（bait 与 target 同一文件）"""
    rc = runner.run(
        runner.tool("gatk",
                    f"gatk --java-options -Xmx{gatk_mem} CollectHsMetrics "
                    f"-I {runner.cpath(bam)} "
                    f"-O {runner.cpath(out_txt)} "
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-BAIT_INTERVALS {runner.cpath(config.TARGETS_INTERVAL_LIST)} "
                    f"-TARGET_INTERVALS {runner.cpath(config.TARGETS_INTERVAL_LIST)}"),
        logger=logger, outputs=[out_txt])
    return rc == 0


def parse_hsmetrics(path):
    """按表头名解析关键列；on-bait/on-target 用 PF_UQ_BASES_ALIGNED（去重后口径）"""
    d = {}
    try:
        with open(path, encoding="utf-8") as f:
            header = None
            for line in f:
                if line.startswith("BAIT_SET"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header and line.startswith("targets\t"):
                    row = dict(zip(header, line.rstrip("\n").split("\t")))

                    def num(key, default=0.0):
                        try:
                            return float(row.get(key, "") or default)
                        except ValueError:
                            return default
                    d["MEAN_TARGET_COVERAGE"] = num("MEAN_TARGET_COVERAGE")
                    d["MEDIAN_TARGET_COVERAGE"] = num("MEDIAN_TARGET_COVERAGE")
                    d["PCT_TARGET_BASES_20X"] = num("PCT_TARGET_BASES_20X")
                    d["PCT_SELECTED_BASES"] = num("PCT_SELECTED_BASES")
                    d["FOLD_ENRICHMENT"] = num("FOLD_ENRICHMENT")
                    d["FOLD_80_BASE_PENALTY"] = num("FOLD_80_BASE_PENALTY")
                    pq = num("PF_UQ_BASES_ALIGNED", 1.0) or 1.0
                    d["ON_BAIT_PCT"] = round(num("ON_BAIT_BASES") / pq * 100, 2)
                    d["ON_TARGET_PCT"] = round(num("ON_TARGET_BASES") / pq * 100, 2)
                    break
    except OSError:
        pass
    return d


def qc_hsmetrics(hs, logger, control=False):
    """MEAN/MED ≥50×、20X ≥90%（数据质量告警级）；on-target≈on-bait（约 8-10%）；
    PCT_SELECTED_BASES（25-28%）是含邻域口径，不得误读为捕获效率。
    control=True（NTC 等对照）：近零覆盖是预期，只记录不判阈值。"""
    mean, med = hs.get("MEAN_TARGET_COVERAGE"), hs.get("MEDIAN_TARGET_COVERAGE")
    p20 = hs.get("PCT_TARGET_BASES_20X")
    level, tag = (logger.info, "对照样本") if control else (logger.warn, "告警")
    if control:
        logger.result(f"HsMetrics[对照]: MEAN={mean}× 20X="
                      f"{p20*100 if p20 is not None else '?'}%（阴性对照近零覆盖属预期）")
        return True
    if mean is not None and mean < config.MEAN_COV_MIN:
        level(f"MEAN_TARGET_COVERAGE {mean}× < {config.MEAN_COV_MIN}×（{tag}）")
    if p20 is not None and p20 * 100 < config.PCT_20X_MIN:
        level(f"PCT_TARGET_BASES_20X {p20*100:.1f}% < {config.PCT_20X_MIN}%（{tag}）")
    logger.result(
        f"HsMetrics: MEAN={mean}× MED={med}× 20X={p20*100 if p20 is not None else '?'}% "
        f"on-bait={hs.get('ON_BAIT_PCT')}% on-target={hs.get('ON_TARGET_PCT')}% "
        f"PCT_SELECTED={hs.get('PCT_SELECTED_BASES')}（含±250bp 邻域口径，非捕获效率）")
    return True


def parse_recal_observations(table_path):
    """RecalTable1 中 EventType=M 行的 Observations 求和——BQSR 校准可信度信号
    （known-sites 覆盖崩坏时观测数骤降，RUN-34）。文件缺失/无有效行返回 None"""
    try:
        with open(table_path, encoding="utf-8") as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    hdr, start = None, None
    for i, ln in enumerate(lines):
        if ln.startswith("#:GATKTable:RecalTable1:"):
            hdr, start = lines[i + 1].split(), i + 2
            break
    if not hdr or "Observations" not in hdr or "EventType" not in hdr:
        return None
    oi, ei = hdr.index("Observations"), hdr.index("EventType")
    total = 0.0
    for ln in lines[start:]:
        if ln.startswith("#"):
            break                      # 下一表开始
        p = ln.split()
        if len(p) > max(oi, ei) and p[ei] == "M":
            try:
                total += float(p[oi])
            except ValueError:
                pass
    return total or None
