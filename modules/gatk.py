#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gatk：MarkDuplicates / BQSR / HaplotypeCaller / CombineGVCFs / GenotypeGVCFs /
SelectVariants / VariantFiltration / CollectHsMetrics 分函数 + interval_list 准备。
命令与笔记《3-去重+校准》《4-变异检测》《5-测序质量》逐字对应；
笔记中写死的 -Xmx2g/-t 8 等资源参数由 resource.py 规划值代入（有意工程化，见 README）。

全部命令固定 -Xrs（DEC-33/RUN-48；定位修订 DEC-43/RUN-59）：JVM 启动时对
SIGHUP 安装自己的处理器，覆盖继承的忽略位——-Xrs 使其不装。v2.29.0 实测证明
这只是纵深防御的第二层：apptainer 容器内部进程链会重置信号继承位，SIG_IGN 与
-Xrs 都护不住容器内 JVM（12 个 GATK 同瞬 "Hangup" exit=129）——根治层是
run_pipeline.ensure_detached() 的 setsid 脱离控制终端（SIGHUP 无从发出）。
代价：SIGQUIT/SIGTERM 优雅停机与 kill -3 线程转储不可用（排障改用
jcmd Thread.print），流程超时兜底本就 SIGKILL，无影响。"""

import os

import config
from runner import nonempty


# ── interval_list 准备（DEC-26 参数 + DEC-28 字典过滤/新鲜度） ────────────
def _dict_contigs(dict_path):
    """genome.dict 的 SN 名集合（@SQ 行第二字段去 SN: 前缀）"""
    try:
        with open(dict_path, encoding="utf-8") as f:
            return {ln.split("\t")[1][3:] for ln in f
                    if ln.startswith("@SQ\t")}
    except OSError:
        return set()


def _bed_contigs(bed_path):
    try:
        with open(bed_path, encoding="utf-8") as f:
            return {ln.split("\t", 1)[0].strip() for ln in f
                    if ln.strip() and not ln.startswith(("#", "track", "browser"))}
    except OSError:
        return set()


def _derived_stale(out_path, src_path):
    """派生靶区文件新鲜度（DEC-28）：缺失/空，或源文件 mtime 更新 → 需重建。
    防"更新 targets.bed 但派生 sorted.bed/interval_list 非空被幂等跳过"的陈旧陷阱
    （RUN-42：sorted.bed 符号链接含 ALT contig 行，HC -L 全军覆没）"""
    if not nonempty(out_path):
        return True
    try:
        return os.path.getmtime(out_path) < os.path.getmtime(src_path)
    except OSError:
        return True


def prep_interval_list(runner, logger):
    if _derived_stale(config.TARGETS_SORTED_BED, config.TARGETS_BED):
        dropped = _bed_contigs(config.TARGETS_BED) - _dict_contigs(config.GENOME_DICT)
        if dropped:
            logger.warn("靶区 bed 含 genome.dict 外 contig，生成 sorted.bed 时过滤丢弃: "
                        + ", ".join(sorted(dropped)))
        logger.info(f"生成 {config.TARGETS_SORTED_BED}"
                    "（awk 前三列 + genome.dict contig 过滤 + 版本排序去重）")
        # outputs 不传 runner：陈旧但非空的派生文件会被 runner 自身幂等误跳过
        # （DEC-28——新鲜度判断上收 prep；字典过滤在 awk 内做，dry-run 零落盘不变）
        rc = runner.run(
            f"awk 'NR==FNR && /^@SQ/ {{sub(/^SN:/, \"\", $2); c[$2]; next}} "
            f"NR>FNR && NF>=3 && ($1 in c) {{print $1\"\\t\"$2\"\\t\"$3}}' "
            f"{config.GENOME_DICT} {config.TARGETS_BED} "
            f"| sort -k1,1V -k2,2n -u > {config.TARGETS_SORTED_BED}",
            logger=logger)
        if rc != 0:
            return False
        if not getattr(runner, "dry_run", False) and not nonempty(config.TARGETS_SORTED_BED):
            logger.error(f"sorted.bed 生成后为空: {config.TARGETS_SORTED_BED}")
            return False
    if _derived_stale(config.TARGETS_INTERVAL_LIST, config.TARGETS_SORTED_BED):
        logger.info(f"生成 {config.TARGETS_INTERVAL_LIST}（BedToIntervalList -SD genome.dict，"
                    "--UNIQUE 去重合并 + --DROP_MISSING_CONTIGS 丢字典外 contig，DEC-26）")
        rc = runner.run(
            runner.tool("gatk",
                        f'gatk --java-options "-Xrs -Xmx2g" BedToIntervalList '
                        f"-I {runner.cpath(config.TARGETS_SORTED_BED)} "
                        f"-O {runner.cpath(config.TARGETS_INTERVAL_LIST)} "
                        f"-SD {runner.cpath(config.GENOME_DICT)} "
                        f"--UNIQUE true --DROP_MISSING_CONTIGS true"),
            logger=logger)
        if rc != 0:
            return False
    return True


# ── Step 3 MarkDuplicates ──
def markdup(runner, in_bam, out_bam, metrics, gatk_mem, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" MarkDuplicates '
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
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" BaseRecalibrator '
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
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" ApplyBQSR '
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
    """-L 用 interval_list 不用 sorted.bed（DEC-28）：GATK 引擎对 -L 区间 contig
    严格校验字典，bed 含字典外 contig（如 ALT）即 USER ERROR 全体失败（RUN-42）；
    interval_list 经 DEC-26 参数保证字典口径+唯一合并。bed 留给 bcftools/mosdepth
    （二者容忍字典外 contig），且 sorted.bed 生成同样已按字典过滤。"""
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" HaplotypeCaller '
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-I {runner.cpath(bqsr_bam)} "
                    f"-O {runner.cpath(out_gvcf)} "
                    f"-ERC GVCF "
                    f"-L {runner.cpath(config.TARGETS_INTERVAL_LIST)} "
                    f"--interval-padding {config.HC_INTERVAL_PADDING} "
                    f"--native-pair-hmm-threads {hmm_threads} "
                    f"-ploidy {config.PLOIDY}"),
        logger=logger, outputs=[out_gvcf, str(out_gvcf) + ".tbi"])
    return rc == 0


def combine_gvcfs(runner, gvcf_list_host, out_host, cohort_mem, logger):
    """gvcf.list 每次重新生成；小 panel 用 CombineGVCFs，不用 GenomicsDB"""
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{cohort_mem}" CombineGVCFs '
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-V {runner.cpath(gvcf_list_host)} "
                    f"-O {runner.cpath(out_host)}"),
        logger=logger, outputs=[out_host, str(out_host) + ".tbi"])
    return rc == 0


def genotype_gvcfs(runner, combined_gvcf, out_host, cohort_mem, logger):
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{cohort_mem}" GenotypeGVCFs '
                    f"-R {runner.cpath(config.GENOME_FA)} "
                    f"-V {runner.cpath(combined_gvcf)} "
                    f"-O {runner.cpath(out_host)}"),
        logger=logger, outputs=[out_host, str(out_host) + ".tbi"])
    return rc == 0


def select_variants(runner, in_vcf, vtype, out_vcf, gatk_mem, logger):
    """SNP / INDEL 分拣（必须先 norm 摊平，防 MIXED 整条丢弃）"""
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" SelectVariants '
                    f"-V {runner.cpath(in_vcf)} --select-type {vtype} "
                    f"-O {runner.cpath(out_vcf)}"),
        logger=logger, outputs=[out_vcf])
    return rc == 0


def variant_filtration(runner, in_vcf, out_vcf, filters, gatk_mem, logger):
    """硬过滤（样本量约 30，不用 VQSR）。filters: [(expr, name), ...]"""
    fargs = " ".join(f'-filter "{expr}" --filter-name {name}' for expr, name in filters)
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" VariantFiltration '
                    f"-V {runner.cpath(in_vcf)} {fargs} "
                    f"-O {runner.cpath(out_vcf)}"),
        logger=logger, outputs=[out_vcf])
    return rc == 0


# ── Step 6 HsMetrics ──
def collect_hsmetrics(runner, bam, out_txt, gatk_mem, logger):
    """BAIT_INTERVALS = TARGET_INTERVALS = targets.sorted.interval_list（bait 与 target 同一文件）"""
    rc = runner.run(
        runner.tool("gatk",
                    f'gatk --java-options "-Xrs -Xmx{gatk_mem}" CollectHsMetrics '
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
    """MEAN/MED ≥50×、20X ≥90%、捕获效率 PCT_SELECTED ≥85%（数据质量告警级）。
    捕获效率口径 v2.19.0/DEC-29：PCT_SELECTED_BASES（on+near bait / 比对碱基，
    与外送报告 pct_selected_bases 同口径）；on-target 仅信息指标——1bp SNP panel
    下 ON_TARGET_BASES 只计区间内碱基（≈0.6%），是 panel 几何产物不反映捕获好坏，
    旧文档"PCT_SELECTED(25-28%) 不误读为捕获效率"是老 panel 时代的结论，已废止。
    control=True（NTC 等对照）：近零覆盖是预期，只记录不判阈值。"""
    mean, med = hs.get("MEAN_TARGET_COVERAGE"), hs.get("MEDIAN_TARGET_COVERAGE")
    p20 = hs.get("PCT_TARGET_BASES_20X")
    psel = hs.get("PCT_SELECTED_BASES")
    level, tag = (logger.info, "对照样本") if control else (logger.warn, "告警")
    if control:
        logger.result(f"HsMetrics[对照]: MEAN={mean}× 20X="
                      f"{p20*100 if p20 is not None else '?'}%（阴性对照近零覆盖属预期）")
        return True
    if mean is not None and mean < config.MEAN_COV_MIN:
        level(f"MEAN_TARGET_COVERAGE {mean}× < {config.MEAN_COV_MIN}×（{tag}）")
    if p20 is not None and p20 * 100 < config.PCT_20X_MIN:
        level(f"PCT_TARGET_BASES_20X {p20*100:.1f}% < {config.PCT_20X_MIN}%（{tag}）")
    if psel is not None and psel * 100 < config.PCT_SELECTED_P1:
        level(f"捕获效率 PCT_SELECTED {psel*100:.1f}% < {config.PCT_SELECTED_P1}%（{tag}）")
    logger.result(
        f"HsMetrics: MEAN={mean}× MED={med}× 20X={p20*100 if p20 is not None else '?'}% "
        f"捕获效率PCT_SELECTED={psel*100 if psel is not None else '?'}% "
        f"on-bait={hs.get('ON_BAIT_PCT')}% on-target={hs.get('ON_TARGET_PCT')}%（信息指标，不告警）")
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
