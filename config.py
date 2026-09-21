#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""配置中心：所有路径 / 线程 / 内存 / 阈值 / 钉钉 webhook。
环境参数三源优先级：进程环境变量 > pipeline/.env > 本文件内置默认
（钉钉地址等环境参数禁止硬编码进代码，唯一文件来源是 .env，DEC-18）。仅标准库。"""

import os
from pathlib import Path

# ── 流程版本（与 design_doc/DESIGN.md design-meta version 同步，
#    check_design.py 校验两处一致；写进 run_summary 与交付 README） ──────
PIPELINE_VERSION = "2.26.0"

# ── 工作区 ──────────────────────────────────────────────────────────────
PIPELINE_DIR = Path(__file__).resolve().parent          # $WORK/pipeline
WORK_DIR = str(PIPELINE_DIR.parent)                     # $WORK（--bind 挂载根）


# ── .env 环境参数文件（默认 pipeline/.env，与代码同目录；GWAS_ENV_FILE 可改址） ──
def parse_env_file(path):
    """解析 .env → {key: value}（纯解析不改环境；文件不存在返回 {}）。
    语法：KEY=VALUE 每行一条，值取首个 = 之后整段（可含 =?& 等字符）；
    支持 # 整行注释、空行、export 前缀、成对单/双引号（引号内原样保留）；
    未加引号的值支持行内 " #" 注释；键须为合法标识符，其余行忽略。"""
    parsed = {}
    if not os.path.isfile(path):
        return parsed
    with open(path, encoding="utf-8-sig") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export ") or line.startswith("export\t"):
                line = line[7:].lstrip()
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            key, val = key.strip(), val.strip()
            if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
                val = val[1:-1]                      # 引号值：原样（保留 # 与空格）
            elif " #" in val:
                val = val.split(" #", 1)[0].rstrip()  # 裸值：剥行内注释
            if key.isidentifier():
                parsed[key] = val
    return parsed


ENV_FILE = os.environ.get("GWAS_ENV_FILE", os.path.join(PIPELINE_DIR, ".env"))
ENV_FILE_KEYS = parse_env_file(ENV_FILE)   # .env 中已定义的键（值可能被进程环境覆盖）
for _k, _v in ENV_FILE_KEYS.items():       # setdefault：进程环境变量优先于 .env
    os.environ.setdefault(_k, _v)

RAW_DATA_DIR = os.environ.get("GWAS_RAW_DATA", os.path.join(WORK_DIR, "0_raw_data"))
RESULTS_ROOT = os.environ.get("GWAS_RESULTS", os.path.join(WORK_DIR, "results"))
# 独立交付目录（v2.3.0 由 delivery/ 改名 Output/）：最终交付文件
# *.PASS.adjudicated.vcf.gz(+.tbi) + MultiQC 报告 + md5sum.txt + MANIFEST
DELIVERY_DIR = os.environ.get("GWAS_DELIVERY_DIR", os.path.join(WORK_DIR, "Output"))
# 归档脚本（archive_results.py）输出目录：超期批次目录 7z 存放处（DEC-25）
ARCHIVE_DIR = os.environ.get("GWAS_ARCHIVE_DIR", os.path.join(WORK_DIR, "archive"))

# ── 容器运行时 ──────────────────────────────────────────────────────────
# 当前机器 singularity 实为 apptainer 1.4.5 别名；可切 apptainer
CONTAINER_RT = os.environ.get("CONTAINER_RT", "singularity")
# 镜像目录默认 pipeline 同级的 $WORK/singularity/（目录约定，不走 .env 默认值；
# GWAS_SIF_DIR 可重定向——测试/多工作区部署用）
SIF_DIR = os.environ.get("GWAS_SIF_DIR", os.path.join(WORK_DIR, "singularity"))
SIF = {
    "fastqc":    os.path.join(SIF_DIR, "fastqc_0.12.1.sif"),
    "fastp":     os.path.join(SIF_DIR, "fastp_1.3.6.sif"),
    "bwa":       os.path.join(SIF_DIR, "bwa-mem2_2.3.sif"),
    "samtools":  os.path.join(SIF_DIR, "samtools_1.24.sif"),
    "gatk":      os.path.join(SIF_DIR, "gatk_4.6.2.0.sif"),
    "bcftools":  os.path.join(SIF_DIR, "bcftools_1.24.sif"),
    "mosdepth":  os.path.join(SIF_DIR, "mosdepth_0.3.14.sif"),
    "multiqc":   os.path.join(SIF_DIR, "multiqc_1.35.sif"),
}

# ── 参考文件（GWAS_REFERENCE_DIR 可重定向，默认 $WORK/reference/） ──────
REF_DIR = os.environ.get("GWAS_REFERENCE_DIR", os.path.join(WORK_DIR, "reference"))
GENOME_FA = os.path.join(REF_DIR, "genome", "genome.fa")          # hg38 + bwa-mem2 索引
GENOME_DICT = os.path.join(REF_DIR, "genome", "genome.dict")
DBSNP_VCF = os.path.join(REF_DIR, "Homo_sapiens_assembly38.dbsnp138.vcf.gz")
MILLS_VCF = os.path.join(REF_DIR, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
GNOMAD_VCF = os.path.join(REF_DIR, "af-only-gnomad.hg38.vcf.gz")  # 备用
# 靶区 bed 属私密文件（绝不入仓库，.gitignore 拦截 *.bed/*.interval_list）：
# GWAS_TARGETS_BED 可经 .env 改址（默认 $WORK/reference/targets.bed），
# 运行时派生的 sorted.bed / interval_list 与 bed 同目录自动生成
TARGETS_BED = os.environ.get("GWAS_TARGETS_BED", os.path.join(REF_DIR, "targets.bed"))
TARGETS_DIR = os.path.dirname(TARGETS_BED)
TARGETS_SORTED_BED = os.path.join(TARGETS_DIR, "targets.sorted.bed")
TARGETS_INTERVAL_LIST = os.path.join(TARGETS_DIR, "targets.sorted.interval_list")

# ── 工具参数（除命令行覆盖外，均由 resource.py 规划，禁止写死） ─────────
FASTP_LENGTH_REQUIRED = 36
HC_INTERVAL_PADDING = 100
PLOIDY = 2
BWA_K = 100000000
# 统计类命令统一超时（flagstat/stats/count_records/view|wc/query 等秒级查询）：
# 挂死不再无限阻塞批次；GATK/bwa 等长任务不设超时（None）
STATS_TIMEOUT_S = 600

# 硬过滤阈值（GATK best practices，样本量约 30，不用 VQSR）——与笔记 4-13/4-14 逐字一致
SNP_HARD_FILTERS = [
    ("QD < 2.0", "QD2"),
    ("QUAL < 30.0", "QUAL30"),
    ("SOR > 3.0", "SOR3"),
    ("FS > 60.0", "FS60"),
    ("MQ < 40.0", "MQ40"),
    ("MQRankSum < -12.5", "MQRankSum"),
    ("ReadPosRankSum < -8.0", "ReadPosRankSum"),
]
INDEL_HARD_FILTERS = [
    ("QD < 2.0", "QD2"),
    ("QUAL < 30.0", "QUAL30"),
    ("SOR > 10.0", "SOR10"),
    ("FS > 200.0", "FS200"),
    ("ReadPosRankSum < -20.0", "ReadPosRankSum"),
]

# ── QC 阈值 ─────────────────────────────────────────────────────────────
FASTP_RETENTION_WARN = 90.0     # 保留率 %，低于提示（与 P1 线 80 之间构成提示带）。
                                # v2.21.0/DEC-31 由 95 下调：260918 批次实测全部样本
                                # 92.18-94.47%，为本 panel + --length_required 36
                                # 修剪口径下的批次正常基线——95 线致提示带全批点名、
                                # 失去区分度；90 后只点名真偏低样本（RUN-45）
MAPPED_MIN_PCT = 90.0           # 比对率 QC 口径下限（<90 → P1 报错不中断，v2.9.0）
DUP_WARN_PCT = 30.0             # 重复率上限（>30 且 ELS 偏小 → 文库复杂度不足）
MEAN_COV_MIN = 50.0             # MEAN/MED_TARGET_COVERAGE 下限
PCT_20X_MIN = 90.0              # PCT_TARGET_BASES_20X 下限 %
DP_MIN = int(os.environ.get("GWAS_DP_MIN", "20"))   # ./. 裁决深度阈值（可配）

# ── 分级告警阈值（v2.9.0 三级体系 DEC-21：P0=阻断级·严重影响分析→中断批次；
#    P1=严重·执行失败隔离/NTC 污染/mapped<90，报错不中断；P2=质量提示·只记录；
#    环境变量可覆盖） ─────────────────────────────────────────────────────
DISK_MIN_FREE_GB = float(os.environ.get("GWAS_DISK_MIN_FREE_GB", "200"))  # 开跑前剩余磁盘
DISK_PER_SAMPLE_GB = float(os.environ.get("GWAS_DISK_PER_SAMPLE_GB", "75"))  # 每样本估算需求
FASTP_RETENTION_P1 = 80.0       # Step1: 保留率 <80% → P2
FASTP_Q30_P1 = 85.0             # Step1: Q30 <85% → P2
MAPPED_NOTIFY_P1 = 95.0         # Step2: mapped <95% → P2（严于 QC 口径 90→P1）
PROPER_PAIR_P1 = 85.0           # Step2: properly paired <85% → P2
DUP_P1 = 30.0                   # Step3: 重复率 >30% → P2
PCT_SELECTED_P1 = 85.0          # Step6: 捕获效率 PCT_SELECTED_BASES <85% → P2（v2.19.0/DEC-29：
                                # 新 panel 本地实测 89.4-90.6%，外送同批 pct_selected_bases
                                # 89.6-90.7% 交叉验证一致；老 panel 26-40% 不适用本阈值）
MEAN_DEPTH_P1 = 50.0            # Step6: mean depth <50× → P2
PCT_20X_P1 = 95.0               # Step6: ≥20x 靶比例 <95% → P2
TITV_P1 = 2.0                   # Step6: Ti/Tv <2.0 → P2（结论口径）
CALL_RATE_P1 = 95.0             # Step6: call rate <95% → P2（非 ./. 基因型比例）
NTC_DEPTH_P0 = 10.0             # Step6: NTC 靶区深度 >10× → P1（污染，报错不中断）
READS_MIN = int(os.environ.get("GWAS_READS_MIN", "1000000"))
                                # Step1: 样本 reads 绝对量 <1M → P1（上样不足）
NTC_READS_PCT_P2 = float(os.environ.get("GWAS_NTC_READS_PCT_P2", "1.0"))
                                # Step1: NTC reads 占批次中位样本 % >1% → P2（污染维度之二）
DEPTH_CV_P2 = float(os.environ.get("GWAS_DEPTH_CV_P2", "0.5"))
                                # Step6: 批次内 mean depth 变异系数 CV >0.5 → P2（混入异常样本）
RECAL_OBS_MIN_P2 = float(os.environ.get("GWAS_RECAL_OBS_MIN_P2", "100000"))
                                # Step4: recal M 事件观测数 <1e5 → P2（known-sites 覆盖崩坏，校准不可信）
MD5_MANIFEST_MAX_BYTES = int(os.environ.get("GWAS_MD5_MANIFEST_MAX_BYTES", str(500 * 1024)))
                                # Step0: md5 清单识别体积上限（TH-37，DEC-36——防同名大文本文件误认）

# ── 钉钉机器人（凭证等环境参数只来自 pipeline/.env 或进程环境变量，禁止硬编码） ──
# v2.15.0（DEC-24）起为企业内部应用机器人（可发文件）；WEBHOOK/KEYWORD 为
# webhook 时代遗留键，保留解析但不再参与发送
DINGTALK_WEBHOOK = os.environ.get("DINGTALK_WEBHOOK", "")   # 遗留（webhook 机器人，已停用）
DINGTALK_KEYWORD = os.environ.get("DINGTALK_KEYWORD", "gwas")   # 遗留（关键词校验，已停用）
DINGTALK_MILESTONES = os.environ.get("DINGTALK_MILESTONES", "1") == "1"  # 里程碑通知开关
DINGTALK_CLIENT_ID = os.environ.get("DINGTALK_CLIENT_ID", "")          # 企业应用 appKey
DINGTALK_CLIENT_SECRET = os.environ.get("DINGTALK_CLIENT_SECRET", "")  # 企业应用密钥（.env，600）
DINGTALK_ROBOT_CODE = os.environ.get("DINGTALK_ROBOT_CODE", "")        # 机器人编码（空=复用 CLIENT_ID）
DINGTALK_CONVERSATION_ID = os.environ.get("DINGTALK_CONVERSATION_ID", "")   # 群 openConversationId
# 未配置企业凭证 → 通知静默跳过（启动 WARN + 配置指引）

# ── 资源规划常量 ────────────────────────────────────────────────────────
RESERVE_CORES = 2                          # 系统预留核数
BWA_INDEX_MEM_GB = 17.0                    # bwa-mem2 索引 mmap 常驻（实测 10+6.2+0.8）
FASTP_BUF_GB = 1.0                         # fastp/杂项进程缓冲
SORT_MEM_MIN, SORT_MEM_MAX = 128, 2048     # MB（钳位）
GATK_MEM_MIN_GB, GATK_MEM_MAX_GB = 1, 8    # -Xmx 钳位
COHORT_MEM_MIN_GB, COHORT_MEM_MAX_GB = 8, 32   # 串行 cohort 步骤（可用内存 60%）
PROFILE_LOW = {"cpu": 16, "mem_gb": 20}    # 低配档基准
PROFILE_HIGH = {"cpu": 100, "mem_gb": 900} # 高配档基准（更大机器 auto 照常用）

# 通知 / 报告
MULTIQC_TITLE = "GWAS Panel 下机数据 QC"


def batch_result_dir(batch, run_date=""):
    """批次结果目录：results/<批次>_<执行日期>/（多次运行相互隔离；
    同日重跑同目录 → 幂等断点续跑；run_date 为 YYYYMMDD）"""
    return os.path.join(RESULTS_ROOT, f"{batch}_{run_date}" if run_date else batch)


def batch_delivery_dir(batch, run_date=""):
    """批次交付目录：Output/<批次>_<执行日期>/（与结果目录同构）"""
    return os.path.join(DELIVERY_DIR, f"{batch}_{run_date}" if run_date else batch)
