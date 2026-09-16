# GWAS 靶向测序 FASTQ→VCF 流程 · 设计文档（人机共写）

```yaml
# ---- design-meta（机器可解析锚点，勿手改格式；版本规则见 §0）----
doc: GWAS-pipeline-design
version: 2.11.0
updated: 2026-09-16
owner_human: gewenlong
owner_machine: ZCode(GLM)
source_of_truth:
  goals_and_scope: pipeline/design_doc/DESIGN.md          # 本文档（人写区为唯一目标来源；自含重建规格）
  thresholds: pipeline/config.py                 # 阈值唯一实现源；本文档 TH 表为镜像
  commands_reference: pipeline/design_doc/notes_code_reference.md   # 思源笔记《下机数据流程》快照
  run_evidence: pipeline/design_doc/RUN_HISTORY.md        # 机器回写运行台账
  usage: pipeline/README.md                       # 使用说明/快速上手（本文档为规格源）
consistency_check: pipeline/check_design.py      # TH 表 ↔ config.py 防漂移 + 版本双源校验
```

> **本文档目标（v2.11.0 起）**：自含——依据本文档即可从零重建整个项目（代码结构、
> 每步行为、全部阈值与分级、配置接口、验收标准）；历史沿革归 §10 CHANGELOG 与
> RUN_HISTORY.md，正文只陈述当前有效规格。

---

## 0. 文档使用协议（人机共写规范）

**目的**：以本文档为中心演进项目，防止反复修改导致目标漂移。

**分区标记**：
- `<!-- HUMAN -->` 区：目标、范围、口径、验收标准——**人**编写/修改，机器只读；
- `<!-- MACHINE -->` 区：状态、证据、运行结果——**机器**回写，人可批注（批注行以 `>>` 开头）。

**编号锚点**（代码注释、测试、日志、通知引用这些编号，形成可追溯链）：
- `REQ-xx` 需求；`DEC-xx` 设计决策；`TH-xx` 阈值；`RUN-xx` 运行轮次（台账在 RUN_HISTORY.md）。

**修改流程**：
1. 人修改 HUMAN 区（目标/阈值/验收）→ 在 §10 CHANGELOG 记一条（版本号+1：口径微调 `minor+1`，目标/范围变更 `major+1`）；
2. 机器解析差异 → 对照实现（代码/测试同步引用编号）；
3. 机器回写 MACHINE 区状态与 RUN_HISTORY 台账；
4. 运行 `python3 pipeline/check_design.py`——TH 表与 config.py 不一致、或 design-meta
   version 与 `config.PIPELINE_VERSION` 漂移即失败；`./pipeline/run_tests.sh` 末尾自动执行。

**冲突裁决顺序**：分析学内容以思源笔记原文为准（可经思源 MCP 回查）→ 工程化要求以本文档为准 → 笔记写死而本文档改为规划值的参数属有意工程化（记录于 pipeline/README.md §8"与笔记差异"）。

**废弃处理**：已废弃的目标/决策/阈值从正文删除，仅保留 CHANGELOG 中的变更记录与编号不复用声明（如 TH-17~20）。

---

## 1. 项目目标与范围 <!-- HUMAN -->

把「人类 GWAS 靶向测序：FASTQ → 质控修剪 → 比对 → 去重 → BQSR 校准 → 变异检测
（gVCF 联合分型 → 硬过滤 → VCF）→ 测序质量汇总」
实现为一套 **Python3 仅标准库、模块化、可断点续跑、资源自适应（16线程/20G ～ 100线程/900G+）、
带钉钉三级分级通知** 的流程控制程序；最终交付物为每样本 `*.PASS.adjudicated.vcf.gz`
（独立 Output 交付目录 + MultiQC 汇总报告 + md5 校验）。

**范围边界**：批次是最小分析单元且完全独立（严禁跨批次合并）；`0_raw_data/`、`back/`
只读；分析产物写 `results/`，交付产物写独立 `Output/`；NTC 等对照纳入 QC、排除出
联合变异检测；不使用 VQSR（样本量约 30 用硬过滤）。

## 2. 需求清单与状态 <!-- REQ 目标列=HUMAN；状态列=MACHINE -->

| 编号 | 需求（目标） | 状态 | 证据 |
| --- | --- | --- | --- |
| REQ-01 | Step 0-6 全流程（笔记 0-5 口径逐字一致，流程规格见 §3.1） | ✅ 完成 | RUN-05/07/08/19 |
| REQ-02 | 模块化架构（run_pipeline/config/logger/runner/resource/scanner/dingtalk/report/alerts/check_design + modules/8 软件模块） | ✅ 完成 | pipeline/ 目录 |
| REQ-03 | 资源自适应：探测(cgroup/WSL2/affinity)→分档推导→内存收紧→钳位；禁写死参数；快速失败；计划透明（推导表见 §5.3） | ✅ 完成 | RUN-10~14；DEC-03 |
| REQ-04 | 幂等断点续跑：产物存在且非空即 SKIP，只补缺失（含 .tbi/.metrics/.table；各步幂等键见 §3.1） | ✅ 完成 | RUN-09（31 SKIP 结果一致） |
| REQ-05 | 批次独立 + 多批次失败隔离；任一批失败退出码非零（退出码语义见 §7.3） | ✅ 完成 | RUN-07/08 |
| REQ-06 | 输入校验：两种布局自动识别；R1/R2 不匹配/单端/0字节/命名无效标记跳过；md5sum.txt 并行校验（失败 P0 阻断） | ✅ 完成 | tests/test_scanner.py |
| REQ-07 | 安全边界：0_raw_data/back 只读；产物仅写 results/ 与 Output/ 交付目录 | ✅ 完成 | 全程无违规写入 |
| REQ-08 | 钉钉通知：markdown 官方子集、启动/全步骤里程碑/完成/失败、**P0/P1/P2 三级分级（§5.2）**、失败降级不中断 | ✅ 完成 | RUN-16/33；DEC-08/21 |
| REQ-10 | CLI 参数全集（§7.2）+ `--version` 版本溯源 | ✅ 完成 | run_pipeline.py argparse |
| REQ-11 | 日志自动落盘（tee 接管 stdio，无需 shell 重定向）+ run_summary.json 结构化（含 pipeline_version） | ✅ 完成 | RUN-16；DEC-09 |
| REQ-12 | 测试与验收体系（§9）：纯标准库回归集 + TH 防漂移校验 + CI，作为改动/迁移/重建的验收闭环 | ✅ 完成 | RUN-20/30；tests/ |
| REQ-13 | 最终交付独立目录 Output/<批次>_<日期>/：VCF+tbi、MultiQC 汇总报告、md5sum.txt（可 -c 校验）、MANIFEST.tsv、README（含版本行），幂等导出；Output/INDEX.md 累积索引 | ✅ 完成 | RUN-18/27/29 |
| REQ-14 | 多次运行隔离：results/<批次>_<执行日期>/，同日同目录续跑、跨日新目录；**执行日期在进程启动时取一次（跨 0 点不切换目录）** | ✅ 完成 | RUN-19；DEC-01；跨午夜锚 |
| REQ-15 | 小测试用例集（纯标准库秒级回归，覆盖历史踩坑），供改动/迁移快速验证 | ✅ 完成 | RUN-20；tests/ |
| REQ-16 | Step6 on-target 告警阈值 8%（小 panel 正常 8-10%，TH-14） | ✅ 完成 | config.ON_TARGET_P1；笔记 5-5 |

> 新增需求：人追加 REQ-17… 行并填目标列，状态列留空由机器回写。

## 3. 系统架构 <!-- MACHINE 维护，人可提议变更（记 DEC） -->

```
$WORK/
├── singularity/       容器镜像目录（pipeline 同级，固定位置；GWAS_SIF_DIR 可重定向）
│   ├── fastqc_0.12.1.sif  fastp_1.3.6.sif  bwa-mem2_2.3.sif  samtools_1.24.sif
│   └── gatk_4.6.2.0.sif   bcftools_1.24.sif mosdepth_0.3.14.sif multiqc_1.35.sif
└── pipeline/
    ├── .env               环境参数文件（DEC-18；钉钉 webhook 等密钥只存于此，600，.env.example 为模板）
    ├── run_pipeline.py   主控：CLI/编排（process_batch→BatchCtx.step0~6 步骤方法）/样本级并行/汇总/通知/交付导出
    ├── config.py         ★ 阈值与路径唯一实现源（TH 表镜像于此；PIPELINE_VERSION；启动加载 .env）
    ├── logger.py         双通道日志 + capture_stdio tee（缓冲模式延迟落盘）
    ├── runner.py         容器命令封装：幂等 SKIP/实时输出/超时 Timer kill/返回码/产物校验/统计超时
    ├── resource.py       ★ 资源探测与规划（cpu=min(os/affinity/cgroup)；mem=min(MemAvailable/cgroup)）
    ├── scanner.py        两种布局扫描 + md5（样本名按布局推导）+ Lane 合并 + samples.tsv
    ├── dingtalk.py       markdown 通知 + 三重规范化（表格降级/\n→\n\n/截断）+ 未配置降级
    ├── alerts.py         P0/P1/P2 三级检查 + 里程碑消息构造
    ├── report.py         run_report/run_summary
    ├── check_design.py   TH 表 ↔ config.py + 版本双源一致性校验（防漂移）
    ├── tests/            回归用例（踩坑映射见 tests/README.md）
    └── modules/          fastqc fastp bwa_mem2 samtools gatk bcftools mosdepth multiqc
```

数据流：`0_raw_data/<批次>/ → results/<批次>_<日期>/{fastq_merged→fastq_clean→bam→gvcf→cohort→matrix→per_sample_vcf→qc} → Output/<批次>_<日期>/`。

参考文件（`$WORK/reference/`，依赖缺失即 P0 阻断）：`genome/genome.fa(+.fai/.dict)`
（hg38 + bwa-mem2 索引）、`targets.bed`（运行时生成 sorted.bed/interval_list）、
`Homo_sapiens_assembly38.dbsnp138.vcf.gz`、`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`。

### 3.1 Step 0-6 流程规格（每步：操作 → 产物〔幂等键〕→ 分级检查）

| Step | 操作（容器工具与关键语义） | 产物〔幂等键〕 | 分级检查 |
| --- | --- | --- | --- |
| **0 清点与合并** | 扫描两种布局（Illumina 平铺 `<样本>_S#_L###_R[12]_001.fastq.gz` / 外送 `<样本>/<样本>_R[12].fastq.gz`）→ `--samples` 白名单过滤 → **样本名白名单 `[A-Za-z0-9_.-]`（DEC-19）** → md5 并行校验 → 磁盘/依赖预检 → 多 Lane `cat` 合并（并行）→ samples.tsv | `fastq_merged/<样本>_R{1,2}.fastq.gz`、`samples.tsv` | P0：样本名非法/依赖缺失/磁盘不足/md5 损坏（统一 raise 中断批次）；P1：合并失败（样本终止）；启动通知含全部预检结论 |
| **1 QC+修剪** | fastqc(raw) → fastp（`-l 36`、adapter 自动检测，线程=fastp_threads）→ fastqc(trim) + **Adapter raw→trim 复检** | `fastq_clean/<样本>_R{1,2}.fastq.gz`〔fastp html/json〕、`qc/fastqc_raw|fastqc_trim|fastp`〔fastqc zip〕 | P1：reads<1M（TH-33）；P2：保留率<80/Q30<85（TH-03/04）、NTC reads 占比>1%（TH-34） |
| **2 比对** | `bwa-mem2 mem -K 100000000 -Y -R '@RG…SM:样本' | samtools sort -@T -m SORT_MEM`（管道流式）→ `samtools index` → flagstat/stats | `bam/<样本>/<样本>.sort.bam(+.bai)`〔sort.bam+.bai〕、`qc/flagstat|stats`〔对应文件〕 | P1：mapped<90 QC 口径（TH-05，报错不中断）；P2：mapped<95/pp<85（TH-06/07） |
| **3 去重** | `gatk MarkDuplicates` → markdup.bam + metrics → flagstat/stats（去重前 duplicates 行不采集） | `bam/<样本>/<样本>.markdup.bam`〔markdup.bam+metrics〕 | P2：dup>30%（TH-09） |
| **4 BQSR** | `gatk BaseRecalibrator`（dbsnp+Mills）→ recal.table → `ApplyBQSR` → **BQSR 前后 flagstat 逐行一致断言**（不一致该样本失败隔离） | `bam/<样本>/<样本>.recal.table`、`<样本>.markdup.BQSR.bam`〔BQSR bam+table〕 | P2：recal M 事件观测数<1e5（TH-36，校准不可信） |
| **5 变异检测** | BedToIntervalList 准备 → 每样本 `gatk HaplotypeCaller -ERC GVCF`（靶区±100bp，样本级并行）→ gvcf.list（sorted，DEC-05）→ 批次内 `CombineGVCFs→GenotypeGVCFs`（串行）→ `bcftools norm -m -any` → SelectVariants SNP/INDEL → VariantFiltration 硬过滤（SNP：QD2/QUAL30/SOR3/FS60/MQ40/MQRankSum/ReadPosRankSum；INDEL：QD2/QUAL30/SOR10/FS200/ReadPosRankSum）→ concat → **FILTER 列出现 `.` 判错** → PASS 提取 → 双口径矩阵导出（`query -R`）+ 每样本 hardfiltered/PASS 拆分（`view -s`） | `gvcf/<样本>.g.vcf.gz`、`cohort/cohort.{raw→split→snp/indel→hardfiltered→PASS}.vcf.gz`、`matrix/genotype_matrix.tsv`、`matrix/genotype_detail_PASS.tsv`、`per_sample_vcf/<样本>.{hardfiltered,PASS}.vcf.gz` | cohort 级失败 raise 中断批次；P2：对账·数量（矩阵行数 vs `view -R` 计数，DEC-11 同口径）、对账·新鲜度（关键 VCF mtime<启动） |
| **6 汇总与交付** | mosdepth×2（markdup/BQSR bam）→ `CollectHsMetrics`（BQSR bam）→ **矩阵 `./.` 裁决**（mosdepth bqsr regions 深度 DP≥TH-15 改判 0/0）→ 每样本 PASS VCF 重建（`view -s` + GT 替换，其余字段原样）→ **MultiQC 最终报告（此时全部流程结束、QC 齐全）** → 交付导出 Output/ | `qc/mosdepth|hs metrics|multiqc`、`matrix/genotype_matrix.adjudicated.tsv`、`per_sample_vcf/<样本>.PASS.adjudicated.vcf.gz(+.tbi)`、`Output/<批次>_<日期>/` 全套 | P1：PASS 裁决 VCF 0 记录（交付为空）、NTC 靶深>10×（TH-21）；P2：on-target<8/depth<50/20X<95/Ti-Tv<2.0/call rate<95（TH-14/11/13/22/23）、批次深度 CV>0.5（TH-35） |
| 每步之后 | `disk_guard`：Step 间磁盘复查 | — | P0：剩余<DISK_MIN_FREE_GB 立即终止批次（TH-24） |

并行模型：Step 1-4 与 Step 6 的样本级任务按 `plan.workers` 线程池并行；Step 5 的
cohort 级与 MultiQC 串行。样本失败即隔离（记入 failed，退出后续步骤，P1 通知）。

## 4. 关键设计决策（DEC）<!-- MACHINE 回写；人可新增决策行（含动机） -->

| 编号 | 决策 | 动机/证据 |
| --- | --- | --- |
| DEC-01 | 批次完全独立；输出 `results/<批次>_<执行日期>/`；**执行日期取自进程启动时间戳一次性派生并传参——运行跨 0 点不切换目录**（同日重跑同目录幂等续跑；跨日中断后重跑进入新日期目录全新计算，属设计行为） | 多次运行隔离 + 同日幂等续跑（RUN-19）；tests/test_misc.py 跨午夜锚锁定 |
| DEC-02 | 交付独立目录 Output/ + 标准 md5sum.txt/MANIFEST；MultiQC 最终报告在**全部流程（含矩阵裁决与每样本 VCF 重建）结束、QC 文件生成完全后**才生成，随后连报告一起交付导出（时序锚锁定） | 与输入侧 md5 约定对称，接收方可 `md5sum -c`（RUN-18/27，用户要求） |
| DEC-03 | 资源模型：索引 17G + 运行时；sort=剩余预算×0.25 钳位 128M-2G（整数 MB）；GATK 1-8g；cohort=可用 60% 钳位 8-32g；峰值>可用→SystemExit | `sort -m 2.0G`→2 字节坑（RUN-02）；cgroup 实测下限 24-26G（RUN-11~14） |
| DEC-04 | 命令/参数语义以笔记为准；资源参数规划化；差异记录于 pipeline/README §8 | 冲突裁决顺序 |
| DEC-05 | 样本列序全链 sorted（calling/gvcf.list/裁决矩阵/重建） | GT 列互换 bug：曾致两样本基因型对调（RUN-07→08） |
| DEC-08 | 钉钉：标题 `[GWAS][P级别]`（关键词实测只校验正文且大小写敏感，正文自动兜底小写 gwas）；表格自动降级列表；DINGTALK_MILESTONES=0 只发异常级 | 官方文档实测（RUN-16 前探测实验） |
| DEC-09 | 日志 tee：capture_stdio **缓冲模式**接管 stdout/stderr，批次发现后落盘到批次目录（单批次=该批次 logs/；多批次=首个批次 logs/）；dry-run/启动即退缓冲丢弃零落盘；results/ 根禁止散文件 | 用户要求"日志自动输出"（RUN-16/24） |
| DEC-10 | run_summary.json **逐批写入该批次目录**（快照含截至该批的 run+资源计划+各批次状态+pipeline_version）；results/ 根不写全局文件；dry-run 不写 | 曾写全局 run_summary 且单批次也落全局（RUN-24 收紧） |
| DEC-11 | 对账同口径：矩阵(query -R) 与核对(view -R) 一致；新鲜度=关键 VCF mtime<本次启动 | -R/-T 对跨界 indel 取舍不同（RUN-19：7399 vs 7390） |
| DEC-12 | 测试集 = 踩坑回归集；迁移验证顺序见 §9.2 | 已抓出 flagstat total 行、runner 静默超时 2 个潜伏 bug（RUN-20） |
| DEC-13 | NTC 纳入 QC（含 mosdepth/HsMetrics）、排除联合检测；对照样本阈值只记录不判级 | NTC 近零覆盖属预期 |
| DEC-14 | 日志严重级区分：流程失败/断言不一致=ERROR；数据质量违例=WARN | "无遗留 ERROR"验收口径需要语义干净的 ERROR |
| DEC-15 | on-target 阈值取 8%（TH-14） | 60% 误用了宽口径，小 panel 正常 on-target≈8-10%（笔记 5-5） |
| DEC-16 | 容器运行时在 Runner 初始化时自解析为**绝对路径**（shutil.which→标准目录探测），子进程 PATH 兜底补齐标准系统目录 | PATH 受限环境曾致 `singularity: not found`(127)（RUN-22） |
| DEC-17 | **流程只在宿主 Python 运行**：启动卫兵检测 `/.singularity.d` 或 `SINGULARITY_*`，容器内立即退出并指引用 /usr/bin/python3；流程纯标准库无需 mamba/conda | 嵌套容器运行全部 127 且 TZ=UTC（RUN-22/23） |
| DEC-18 | **环境参数统一 .env**：钉钉 webhook（含 access_token）等环境参数只存 `pipeline/.env`（模板 .env.example；`GWAS_ENV_FILE` 可改址），config.py 启动解析并入（setdefault），**三源优先级=进程环境变量 > .env > 内置默认**；webhook 未配置→通知静默跳过（启动 WARN+配置指引） | 密钥硬编码随代码泄露且换环境必改源码（RUN-26）；.env 600 |
| DEC-19 | **样本名白名单 P0 阻断**：含 `[A-Za-z0-9_.-]` 外字符 → P0 异常项入启动通知后直接中断批次（dry-run 同样触发） | 样本名进入 shell 命令拼接与 @RG 头，非常规字符属注入面（RUN-31→32） |
| DEC-20 | **process_batch 步骤方法化**：Step 0-6 拆为 `BatchCtx.step0_scan~step6_summary_delivery` 七方法，编排层仅 ~40 行；跨步骤状态挂 ctx（bdata/merged/excluded/cohort_stats/notify_on/run_date/_t0） | 原单函数 ~700 行难读难测（RUN-31）；MultiQC 时序锚锚定 step6 方法源码 |
| DEC-21 | **三级分级体系**：P0=阻断级（严重影响分析→raise 中断批次，退出码 1）；P1=严重（执行失败样本隔离/NTC 污染/mapped<90，报错不中断）；P2=质量提示（阈值越界/口径存疑，只记录）。质量类一律不中断 | 用户决定：P0 应为严重影响分析的阻断级；质量只报错不中断（RUN-33） |
| DEC-22 | **审计项落地**：P0 Step 间磁盘复查（disk_guard）；P1 reads<1M、PASS 裁决 VCF 空结果；P2 NTC reads 占比/批次深度 CV/recal 观测数（TH-33~36） | RUN-33 审计清单经用户圈选实施（RUN-34） |

## 5. 统一口径与阈值总表（TH = config.py 镜像）<!-- HUMAN 可改值；MACHINE 同步 config 后过校验 -->

> **单一事实来源是 `pipeline/config.py`**；下表改动必须同步 config.py（或让机器代改），
> `python3 pipeline/check_design.py` 校验两处一致，`run_tests.sh` 末尾自动执行。
> TH-17~20 编号已废弃不复用；TH-21 起编号保持不变以维持引用稳定。

### 5.1 阈值表

| 编号 | config 键 | 值 | 语义 |
| --- | --- | --- | --- |
| TH-01 | FASTP_LENGTH_REQUIRED | 36 | fastp 最短读长（笔记） |
| TH-02 | FASTP_RETENTION_WARN | 95.0 | 保留率<95% 提示 |
| TH-03 | FASTP_RETENTION_P1 | 80.0 | 保留率<80% → P2 |
| TH-04 | FASTP_Q30_P1 | 85.0 | Q30<85% → P2 |
| TH-05 | MAPPED_MIN_PCT | 90.0 | mapped<90% → P1 报错不中断（QC 口径） |
| TH-06 | MAPPED_NOTIFY_P1 | 95.0 | mapped<95% → P2（通知从严） |
| TH-07 | PROPER_PAIR_P1 | 85.0 | properly paired<85% → P2 |
| TH-08 | DUP_WARN_PCT | 30.0 | dup>30% 且 ELS 小 → 复杂度告警 |
| TH-09 | DUP_P1 | 30.0 | dup>30% → P2 |
| TH-10 | MEAN_COV_MIN | 50.0 | MEAN 靶深度<50× WARN |
| TH-11 | MEAN_DEPTH_P1 | 50.0 | mean depth<50× → P2 |
| TH-12 | PCT_20X_MIN | 90.0 | 20X<90% WARN |
| TH-13 | PCT_20X_P1 | 95.0 | ≥20x 靶比例<95% → P2 |
| TH-14 | ON_TARGET_P1 | 8.0 | on-target<8% → P2（小 panel 正常 8-10%，笔记 5-5） |
| TH-15 | DP_MIN | 20 | `./.` 裁决深度阈值 |
| TH-16 | HC_INTERVAL_PADDING | 100 | HC 靶区外扩 bp（笔记） |
| TH-21 | NTC_DEPTH_P0 | 10.0 | NTC 靶深>10× → P1 污染（报错不中断） |
| TH-22 | TITV_P1 | 2.0 | Ti/Tv<2.0 → P2 |
| TH-23 | CALL_RATE_P1 | 95.0 | call rate<95% → P2 |
| TH-24 | DISK_MIN_FREE_GB | 200 | 磁盘剩余 P0 阻断线（开跑前+Step 间复查） |
| TH-25 | DISK_PER_SAMPLE_GB | 75 | 每样本磁盘估算 |
| TH-26 | BWA_INDEX_MEM_GB | 17.0 | bwa-mem2 索引 mmap 常驻 |
| TH-27 | SORT_MEM_MIN | 128 | sort 每线程下限 MB |
| TH-28 | SORT_MEM_MAX | 2048 | sort 每线程上限 MB |
| TH-29 | GATK_MEM_MIN_GB | 1 | GATK -Xmx 下限 g |
| TH-30 | GATK_MEM_MAX_GB | 8 | GATK -Xmx 上限 g |
| TH-31 | COHORT_MEM_MIN_GB | 8 | cohort -Xmx 下限 g |
| TH-32 | COHORT_MEM_MAX_GB | 32 | cohort -Xmx 上限 g |
| TH-33 | READS_MIN | 1000000 | 样本 reads 绝对量 <1M → P1（上样不足） |
| TH-34 | NTC_READS_PCT_P2 | 1.0 | NTC reads 占批次中位样本 %>1% → P2 |
| TH-35 | DEPTH_CV_P2 | 0.5 | 批次内 mean depth 变异系数 CV>0.5 → P2 |
| TH-36 | RECAL_OBS_MIN_P2 | 100000 | BQSR recal M 事件观测数 <1e5 → P2（校准不可信） |

### 5.2 分级告警体系（DEC-21/22，三级）

**P0=阻断级**（严重影响分析 → 直接中断批次，退出码 1，多批次运行其余批次照常隔离）；
**P1=严重**（执行失败按样本隔离 / 严重数据异常，报错不中断）；
**P2=质量提示**（阈值越界 / 口径存疑，只记录，不影响运行）。

| 时机 | 规则 | 级别 |
| --- | --- | --- |
| 开跑前 · 依赖文件 | sif / genome(.fai/.dict) / targets / known-sites 缺失或 0 字节 → 中断批次 | P0 |
| 开跑前 · 样本名 | 含非常规字符（A-Za-z0-9_.- 外，注入面）→ 中断批次（DEC-19） | P0 |
| 开跑前 · 磁盘空间 | 剩余 < TH-24（每样本估算 TH-25）→ 中断批次 | P0 |
| 开跑前 · md5 | 输入校验失败（数据损坏）→ 中断批次 | P0 |
| Step 1-6 间 · 磁盘复查 | Step 结束时剩余 < TH-24 → 终止批次 | P0 |
| Step 0 · 合并失败 | cat 非 0（样本终止，其余照常） | P1 |
| Step 1-6 · 执行失败 | 任一步 exit≠0 / 产物 0 字节（样本级隔离，不中断批次） | P1 |
| Step 1 · 上样量 | 样本 reads 绝对量 < TH-33 | P1 |
| Step 2 · QC 口径 | mapped < TH-05（报错不中断） | P1 |
| Step 6 · 空交付 | per-sample PASS 裁决 VCF 0 条记录 | P1 |
| Step 6 · 对照污染 | NTC 靶区深度 > TH-21 | P1 |
| 全流程 · 失败样本 | 任一样本在某步失败（汇总报错） | P1 |
| Step 1 · 修剪 | fastp 保留率 <TH-03 或 Q30 <TH-04（TH-02~03 区间为提示行） | P2 |
| Step 1 · NTC reads | NTC reads 占批次中位样本 > TH-34 | P2 |
| Step 2 · 比对 | mapped <TH-06 或 properly paired <TH-07 | P2 |
| Step 3 · 去重 | 重复率 > TH-09（建库复杂度告急） | P2 |
| Step 4 · 校准可信 | BQSR recal M 事件观测数 < TH-36 | P2 |
| Step 5 · 对账·数量 | 矩阵行数 ≠ 靶区记录数（view -R 同口径，DEC-11） | P2 |
| Step 5 · 对账·新鲜度 | 关键 VCF mtime < 本次启动（断点续跑复用旧产物） | P2 |
| Step 6 · 捕获/覆盖/口径 | on-target <TH-14 / depth <TH-11 / 20X <TH-13 / Ti-Tv <TH-22 / call rate <TH-23 | P2 |
| Step 6 · 深度离散 | 批次内 mean depth CV > TH-35（疑似混入异常样本） | P2 |

### 5.3 资源规划推导（DEC-03）

- 探测：CPU = min(os.cpu_count, sched_getaffinity, cgroup cpu.max 配额)；内存 = min(MemAvailable, cgroup memory.max)。
- 预留：2 核 + max(2G, 可用内存 5%)。
- 每样本线程 T 分档：可用核 ≥64→24；32-63→12；16-31→8；<16→4。
- workers = min(⌊可用核/T⌋, ⌊可用内存/单样本峰值⌋)；单样本峰值 = 17G 索引 + T×sort缓冲 + GATK + 1G 杂项。
- sort -m = 剩余预算×0.25，钳位 TH-27~28（整数+单位，如 `512M`）；GATK -Xmx 钳位 TH-29~30；cohort（串行）= 可用 60% 钳位 TH-31~32。
- 快速失败：workers=1 时峰值仍 > 可用 → 启动即 SystemExit。
- 档位：`auto`=实测（超过 100 线程/900G 照常用）；`low`=16 线程/20G 强制规划；`high`=100 线程/900G。
- 优先级：`--workers` > `--threads`/`--max-memory` > `--resource-profile` > 自动探测。

### 5.4 通知时机与模板

每批次：启动（样本数/输入体量/资源计划+预检结论）→ Step 0-6 每步里程碑 → 完成/失败；
多批次另有总览。`DINGTALK_MILESTONES=0` 只发异常级（P0/P1）。里程碑模板：

```
[GWAS][P1] 20260720批次 · Step 2 比对完成
样本: 4/4 成功 | Lane 合并 16/16
指标: mapped 98.7% | proper pair 94.2%
异常: [P1] L20260615001 mapped 91.3%（阈值 95%）← 需确认
产物: bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22
日志: tail -f logs/sample_L20260615001.log
```

## 6. 输出与目录规约 <!-- MACHINE -->

- **`results/` 根下只允许 `<批次>_<执行日期>/` 目录，无任何散文件/子目录**
- `results/<批次>_<执行日期>/`：samples.tsv、fastq_merged/fastq_clean/bam/gvcf/cohort/matrix/per_sample_vcf/qc、run_report.md、run_summary.json、logs/（执行日期启动时固定，跨 0 点不切换，DEC-01）
- `Output/<批次>_<执行日期>/`：`*.PASS.adjudicated.vcf.gz(+.tbi)`、`*multiqc_report.html`、md5sum.txt、MANIFEST.tsv、README.md（含版本行）；`Output/INDEX.md` 跨批次**累积**索引（扫描全部历史交付目录 ∪ 本次运行）
- 日志层级：`<批次>/logs/run_<ts>.log`（入口 tee，多批次写首个批次目录）→ `logs/pipeline_<ts>.log` → `sample_<样本>.log`；dry-run 全部不落盘
- run_summary.json：逐批写批次目录快照（DEC-10）；`--out` 显式覆盖时写该目录

```
results/260422_20260914/                     Output/260422_20260914/
├── samples.tsv                               ├── INDEX.md（在 Output/ 根，累积）
├── fastq_merged/  fastq_clean/               └── NA12878.PASS.adjudicated.vcf.gz (+.tbi)
├── bam/<样本>/ （sort/markdup/BQSR+metrics+table） ├── GWAS-Panel-下机数据-QC_multiqc_report.html
├── gvcf/<样本>.g.vcf.gz                      ├── md5sum.txt（md5sum -c 可校验）
├── cohort/（raw→split→过滤→PASS）            ├── MANIFEST.tsv（样本/文件/大小/md5/记录数/来源）
├── matrix/（全口径+裁决矩阵+PASS 详情）       └── README.md（口径+版本行）
├── per_sample_vcf/（hardfiltered/PASS/裁决）
├── qc/（fastqc_raw|trim、fastp、flagstat、stats、
│        hsmetrics、mosdepth、bcftools_stats、multiqc）
├── run_report.md  run_summary.json
└── logs/（run_* / pipeline_* / sample_*）
```

## 7. 配置与接口 <!-- MACHINE -->

### 7.1 环境参数（.env，DEC-18）

**三源优先级：进程环境变量 > `pipeline/.env` > config.py 内置默认**。
语法：`KEY=VALUE`；`#` 整行注释，裸值支持行内 ` #` 注释，引号值原样保留。

| 键 | 用途 | 默认（未配置时） |
| --- | --- | --- |
| `DINGTALK_WEBHOOK` | 钉钉机器人地址（**密钥，600 权限**） | 空 → 通知静默跳过（启动 WARN） |
| `DINGTALK_KEYWORD` / `DINGTALK_MILESTONES` | 关键词 / 里程碑开关（0=只发异常级） | `gwas` / `1` |
| `CONTAINER_RT` | 容器运行时 | `singularity` |
| `GWAS_ENV_FILE` | .env 文件位置 | `pipeline/.env` |
| `GWAS_RAW_DATA` / `GWAS_RESULTS` / `GWAS_DELIVERY_DIR` | 输入/结果/交付目录 | `$WORK/0_raw_data`、`$WORK/results`、`$WORK/Output` |
| `GWAS_REFERENCE_DIR` / `GWAS_SIF_DIR` | 参考文件/镜像目录重定向 | `$WORK/reference`、`$WORK/singularity` |
| `GWAS_DP_MIN` …（TH-15/24/25/33~36 同名键） | 阈值覆盖 | 见 §5.1 |

### 7.2 CLI 参数

| 参数 | 说明 | 默认 |
| --- | --- | --- |
| `--input <目录>` | 多批次输入根目录（相对名按 RAW_DATA_DIR 语境解析，不回落 $WORK） | `0_raw_data` |
| `--batch <批次名>` | 指定单批次（优先于 --input 遍历） | - |
| `--step <0-6>` | 执行到第几步（§3.1） | `6` |
| `--samples a,b` / `--exclude-samples NTC` | 样本白名单 / 联合检测排除对照（空串关闭） | 全部 / `NTC` |
| `--workers N` / `--threads N` / `--max-memory NG` | 覆盖资源规划（优先级见 §5.3） | 规划值 |
| `--resource-profile auto\|low\|high` | 资源档位（§5.3） | `auto` |
| `--serial` | 强制单样本串行 | - |
| `--dry-run` | 打印命令与资源计划不执行（零落盘，含 P0 预检暴露） | - |
| `--notify on\|off` / `--notify-test` | 通知开关 / 发送测试后退出 | `on` / - |
| `--version` | 显示版本（config.PIPELINE_VERSION，与本文档同步） | - |
| `--out <目录>` | 覆盖批次结果目录 | `results/<批次>_<日期>` |

### 7.3 退出码

- `0`：全部批次 success 或 skipped（无有效样本跳过不算失败）
- `1`：任一批 failed——P0 阻断（依赖/样本名/磁盘/md5/Step 间磁盘）、全部样本合并失败、
  无 calling 样本、步骤异常；容器内运行卫兵拦截（DEC-17）同样退出 1
- `2`：输入目录不存在 / argparse 参数错误

## 8. 已知边界与风险 <!-- MACHINE -->

1. **cgroup 硬上限**：bwa-mem2 17G 索引 + 运行时/页缓存，硬上限 <24GiB 会被 OOM-kill
   （实测修正 sort 分配后 24G 通过，32G 峰值 25.59GB）；裸机 20G 属临界。规划器对
   装不下的机器启动即快速失败。
2. **CPUQuota 在 WSL2 用户 slice 不生效**（cpu.max 缺失）：限核用 taskset（规划器已识别
   sched_getaffinity）；systemd 托管机 cpu.max 可被识别。
3. 思源笔记为活文档：手册只是快照，分析学口径冲突时以笔记原文（MCP 回查）为准。
4. 跨日中断重跑进入新日期目录从头计算（同日才续跑）：长时间批次建议白天启动或
   次日同日补跑（DEC-01 设计行为）。
5. 统计类命令统一超时 STATS_TIMEOUT_S=600（flagstat/stats/count/query/对账）；
   GATK/bwa 等长任务不设超时。

## 9. 测试与验收（重建验收标准）<!-- MACHINE -->

### 9.1 测试集构成（`pipeline/tests/`，纯标准库、无网络无容器，约 5 秒）

| 文件 | 覆盖 |
| --- | --- |
| test_runner.py | tool 拼装/binds、cpath 映射、rt 绝对路径自包含解析、PATH 兜底、幂等 SKIP、dry-run 不执行、超时 kill、统计命令 timeout=600 传递 |
| test_resource.py | sort -m 整数 MB、-Xmx 格式、low/high 档精确值、快速失败（mock 探测）、计划表、workers 覆盖 |
| test_scanner.py | 两种布局、无效输入（R1R2 不匹配/0 字节）、md5、合并（真实与 dry-run）、samples.tsv 列 |
| test_dingtalk.py | 表格降级、换行规范化、截断、发送成功落日志（mock urlopen） |
| test_alerts.py | P0/P1/P2 三级判定全集（含 TH-33~36 新检查）、最差级别、里程碑模板 |
| test_parsers.py | fastqc/fastp(json 真实结构)/markdup/hsmetrics/flagstat/stats/bcftools/mosdepth/recal 观测数 |
| test_variant_post.py | 矩阵 `./.` 裁决、重建只换 GT、GT 列显式映射、norm outputs 口径 |
| test_envfile.py | .env 解析语法、三源优先级、webhook 默认空、防回潮锚（源码无 access_token=）、降级 |
| test_misc.py | tee 自动落盘、目录命名、跨午夜锚、全流程 dry-run success+零落盘锚、样本名 P0 阻断锚、--input 覆盖锚、--version、MultiQC 时序锚（step6 方法）、交付导出（Output 命名/md5sum/MANIFEST/MultiQC extra/幂等）、INDEX 累积锚、disk_guard |

`./run_tests.sh` = unittest 全量 + `check_design.py`（TH↔config + 版本双源）；
CI（.github/workflows/ci.yml）在 py3.10/3.12 矩阵执行。

### 9.2 迁移/重建验证顺序（DEC-12）

1. `./run_tests.sh` 全绿（全部用例数见 §9.1 各文件，当前共 103）
2. `cp .env.example .env` 填 webhook → `python3 run_pipeline.py --notify-test`（连通性）
3. `python3 run_pipeline.py --dry-run --batch <小批次>`（容器/参考文件/路径与资源计划）
4. `python3 run_pipeline.py --resource-profile low --dry-run`（低配档口径）
5. 单批次冒烟（2 样本）→ 指标对照 RUN_HISTORY §三 → 全量

### 9.3 重建完成判据

103 用例 + check_design 全绿；`--version` 输出与本文档 version 一致；dry-run 零落盘；
单批次实跑 success 且 Step 6 交付目录含 VCF+tbi+MultiQC，`md5sum -c` 全过。

## 10. 变更日志（CHANGELOG）<!-- 人机共写：每方改动各记一行 -->

| 版本 | 日期 | 角色 | 变更 |
| --- | --- | --- | --- |
| 1.0.0 | 2026-09-14 | 人 | 依提示词确立目标与验收（REQ-01~12） |
| 1.1.0 | 2026-09-14 | 人 | 新增 REQ-13 交付目录 / REQ-14 多次运行目录隔离 / REQ-15 测试集 |
| 1.1.0 | 2026-09-14 | 机 | 完成实现与验证（RUN-01~20），建立本文档与 RUN_HISTORY、check_design 防漂移校验 |
| 2.0.0 | 2026-09-14 | 人 | 确认 REQ-14 执行日期启动时固定；on-target 阈值 60→8（REQ-16/TH-14，DEC-15）；移除内置比对（比对不属于常规运行） |
| 2.0.0 | 2026-09-14 | 机 | 删除 compare 模块与相关配置/告警/用例（75→69）；tests+check_design 全绿；RUN-21 |
| 2.0.1 | 2026-09-14 | 机 | 健壮性修复（DEC-16）：Runner 自解析运行时绝对路径 + PATH 兜底；测试 69→71 |
| 2.0.2 | 2026-09-14 | 机 | 根因修正（DEC-17）：嵌套容器卫兵拒绝并指引宿主 python3；测试 71→73 |
| 2.1.0 | 2026-09-14 | 人 | 目录规约收紧：results/ 根只允许批次_日期目录；运行日志归入批次目录 |
| 2.1.0 | 2026-09-14 | 机 | tee 缓冲模式延迟落盘；run_summary 逐批写批次目录、dry-run 不写；测试 73→76 |
| 2.1.1 | 2026-09-15 | 机 | Q30 解析修复（RUN-25）：fastp json 键名按真实产物校准，统一百分数口径 |
| 2.2.0 | 2026-09-16 | 人 | 环境参数不得硬编码，改 .env；singularity 位置不变 |
| 2.2.0 | 2026-09-16 | 机 | DEC-18：.env 三源加载、webhook 密钥迁移、未配置降级；测试 75→86（RUN-26） |
| 2.3.0 | 2026-09-16 | 人 | 交付目录 delivery→Output；MultiQC 报告随交付 |
| 2.3.0 | 2026-09-16 | 机 | DEC-02：MultiQC 移至全流程结束后、export_delivery 增 extra_files；修 3 处 dry-run 守卫缺失；测试 86→89（RUN-27） |
| 2.4.0 | 2026-09-16 | 人 | 文档清理：删除全部比对表述与废弃内容（删除项经确认） |
| 2.4.0 | 2026-09-16 | 机 | DESIGN 删比对数据节（节号重排）/REQ-09/DEC-06/07；RUN_HISTORY 清痕迹保留修复链；README/tests/config 同步；RUN-28 |
| 2.5.0 | 2026-09-16 | 人 | 依据实跑日志分析执行 4 项修复（norm 误报/对账文案/通知留痕/INDEX 累积） |
| 2.5.0 | 2026-09-16 | 机 | norm outputs 去 .tbi；-T→-R 文案；发送成功落 INFO；INDEX.md 累积合并；测试 89→92（RUN-29） |
| 2.6.0 | 2026-09-16 | 人 | design_doc 入库（PROMPT 与内部版 README 不入）；增加 GitHub Actions CI |
| 2.6.0 | 2026-09-16 | 机 | design_doc 迁入 pipeline/；ci.yml（py3.10/3.12）；测试 CI 兼容化（rt 自包含/low 档/哑依赖/资源 mock）；CI 两轮根因修复后 success；测试 92→96（RUN-30） |
| 2.7.0 | 2026-09-16 | 人 | 五项修改：样本名阻断（H1）/Step 方法化（M3）/版本溯源（M4）/统计超时（M5）/MIT 协议 |
| 2.7.0 | 2026-09-16 | 机 | DEC-19/20；PIPELINE_VERSION+--version+run_summary+交付 README 版本；STATS_TIMEOUT_S=600；LICENSE；测试 92→95（RUN-31） |
| 2.8.0 | 2026-09-16 | 人 | 样本名非法收紧为 P0 直接中断 |
| 2.8.0 | 2026-09-16 | 机 | DEC-19 更新为阻断；顺带修相对 --input 解析隐患（RAW_DATA_DIR 语境+锚）；测试 95→96（RUN-32，含验证事故记录） |
| 2.9.0 | 2026-09-16 | 人 | 分级体系修正：P0=严重影响并中断；原 P0/P1 降 P1/P2；质量只报错不中断；要求审计清单 |
| 2.9.0 | 2026-09-16 | 机 | DEC-21 三级体系（磁盘/md5 新增 P0 中断；mapped<90 不再判失败；alerts 全量降级）；顺带修 verify_md5 样本名推导（平铺布局 md5 失败样本从未被剔除的历史 bug）；RUN-33 交付审计清单 |
| 2.10.0 | 2026-09-16 | 人 | 实施 RUN-33 审计项：P0 磁盘复查 / P1 reads+空VCF / P2 NTC reads+CV+recal |
| 2.10.0 | 2026-09-16 | 机 | DEC-22：disk_guard/四新检查/TH-33~36；CI 修复 disk_guard 路径回溯；测试 96→103（RUN-34） |
| 2.11.0 | 2026-09-16 | 人 | 设计文档合并为**自含重建规格**：当前实际/修改情况/关键点全部并入；废弃目标删除（正文仅陈述现行为，沿革归 CHANGELOG）；目标=可据本文档从头重建项目 |
| 2.11.0 | 2026-09-16 | 机 | 新增 §3.1 Step 0-6 流程规格表（操作/幂等键/分级）、§5.2 分级告警体系全表、§5.3 资源推导、§7 配置与接口（.env/CLI/退出码）、§9 测试与验收（含 9.3 重建完成判据）、§6 目录树；REQ-08/12 更新为现口径；DEC 表按编号重排；TH/REQ/DEC 编号锚全部稳定不变；RUN-35 |

## 11. 证据索引 <!-- MACHINE -->

| 内容 | 路径 |
| --- | --- |
| 运行台账（每轮） | pipeline/design_doc/RUN_HISTORY.md |
| 验收运行（0_raw_data_test） | results/260422_20260914/、results/260422_20260916/ 等（运行日志在各批次 logs/） |
| 交付 | Output/<批次>_<日期>/ + INDEX.md（累积） |
| 测试集与 CI | pipeline/tests/（103 用例）+ run_tests.sh + .github/workflows/ci.yml |
| 使用说明/与笔记差异 | pipeline/README.md |
| 环境参数 | pipeline/.env（密钥，600）+ pipeline/.env.example（模板） |
| 命令参考快照 | pipeline/design_doc/notes_code_reference.md |
