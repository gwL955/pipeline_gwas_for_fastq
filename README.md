# GWAS 靶向测序 FASTQ→VCF 通用流程（pipeline/）

人类 GWAS 靶向测序全流程控制程序：**FASTQ → 质控修剪 → 比对 → 去重 → BQSR 校准 →
变异检测（gVCF 联合分型 → 硬过滤 → VCF）→ 测序质量汇总**。
Python3 **仅标准库**；所有生信工具经 singularity(apptainer) 容器执行；
通用、模块化、可断点续跑、资源自适应（16 线程/20G ～ 100 线程/900G 及以上）、带钉钉机器人通知。

命令与参数以思源笔记《下机数据流程》0-5 为准（备查手册
`design_doc/notes_code_reference.md`，已回查笔记原文核对一致）；
工程化约定以 `design_doc/DESIGN.md` 为准（`PROMPT_GWAS_pipeline.md` 为历史提示词，不入仓库）。

---

## 1. 快速开始

> **运行环境（重要）**：必须用**宿主系统 Python**（`/usr/bin/python3`）运行。
> 若 `python` 是 singularity 容器别名（mamba 解释器）——容器内看不到宿主的 singularity，
> 嵌套运行会全部 `singularity: not found`(127) 且时区为 UTC；程序启动会自检并拒绝。
> 本流程纯标准库，不依赖任何第三方/conda 包。

```bash
cd $WORK/pipeline    # $WORK = 工作区根目录（pipeline/ singularity/ 0_raw_data/ results/ 所在处）

# ① 干跑：打印全流程命令与资源计划，不执行
python3 run_pipeline.py --dry-run --batch 260422

# ② 冒烟批次（Step 0-6）
python3 run_pipeline.py --batch 260422

# ③ 全量：遍历 0_raw_data 下全部批次，逐批独立执行
#    日志全自动落盘（tee 进 results/<批次>_<日期>/logs/run_<时间戳>.log），无需 shell
#    重定向；实跑控制台在打印日志路径提示行后即静默（v2.13.0），nohup 仅用于后台
#    防断线，nohup.out 不会再堆积运行日志
nohup python3 run_pipeline.py --input 0_raw_data &

# ④ 低配档核对 / 钉钉测试
python3 run_pipeline.py --resource-profile low --dry-run
python3 run_pipeline.py --notify-test
```

中断后**直接重跑同一命令**即可续跑：每个产物（含 .tbi/.metrics/.table）落盘前检查
存在且非空，已存在即 `[SKIP]`，只补缺失部分。

## 1.1 日志说明（自动输出，无需重定向）

程序启动即接管自身 stdout/stderr（tee 镜像），包括资源计划表、快速失败报错与任何
未捕获输出，全部自动落盘。**实跑（v2.13.0 起）取消控制台外显**：控制台只显示启动段
（资源计划、参数、批次清单），到"运行日志: <路径>"提示行为止——此后全量日志只进
run_<ts>.log，`nohup … &` 后台运行不再向 nohup.out 倾倒（启动错误如输入目录不存在
发生在外显期，仍直接可见）；**dry-run 保持控制台全程可见**且零落盘。
**`results/` 根下只有 `<批次>_<执行日期>/` 目录**，所有日志归属批次目录：

- `results/<批次>_<日期>/logs/run_<ts>.log`：本次运行的入口日志（tee 全量镜像，
  含资源计划与主控输出；多批次运行写在首个批次目录下）
- `results/<批次>_<日期>/logs/pipeline_<ts>.log`：批次主日志
- `results/<批次>_<日期>/logs/sample_<样本>.log`：每样本独立日志（时间戳+耗时）
- `results/<批次>_<日期>/run_summary.json`：运行清单（主机/命令行/资源计划/
  各批次状态与指标；逐批增量快照，最后完成的批次文件即全貌）

## 1.2 环境配置（.env，v2.2.0 起）

**环境参数（钉钉企业机器人凭证、目录覆盖等）不写在代码里**，统一放
`pipeline/.env`（与代码同目录，模板 `pipeline/.env.example`；
`singularity/` 镜像目录仍在 pipeline 同级，属目录约定不走 `.env`）。
程序启动时自动加载，**优先级：进程环境变量 > `.env` > `config.py` 内置默认**。

```bash
cp .env.example .env    # 在 pipeline/ 下；填入 DINGTALK_CLIENT_ID 等企业机器人凭证
python3 run_pipeline.py --notify-test    # 验证钉钉链路（markdown + 文件）
```

| 键 | 用途 | 默认（未配置时） |
| --- | --- | --- |
| `DINGTALK_CLIENT_ID` / `DINGTALK_CLIENT_SECRET` | 企业内部机器人凭证（appKey/appSecret，**密钥勿外传**） | 空 → 通知整体静默跳过（启动日志提示） |
| `DINGTALK_CONVERSATION_ID` / `DINGTALK_ROBOT_CODE` | 群 openConversationId / 机器人编码（空=复用 CLIENT_ID） | 空 |
| `DINGTALK_MILESTONES` | `0` 时里程碑只发 P0/P1，OK 级静默（webhook 时代 WEBHOOK/KEYWORD 键已停用） | `1` |
| `GWAS_ARCHIVE_DIR` | 归档脚本 7z 输出目录 | `$WORK/archive` |
| `CONTAINER_RT` | 容器运行时（本机 singularity 实为 apptainer 别名） | `singularity` |
| `GWAS_ENV_FILE` | `.env` 文件本身的位置 | `pipeline/.env` |
| `GWAS_RAW_DATA` / `GWAS_RESULTS` / `GWAS_DELIVERY_DIR` | 输入/结果/交付目录覆盖 | 基于 `$WORK` 推导 |
| `GWAS_REFERENCE_DIR` / `GWAS_SIF_DIR` | 参考文件/容器镜像目录重定向（测试与多工作区部署用） | `$WORK/reference`、`$WORK/singularity` |
| `GWAS_TARGETS_BED` | 靶区 bed 改址（私密文件绝不入仓库；只指 bed 本体，派生 sorted.bed/interval_list 随其同目录生成） | `$WORK/reference/targets.bed` |
| `GWAS_DP_MIN` / `GWAS_DISK_MIN_FREE_GB` / `GWAS_DISK_PER_SAMPLE_GB` | 可覆盖阈值（TH-15/24/25） | 20 / 200 / 75 |

语法：`KEY=VALUE`；`#` 整行注释，裸值支持行内 ` #` 注释，引号值原样保留。
镜像目录**固定**为 `$WORK/singularity/`（目录约定，不走 `.env`）。

## 2. 输入与安全边界

- 输入：`0_raw_data/<批次>/`，批次是最小分析单元，**各批次完全独立**
  （独立扫描/QC/比对/去重/BQSR/批次内联合分型/过滤，结果隔离在
  `results/<批次>_<执行日期>/`，严禁跨批次合并样本或 gVCF）。三种布局自动识别：
  1. Illumina 平铺多 Lane：`260422/NA12878_S46_L001_R1_001.fastq.gz`
     （样本名 = 去掉 `_S\d+_L\d+_R[12]_001.fastq.gz` 的前缀）
  2. 外送子目录：`20260720/<样本名>/<样本名>_R1.fastq.gz`（样本名 = 子目录名）
  3. 外送平铺（v2.12.0/RUN-36）：`260917/<样本名>_R1.fastq.gz` 直接放批次目录
     （样本名 = 去掉 `_R[12].fastq.gz` 的前缀）
  每种布局一个独立识别器函数（`scanner.LAYOUT_SCANNERS`，DEC-23），
  依序应用、先认者优先，新增输入格式只需追加函数
- 批次内 `md5sum.txt` 存在时先做并行 md5 校验，失败样本终止分析并进入通知
- **只读**：`0_raw_data/`、`back/` 一律只读；一切分析产物写入 `results/<批次>_<执行日期>/`；
  交付文件仅由导出步骤写入 `Output/<批次>_<执行日期>/`（v2.3.0 前为 `delivery/`）
- 输入校验：R1/R2 文件数不一致、单端、0 字节、命名不匹配 → 该样本标记无效并跳过；
  **批次名/样本名含中文或非常规字符（DEC-19）→ P0 级错误直接中断该批次分析**
  （注入面 + 中文路径在 singularity 容器内因 locale 报错；多批次运行其余批次照常）
  （列入日志/通知/报告，不影响其余样本）；批次内无有效样本 → 跳过该批次并告警（退出码 0）
- NTC 等对照样本：纳入 QC，默认排除出联合变异检测（`--exclude-samples` 可改）
- 多批次失败隔离：某批失败只记录该批状态并通知，不阻断其余批次；
  任一批次失败（含部分样本失败）最终退出码非零

## 3. 输出目录命名（多次运行隔离）与交付目录

**批次结果目录 = `results/<批次名>_<执行日期>/`**（如 `results/260422_20260914/`）：
每次运行的产物相互隔离、可追溯；**同一天重跑命中同一目录 → 幂等断点续跑**（已有产物
`[SKIP]` 只补缺失），跨日运行自动开新目录全新计算。执行日期在启动时取一次（跨午夜
不切换）。`--out` 可显式覆盖。交付目录同构：`Output/<批次名>_<执行日期>/`。

分析中间产物（可整目录删除后重算）：

```
results/260422_20260914/
├── samples.tsv                     # 样本名/R1/R2 + S 号与 Lane 追溯列
├── fastq_merged/  fastq_clean/     # Step 0 多 Lane cat 合并 / Step 1 fastp
├── bam/<样本>/                     # sort/markdup/BQSR bam + metrics + recal.table
├── gvcf/<样本>.g.vcf.gz            # Step 5 单样本 gVCF（-ERC GVCF, 靶区±100bp）
├── cohort/                         # 批次内 CombineGVCFs→GenotypeGVCFs→norm→
│                                   #   SNP/INDEL 硬过滤→concat→PASS（小 panel 不用 VQSR）
├── matrix/                         # 基因型矩阵（全口径 + ./. 裁决后）与 PASS 位点详情表
├── per_sample_vcf/                 # 每样本 hardfiltered/PASS/裁决后 PASS VCF（工作副本）
├── qc/                             # fastqc_raw/fastqc_trim/fastp/flagstat/stats/
│   │                               #   hsmetrics/mosdepth/bcftools_stats/multiqc
├── run_report.md  run_summary.json
└── logs/                           # 入口日志(run_*) + 批次主日志(pipeline_*) + 每样本日志
```

**`results/` 根目录下只有 `<批次>_<执行日期>/` 目录**——不产生全局 run_summary、
全局 logs 等任何散文件（run_summary 逐批写入批次目录，最后完成的批次文件即全貌）。

**最终交付文件**是每样本 `*.PASS.adjudicated.vcf.gz`。Step 6 中**全部流程（含矩阵裁决与
每样本 VCF 重建）结束、QC 文件生成完全后**，MultiQC 生成最终汇总报告，随后自动导出到
**独立交付目录**（与分析工作区分离，`GWAS_DELIVERY_DIR` 可覆盖，默认
`$WORK/Output/<批次名>_<执行日期>/`；v2.3.0 前为 `delivery/`，MultiQC 报告同批交付）：

```
Output/
├── INDEX.md                            # 交付索引（累积：扫描全部历史交付目录 ∪ 本次运行）
└── 260422_20260914/
    ├── NA12878.PASS.adjudicated.vcf.gz (+.tbi)   # 交付主体
    ├── NA18544.PASS.adjudicated.vcf.gz (+.tbi)
    ├── GWAS-Panel-下机数据-QC_multiqc_report.html  # MultiQC 全流程 QC 汇总（v2.3.0 起）
    ├── md5sum.txt                      # 标准格式，接收方 md5sum -c 直接校验
    ├── MANIFEST.tsv                    # 样本/文件/大小/md5/记录数/来源（MultiQC 行 sample 列=multiqc）
    └── README.md                       # 交付口径说明（硬过滤→PASS→DP≥20 裁决重建）
```

导出幂等：源 VCF 未更新则不复制（mtime 比较），manifest/md5 每次重新生成。

### 3.1 results/ 归档（archive_results.py，v2.15.0/DEC-25）

批次结果目录随运行累积，`archive_results.py` 把**执行日期超过 30 天**（按目录名
`<批次>_<YYYYMMDD>` 后缀判定）的批次目录打包为 7z 并删除源目录（离线维护工具，
不属于流程运行时）：

```bash
python3 archive_results.py --dry-run   # 只列出将归档的目录，不动任何文件
python3 archive_results.py             # 归档 results/ 下超 30 天批次目录 → $WORK/archive/
python3 archive_results.py --days 60   # 自定义保留天数；--results/--out 可改输入输出目录
```

- 压缩口径（固定）：`7z a -t7z -mx=9 -mfb=192 -ms=on -md=256m -snl -mmt -sdel`
  ——**`-sdel` 由 7z 保证"压缩成功后才删源"**，失败自动保留；
- 幂等：输出目录已存在同名 `.7z` → 跳过（不重压不覆盖）；
- 非批次命名目录（无 `_YYYYMMDD` 后缀）与非法日期一律跳过；需要 `7z` 在 PATH。

## 4. CLI 参数

| 参数 | 说明 | 默认 |
| --- | --- | --- |
| `--input <目录>` | 多批次输入根目录，遍历其下全部批次子目录逐批执行 | `0_raw_data` |
| `--batch <批次名>` | 指定单个批次（与 `--input` 同时给出时**优先**） | - |
| `--step <0-6>` | 执行到第几步（6=质量汇总） | 6 |
| `--samples a,b` | 仅处理指定样本 | 全部 |
| `--exclude-samples NTC` | 从联合变异检测排除的对照（逗号分隔，空串关闭） | `NTC` |
| `--workers N` | 并行样本数（**最高优先**，覆盖规划） | 规划值 |
| `--threads N` | 每样本线程（覆盖规划） | 规划值 |
| `--max-memory NG` | 覆盖可用内存探测值 | 探测值 |
| `--resource-profile auto\|low\|high` | 资源档位（见下节） | `auto` |
| `--serial` | 强制单样本串行（= `--workers 1`） | - |
| `--dry-run` | 打印命令与资源计划不执行 | - |
| `--notify on\|off` | 钉钉通知开关 | `on` |
| `--notify-test` | 发送测试 markdown 后退出 | - |
| `--version` | 显示流程版本（与 DESIGN.md 同步）后退出 | - |
| `--out <目录>` | 覆盖批次结果目录 | `results/<批次>_<执行日期>` |

资源参数优先级：`--workers` > `--threads`/`--max-memory` > `--resource-profile` > 自动探测。

## 5. 资源自适应（resource.py，硬性要求）

探测：CPU = min(os.cpu_count, cgroup cpu.max 配额)（容器/cgroup 环境识别）；
内存 = min(MemAvailable, cgroup 限额)。预留 2 核 + max(2G, 可用 5%)。

推导：按可用核分档定每样本线程 T（≥64→24；32-63→12；16-31→8；<16→4），
workers(核)=⌊可用核/T⌋；再按内存收紧 workers（单样本峰值 = bwa-mem2 索引
17G(mmap 常驻) + T×sort缓冲 + GATK + 杂项 1G），取较小值；随后反推
SORT_MEM（128M-2G 钳位）与 GATK_MEM（1g-8g 钳位）。串行 cohort 步骤取可用内存
60%（-Xmx8g-32g 钳位）。**快速失败**：单样本峰值 > 可用内存（如 <20G 装不下人类
bwa-mem2 索引）启动即报错退出，禁止跑到一半 OOM。计划透明：启动打印完整推导表、
写入 run_summary.json、进入钉钉启动通知。

| 档位 | workers | T | sort -m | GATK -Xmx | 说明 |
| --- | --- | --- | --- | --- | --- |
| auto（自动探测） | ⌊可用核/T⌋ | 24 | 2G | 8g | 探测推导 |
| high（100 线程/900G） | 4 | 24 | 2G | 8g | 强制规划 |
| low（16 线程/20G） | 1 | 4 | 128M | 1g | 强制规划；峰值 19.5G≤20G |

auto 档超过 100 线程/900G 照常使用（不设人为上限），只扣系统预留。

**低配边界（要点）**：bwa-mem2 人类索引常驻 17G，mmap 之外运行时+缓存开销实测
≈2-4G；cgroup 硬上限场景实测下限约 24-26G，裸机 20G 属临界（依赖内核对 mmap
干净页的回收）。规划器对"单样本峰值 > 可用内存"的机器启动即快速失败，防跑一半
OOM。完整实测记录（复现命令与数据）见 `design_doc/RUN_HISTORY.md`（本仓库）。

## 6. 钉钉通知（企业机器人 · 分级告警 + 全步骤里程碑 + 交付文件推送）

### 6.1 渠道与消息结构（企业内部机器人，v2.15.0/DEC-24）

**企业内部应用机器人**（webhook 机器人不能发文件，已弃用）：markdown 群消息走
v1.0 `robot/groupMessages/send`，文件卡片走 oapi `media/upload` + `sampleFile`；
accessToken 进程内缓存（7200s，提前 5 分钟过期）。官方 markdown 子集仅：
标题/引用/文字效果/链接/图片/有序无序列表——**不支持表格**，换行需 `\n\n`，
发送前自动三重规范化（表格降级为列表、单换行提升、超长截断）。
`--notify off` 关闭；`DINGTALK_MILESTONES=0` 时只发异常级（P0/P1）里程碑、
OK 级静默。发送失败仅降级写日志、不影响退出码。
**凭证四键（CLIENT_ID/CLIENT_SECRET/ROBOT_CODE/CONVERSATION_ID）来自
`pipeline/.env`（见 §1.2），代码不存密钥**；未配置时通知静默跳过并提示配置方法。
`--notify-test` 实发一条 markdown + 一个小 zip 文件，验证双链路连通性。

### 6.2 通知时机

每批次：启动（含开跑前检查）→ Step 0-6 每步里程碑 → 全流程完成；多批次另有总览。
**交付文件推送（v2.15.0）：批次全部结束后**（多批次不逐批推，防通知淹没）逐成功
批把 `Output/<批次>_<日期>/` 打包临时 zip → 说明消息 + 文件卡片 → 清理临时包；
zip 超钉钉 20MB 上限时自动分卷压缩发送（每卷独立合法 zip ≤20MB、partNNofMM 命名，
全部下载后解压到同一目录即还原；单卷装不下的文件或卷数超上限时回落说明消息）。
模板：

```
[GWAS][P1] 20260720批次 · Step 2 比对完成
样本: 4/4 成功 | Lane 合并 16/16
指标: mapped 98.7% | proper pair 94.2%
异常: [P1] L20260615001 mapped 91.3%（阈值 95%）← 需确认
产物: bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22
日志: tail -f logs/sample_L20260615001.log
```

### 6.3 分级告警体系（v2.9.0 三级，DEC-21）

**P0=阻断级**（严重影响分析 → 直接中断批次，退出码 1，多批次运行其余批次照常隔离）；
**P1=严重**（执行失败按样本隔离 / 严重数据异常，报错不中断）；
**P2=质量提示**（阈值越界 / 口径存疑，只记录）。

| 时机 | 规则 | 级别 |
| --- | --- | --- |
| 开跑前 · 依赖文件 | sif / genome(.fai/.dict) / targets / known-sites 缺失或 0 字节 → 中断批次 | P0 |
| 开跑前 · 批次名/样本名 | 含非常规字符（A-Za-z0-9_.- 外，含中文；注入面+容器 locale 报错）→ 中断批次（DEC-19） | P0 |
| 开跑前 · 磁盘空间 | 剩余 < 200GB（每样本估算约需 75G）→ 中断批次 | P0 |
| 开跑前 · md5 | 输入校验失败（数据损坏）→ 中断批次（v2.9.0 前仅剔除该样本） | P0 |
| Step 0 · 合并失败 | cat 非 0（样本终止，其余照常） | P1 |
| Step 1-6 间 · 磁盘复查 | Step 结束时剩余 < 200GB → 终止批次（RUN-34，开跑前检查的运行中补位） | P0 |
| Step 1 · 上样量 | 样本 reads 绝对量 < 100 万 → 上样不足 | P1 |
| Step 6 · 空交付 | per-sample PASS 裁决 VCF 0 条记录（交付结果为空） | P1 |
| Step 1-6 · 执行失败 | 任一步 exit≠0 / 产物 0 字节（样本级隔离，不中断批次） | P1 |
| Step 2 · QC 口径 | mapped <90%（报错不中断；v2.9.0 前判样本失败） | P1 |
| Step 6 · 对照污染 | NTC 靶区深度 >10×（报错不中断；v2.9.0 前为 P0） | P1 |
| 全流程 · 失败样本 | 任一样本在某步失败（汇总报错） | P1 |
| Step 1 · 修剪 | fastp 保留率 <80% 或 Q30 <85%（95-80% 为提示行） | P2 |
| Step 2 · 比对 | mapped <95% 或 properly paired <85% | P2 |
| Step 3 · 去重 | 重复率 >30%（建库复杂度告急） | P2 |
| Step 5 · 对账·数量 | 矩阵行数 ≠ 靶区记录数（view -R 同口径对账，DEC-11） | P2 |
| Step 5 · 对账·新鲜度 | 关键 VCF mtime < 本次启动（断点续跑复用旧产物） | P2 |
| Step 6 · 捕获/覆盖/口径 | 捕获效率 PCT_SELECTED <85% / mean depth <50× / 20X <95% / Ti/Tv <2.0 / call rate <95% | P2 |
| Step 6 · 深度离散 | 批次内 mean depth 变异系数 CV >0.5（疑似混入异常样本） | P2 |
| Step 1 · NTC reads | NTC reads 占批次中位样本 >1%（污染维度之二，与深度互补） | P2 |
| Step 4 · 校准可信 | BQSR recal M 事件观测数 <10 万（known-sites 覆盖异常，校准不可信） | P2 |

> **样本点名规则（v2.20.0/DEC-30）**：逐指标**越界样本全点名**，每样本一行（名字序）——
> 此前每指标只点名最差一个样本（260918 自测三样本 on-target 全部越界只报一行，被误读为
> "只有一个样本坏"）。例外：fastp 保留率 80-95% 提示带为 OK 级聚合一行；批级指标
> （Ti/Tv、call rate、深度 CV）单值无点名问题。逐样本完整数值以 run_summary.json 为准。


## 7. 质检口径（自动判定，与笔记一致）

- fastp 保留率 <95% 告警；Adapter Content PASS<5%/WARN 5-20%/FAIL>20%，修剪后复检
- mapped>90%（<90 报 P1 不中断，v2.9.0 前曾判样本失败）；singletons 应很低；properly paired ≈97-98%（偏低伴跨染色体配对升高
  → 节段重复区正常现象，告警不报错）
- dup>30% 且 ELS 偏小 → 文库复杂度不足告警；去重前 flagstat 的 duplicates 行不采集
- **BQSR 前后 flagstat 逐行 diff 必须一致，不一致判 FAIL**（该样本退出变异检测）
- MEAN/MED_TARGET_COVERAGE ≥50×、PCT_TARGET_BASES_20X ≥90%；捕获效率 = PCT_SELECTED_BASES
  ≥85%（v2.19.0/DEC-29：(on+near bait)/比对碱基，与外送 pct_selected_bases 同口径；新 panel
  实测 89-91%）；on-target 降为信息指标不告警（1bp SNP panel 下 ≈0.6% 为几何产物）
- `*`（spanning deletion）不计入 SNP/INDEL；MIXED/多等位先 norm 摊平再分拣
- FILTER 列出现 `.` = 过滤漏跑，判错；Ti/Tv raw→PASS 应上升（看趋势不看绝对值）
- 矩阵 `./.`：mosdepth bqsr regions 深度 DP≥20（`GWAS_DP_MIN` 可配）改判 0/0，不足保留 `./.`；
  每样本 PASS VCF 重建只替换 GT，其余字段原样保留

## 8. 与笔记差异（冲突裁决第 ③ 类：笔记写死 → 工程化规划值）

以下均为**有意工程化**（提示词第四节硬性要求"禁止写死资源参数"），命令语义不变：

| # | 笔记原文 | 本程序 | 说明 |
| --- | --- | --- | --- |
| 1 | `bwa-mem2 mem -t 8`、`samtools sort -@ 8` | `-t/-@ = plan.threads` | 每样本线程由分档表推导（104 核机 24） |
| 2 | `samtools sort -m 4G` | `-m plan.sort_mem` | 128M-2G 钳位（低配 128M） |
| 3 | `fastp --thread 4`、`fastqc -t 4`、`mosdepth -t 4` | `plan.fastp_threads` 等 | 同上，fastp 上限 16 |
| 4 | GATK `--java-options -Xmx2g`（MarkDuplicates/BQSR/HC/Select/Filter） | `-Xmx plan.gatk_mem` | 1g-8g 钳位 |
| 5 | GATK `-Xmx4g`（CombineGVCFs/GenotypeGVCFs/CollectHsMetrics） | cohort 用 `-Xmx plan.cohort_mem`；HsMetrics 随每样本 `plan.gatk_mem` | cohort 串行放宽至 8g-32g |
| 6 | HC 3 个一批（`i%3` 控制） | workers 由规划器推导 | 104 核机 4 样本并行 |
| 7 | 笔记 5-4 每样本矩阵列号硬编码 `declare -A COL` | 按 calling 样本顺序动态计算列 | 消除对 12 样本硬编码 |
| 8 | 笔记散述"bcftools isec 求交集/差集" | norm 后 (CHROM,POS) 键集合 Python 交/差 | 等价实现：norm 后键唯一；GT 级比较仍需逐条查询 |
| 9 | 笔记 CollectHsMetrics 先用 markdup.bam（3-9）后用 BQSR bam（5-5） | 统一用 BQSR bam（5-5/Step 6 口径） | mosdepth 双跑（md+bqsr）保留对照 |
| 10 | 笔记输出至工作区根部（qc/ bam/ gvcf/ …） | 全部收进 `results/<批次>_<执行日期>/` | 批次隔离与安全边界要求 |
| 11 | 笔记散述 HsMetrics/捕获口径 | 捕获效率告警 = PCT_SELECTED ≥85%，on-target 不告警（信息指标）（v2.19.0/DEC-29） | 老结论"PCT_SELECTED(25-28%) 不误作捕获效率"系老 panel 观测，已废止——新 panel 下与外送 pct_selected_bases 交叉验证一致（<0.2pp 偏差） |
| 12 | 笔记 BedToIntervalList 无 `--UNIQUE`/`--DROP_MISSING_CONTIGS` | 固定 `--UNIQUE true --DROP_MISSING_CONTIGS true`（DEC-26） | 新 panel bed 含字典外 ALT contig（否则 PicardException 中断）；重叠/相邻探针区间合并为唯一区间，靶区碱基按唯一口径（实测 31310→18339bp） |

其余分析学内容（命令、参数语义、阈值：fastp length_required 36、bwa `-K 100000000 -Y`、
HC interval-padding 100、硬过滤 QD2/QUAL30/SOR3/FS60/MQ40/MQRankSum/ReadPosRankSum 与
INDEL 口径、DP_MIN 20、±100bp padding 等）与笔记**逐字一致**。

## 9. 模块结构

```
$WORK/
├── singularity/       # 容器镜像目录（pipeline 同级，固定位置，不走 .env）
└── pipeline/
    ├── .env               # 环境参数（webhook 等，见 §1.2；.env.example 为模板）
    ├── run_pipeline.py   # 主控：CLI、步骤编排、样本级并行、汇总、通知触发、交付导出
    ├── config.py         # 路径/线程/内存/阈值（.env + 环境变量 + 内置默认三源）
    ├── logger.py         # 主日志 + 每样本日志（线程安全、时间戳+耗时、stdio tee）
    ├── runner.py         # singularity 封装、实时输出、超时、返回码、幂等 SKIP
    ├── resource.py       # ★ 资源探测与规划（cgroup/WSL2 识别、workers×线程×内存推导）
    ├── scanner.py        # 三种输入布局扫描（识别器独立函数）+ md5 校验 + Lane 合并 + samples.tsv
    ├── dingtalk.py       # 钉钉企业机器人（markdown 群消息 + 文件卡片 + 交付 zip 推送）
    ├── alerts.py         # ★ 分级告警（P0/P1 阈值检查）+ 步骤里程碑消息构造
    ├── report.py         # 运行报告 markdown + run_summary.json
    ├── archive_results.py # results/ 超期批次目录归档 7z（-sdel，DEC-25）
    ├── check_design.py   # DESIGN.md TH 表 ↔ config.py 一致性校验（防漂移）
    ├── tests/            # 回归测试集（run_tests.sh 末尾自动执行 check_design）
    └── modules/          # 每软件一模块：fastqc fastp bwa_mem2 samtools gatk bcftools mosdepth multiqc
```
