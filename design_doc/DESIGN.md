# GWAS 靶向测序 FASTQ→VCF 流程 · 设计文档（人机共写）

```yaml
# ---- design-meta（机器可解析锚点，勿手改格式；版本规则见 §0）----
doc: GWAS-pipeline-design
version: 2.6.0
updated: 2026-09-16
owner_human: gewenlong
owner_machine: ZCode(GLM)
source_of_truth:
  goals_and_scope: pipeline/design_doc/DESIGN.md          # 本文档（人写区为唯一目标来源）
  thresholds: pipeline/config.py                 # 阈值唯一实现源；本文档 TH 表为镜像
  commands_reference: pipeline/design_doc/notes_code_reference.md   # 思源笔记《下机数据流程》快照
  run_evidence: pipeline/design_doc/RUN_HISTORY.md        # 机器回写运行台账
consistency_check: pipeline/check_design.py      # TH 表 ↔ config.py 防漂移校验
```

---

## 0. 文档使用协议（人机共写规范）

**目的**：以本文档为中心演进项目，防止反复修改导致目标漂移。

**分区标记**：
- `<!-- HUMAN -->` 区：目标、范围、口径、验收标准——**人**编写/修改，机器只读；
- `<!-- MACHINE -->` 区：状态、证据、运行结果——**机器**回写，人可批注（批注行以 `>>` 开头）。

**编号锚点**（代码注释、测试、日志、通知引用这些编号，形成可追溯链）：
- `REQ-xx` 需求；`DEC-xx` 设计决策；`TH-xx` 阈值；`RUN-xx` 运行轮次（台账在 RUN_HISTORY.md）。

**修改流程**：
1. 人修改 HUMAN 区（目标/阈值/验收）→ 在 §8 CHANGELOG 记一条（版本号+1：口径微调 `minor+1`，目标/范围变更 `major+1`）；
2. 机器解析差异 → 对照实现（代码/测试同步引用编号）；
3. 机器回写 MACHINE 区状态与 RUN_HISTORY 台账；
4. 运行 `python3 pipeline/check_design.py`——TH 表与 config.py 不一致即失败（防两处漂移）；`./pipeline/run_tests.sh` 末尾也会自动执行该校验。

**冲突裁决顺序**（继承自提示词）：分析学内容以思源笔记原文为准（可经思源 MCP 回查）→ 工程化要求以本文档为准 → 笔记写死而本文档改为规划值的参数属有意工程化（记录于 pipeline/README.md §8"与笔记差异"）。

---

## 1. 项目目标与范围 <!-- HUMAN -->

把「人类 GWAS 靶向测序：FASTQ → 质控修剪 → 比对 → 去重 → BQSR 校准 → 变异检测
（gVCF 联合分型 → 硬过滤 → VCF）→ 测序质量汇总」
实现为一套 **Python3 仅标准库、模块化、可断点续跑、资源自适应（16线程/20G ～ 100线程/900G+）、
带钉钉分级通知** 的流程控制程序；最终交付物为每样本 `*.PASS.adjudicated.vcf.gz`
（独立 Output 交付目录 + MultiQC 汇总报告 + md5 校验）。

**范围边界**：批次是最小分析单元且完全独立（严禁跨批次合并）；`0_raw_data/`、`back/`
只读；交付产物写入独立 `Output/` 目录（v2.3.0 前为 `delivery/`）；NTC 等对照纳入 QC、
排除出联合变异检测；不使用 VQSR（样本量约 30 用硬过滤）。

## 2. 需求清单与状态 <!-- REQ 目标列=HUMAN；状态列=MACHINE -->

| 编号 | 需求（目标） | 状态 | 证据 |
| --- | --- | --- | --- |
| REQ-01 | Step 0-6 全流程（笔记 0-5 口径逐字一致） | ✅ 完成 | RUN-05/07/08/19 |
| REQ-02 | 模块化架构（run_pipeline/config/logger/runner/resource/scanner/dingtalk/report/alerts + modules/8 软件模块） | ✅ 完成 | pipeline/ 目录 |
| REQ-03 | 资源自适应：探测(cgroup/WSL2/affinity)→分档推导→内存收紧→钳位；禁写死参数；快速失败；计划透明 | ✅ 完成 | RUN-10~14；DEC-03 |
| REQ-04 | 幂等断点续跑：产物存在即 SKIP，只补缺失（含 .tbi/.metrics/.table） | ✅ 完成 | RUN-09（31 SKIP 结果一致） |
| REQ-05 | 批次独立 + 多批次失败隔离；任一批失败退出码非零 | ✅ 完成 | RUN-07/08 |
| REQ-06 | 输入校验：两种布局、R1/R2 不匹配/单端/0字节/命名无效标记跳过；md5sum.txt 并行校验 | ✅ 完成 | tests/test_scanner.py |
| REQ-07 | 安全边界：0_raw_data/back 只读；产物仅写 results/（与 Output/ 交付目录；v2.3.0 前为 delivery/，Output 只读旧口径废止） | ✅ 完成 | 全程无违规写入 |
| REQ-08 | 钉钉通知：markdown 官方子集、启动/全步骤里程碑/完成/失败、P0/P1 分级（§5 阈值表）、失败降级不中断 | ✅ 完成 | RUN-16；DEC-08 |
| REQ-10 | CLI 参数全集（--input/--batch/--step/--samples/--exclude-samples/--workers/--threads/--max-memory/--resource-profile/--serial/--dry-run/--notify/--notify-test/--out） | ✅ 完成 | run_pipeline.py argparse |
| REQ-11 | 日志自动落盘（tee 接管 stdio，无需 shell 重定向）+ run_summary.json 结构化 | ✅ 完成 | RUN-16；DEC-09 |
| REQ-12 | 提示词第七节验收测试 1-7 全部执行并留存证据 | ✅ 完成 | RUN-01~15；tests/ |
| REQ-13 | 最终交付独立目录 Output/<批次>_<日期>/（v2.3.0 前为 delivery/）：VCF+tbi、**MultiQC 汇总报告（v2.3.0 起）**、md5sum.txt（可 -c 校验）、MANIFEST.tsv、README，幂等导出 | ✅ 完成 | RUN-18/27 |
| REQ-14 | 多次运行隔离：results/<批次>_<执行日期>/，同日同目录续跑、跨日新目录；**执行日期在进程启动时取一次（跨 0 点不切换目录）** | ✅ 完成 | RUN-19；DEC-01；tests/test_misc.py 跨午夜锚 |
| REQ-15 | 小测试用例集（纯标准库秒级回归，覆盖历史踩坑），供改动/迁移快速验证 | ✅ 完成 | RUN-20；tests/ |
| REQ-16 | Step6 on-target 告警阈值 8%（小 panel 正常 8-10%，TH-14；原 60% 为口径误用） | ✅ 完成 | config.ON_TARGET_P1；笔记 5-5 |

> 新增需求：人追加 REQ-16… 行并填目标列，状态列留空由机器回写。

## 3. 系统架构 <!-- MACHINE 维护，人可提议变更（记 DEC） -->

```
$WORK/
├── singularity/       容器镜像目录（pipeline 同级，固定位置，不走 .env）
└── pipeline/
    ├── .env               环境参数文件（DEC-18；钉钉 webhook 等密钥只存于此，.env.example 为模板）
    ├── run_pipeline.py   主控：CLI/编排/样本级并行/汇总/通知触发/交付导出
    ├── config.py         ★ 阈值与路径唯一实现源（TH 表镜像于此；启动时加载 .env）
    ├── logger.py         双通道日志 + capture_stdio tee
    ├── runner.py         容器命令封装：幂等 SKIP/实时输出/超时 Timer kill/返回码/产物校验
    ├── resource.py       ★ 资源探测与规划（cpu=min(os/affinity/cgroup)；mem=min(MemAvailable/cgroup)）
    ├── scanner.py        两种布局扫描 + md5 + Lane 合并 + samples.tsv
    ├── dingtalk.py       markdown 通知 + 三重规范化（表格降级/\n→\n\n/截断）+ 未配置降级
    ├── alerts.py         P0/P1 阈值检查 + 里程碑消息构造（模板见 §5 后）
    ├── report.py         run_report/run_summary
    ├── check_design.py   DESIGN.md TH 表 ↔ config.py 一致性校验（防漂移）
    ├── tests/            回归用例（踩坑映射见 tests/README.md）
    └── modules/          fastqc fastp bwa_mem2 samtools gatk bcftools mosdepth multiqc
```

数据流：`0_raw_data/<批次>/ → results/<批次>_<日期>/{fastq_merged→fastq_clean→bam→gvcf→cohort→matrix→per_sample_vcf→qc} → Output/<批次>_<日期>/`。

## 4. 关键设计决策（DEC）<!-- MACHINE 回写；人可新增决策行（含动机） -->

| 编号 | 决策 | 动机/证据 |
| --- | --- | --- |
| DEC-01 | 批次完全独立；输出 `results/<批次>_<执行日期>/`；**执行日期取自进程启动时间戳 `run_date = ts.split("_")[0]` 一次性派生并传参到各批次——运行跨 0 点不切换目录**（同日重跑同目录幂等续跑；跨日中断后重跑进入新日期目录全新计算，属设计行为） | 多次运行隔离 + 同日幂等续跑（RUN-19）；跨午夜语义由 tests/test_misc.py::test_run_date_fixed_at_startup 锁定 |
| DEC-02 | 交付独立目录 + 标准 md5sum.txt/MANIFEST；v2.3.0 起目录名 delivery→**Output**、MultiQC 最终报告在**全部流程（含矩阵裁决与每样本 VCF 重建）结束、QC 文件生成完全后**才生成，随后连报告一起交付导出（时序锚 tests/test_misc.py） | 与输入侧 md5 约定对称，接收方可 `md5sum -c`（RUN-18）；改名、MultiQC 纳入与时序均为用户要求（RUN-27） |
| DEC-03 | 资源模型：索引 17G + 运行时；sort=剩余预算×0.25 钳位 128M-2G（整数 MB）；GATK 1-8g；cohort=可用 60% 钳位 8-32g；峰值>可用→SystemExit | `sort -m 2.0G`→2 字节坑（RUN-02）；cgroup 实测下限 24-26G（RUN-11~14） |
| DEC-04 | 命令/参数语义以笔记为准；资源参数规划化；差异记录于 pipeline/README §8 | 冲突裁决顺序（提示词第三节） |
| DEC-05 | 样本列序全链 sorted（calling/gvcf.list/裁决矩阵/重建） | GT 列互换 bug：曾致两样本基因型对调（RUN-07→08） |
| DEC-08 | 钉钉：标题 `[GWAS][P级别]`（关键词实测只校验正文且大小写敏感，正文自动兜底小写 gwas）；表格自动降级列表；DINGTALK_MILESTONES=0 只发异常级 | 官方文档实测（RUN-16 前的探测实验） |
| DEC-09 | 日志 tee：capture_stdio **缓冲模式**接管 stdout/stderr，批次发现后 tee_set_log_path 落盘到批次目录（单批次=该批次 logs/；多批次=首个批次 logs/）；dry-run/启动即退缓冲丢弃零落盘 | 用户要求"日志自动输出"（RUN-16）；v2.1.0 起 results/ 根禁止散文件（RUN-24） |
| DEC-10 | run_summary.json **逐批写入该批次目录**（快照含截至该批的 run+资源计划+各批次状态，最后完成的批次文件即全貌）；results/ 根不再写全局文件；dry-run 不写（防覆盖实跑） | 曾写全局 results/run_summary.json 且 --input 单批次也落全局（RUN-24 用户指出多余）；v2.1.0 收紧 |
| DEC-11 | 对账同口径：矩阵(query -R) 与核对(view -R) 一致；新鲜度=关键 VCF mtime<本次启动 | -R/-T 对跨界 indel 取舍不同（RUN-19 发现 7399 vs 7390） |
| DEC-12 | 测试集 = 踩坑回归集；迁移验证顺序：tests→notify-test→dry-run→low dry-run→冒烟→全量 | 已抓出 flagstat total 行、runner 静默超时 2 个潜伏 bug（RUN-20） |
| DEC-13 | NTC 纳入 QC（含 mosdepth/HsMetrics）、排除联合检测；对照样本阈值只记录不判级 | 提示词；NTC 近零覆盖属预期 |
| DEC-14 | 严重级区分：流程失败/断言不一致=ERROR；数据质量违例（覆盖度/保留率）=WARN | "无遗留 ERROR"验收口径需要语义干净的 ERROR |
| DEC-15 | on-target P1 阈值 60%→8% | 60% 误用了宽口径，小 panel 正常 on-target≈8-10%（笔记 5-5；用户 v2.0.0 决定） |
| DEC-16 | 容器运行时在 Runner 初始化时自解析为**绝对路径**（shutil.which → 标准目录探测），且子进程 PATH 兜底补齐标准系统目录——不依赖调用方 PATH 完整性 | 实测：从 PATH 受限环境（PATH=/usr/bin:/bin，IDE 面板/精简 env 类）启动曾致 `singularity: not found` exit=127（RUN-22）；singularity 实装于 /usr/local/bin/singularity |
| DEC-17 | **流程只在宿主 Python 运行**：main() 启动卫兵检测 /.singularity.d 或 SINGULARITY_* 环境变量，容器内运行立即退出并指引用 /usr/bin/python3；本机 `python` 是容器别名（mamba 解释器），嵌套容器既看不到宿主 singularity、时区还是 UTC | RUN-22/23：alias python 嵌套运行全部 127；流程纯标准库无需 mamba/conda |
| DEC-18 | **环境参数去硬编码，统一 .env**：钉钉 webhook（含 access_token）等环境参数只存 `pipeline/.env`（与代码同目录；模板 .env.example；`GWAS_ENV_FILE` 可改址），config.py 启动时解析并入（setdefault），三源优先级=进程环境变量 > .env > 内置默认；webhook 未配置→通知静默跳过（启动 WARN + 首次发送返回配置指引）；`singularity/` 镜像目录仍为 pipeline 同级目录约定，不走 .env | 用户要求（RUN-26）：密钥硬编码在代码里会随代码分发泄露且换环境必改源码；.env 权限 600 |

## 5. 统一口径与阈值总表（TH = config.py 镜像）<!-- HUMAN 可改值；MACHINE 同步 config 后过校验 -->

> **单一事实来源是 `pipeline/config.py`**；下表改动必须同步 config.py（或让机器代改），
> `python3 pipeline/check_design.py` 校验两处一致，`run_tests.sh` 末尾自动执行。

| 编号 | config 键 | 值 | 语义 |
| --- | --- | --- | --- |
| TH-01 | FASTP_LENGTH_REQUIRED | 36 | fastp 最短读长（笔记） |
| TH-02 | FASTP_RETENTION_WARN | 95.0 | 保留率<95% 提示 |
| TH-03 | FASTP_RETENTION_P1 | 80.0 | 保留率<80% → P1 |
| TH-04 | FASTP_Q30_P1 | 85.0 | Q30<85% → P1 |
| TH-05 | MAPPED_MIN_PCT | 90.0 | mapped<90% ERROR（QC 口径） |
| TH-06 | MAPPED_NOTIFY_P1 | 95.0 | mapped<95% → P1（通知从严） |
| TH-07 | PROPER_PAIR_P1 | 85.0 | properly paired<85% → P1 |
| TH-08 | DUP_WARN_PCT | 30.0 | dup>30% 且 ELS 小 → 复杂度告警 |
| TH-09 | DUP_P1 | 30.0 | dup>30% → P1 |
| TH-10 | MEAN_COV_MIN | 50.0 | MEAN 靶深度<50× WARN |
| TH-11 | MEAN_DEPTH_P1 | 50.0 | mean depth<50× → P1 |
| TH-12 | PCT_20X_MIN | 90.0 | 20X<90% WARN |
| TH-13 | PCT_20X_P1 | 95.0 | ≥20x 靶比例<95% → P1 |
| TH-14 | ON_TARGET_P1 | 8.0 | on-target<8% → P1（小 panel 正常 8-10%，笔记 5-5；v2.0.0 由 60 修正） |
| TH-15 | DP_MIN | 20 | `./.` 裁决深度阈值 |
| TH-16 | HC_INTERVAL_PADDING | 100 | HC 靶区外扩 bp（笔记） |
| TH-21 | NTC_DEPTH_P0 | 10.0 | NTC 靶深>10× → P0 污染 |
| TH-22 | TITV_P1 | 2.0 | Ti/Tv<2.0 → P1 |
| TH-23 | CALL_RATE_P1 | 95.0 | call rate<95% → P1 |
| TH-24 | DISK_MIN_FREE_GB | 200 | 开跑前剩余磁盘 P0 线 |
| TH-25 | DISK_PER_SAMPLE_GB | 75 | 每样本磁盘估算 |
| TH-26 | BWA_INDEX_MEM_GB | 17.0 | bwa-mem2 索引 mmap 常驻 |
| TH-27 | SORT_MEM_MIN | 128 | sort 每线程下限 MB |
| TH-28 | SORT_MEM_MAX | 2048 | sort 每线程上限 MB |
| TH-29 | GATK_MEM_MIN_GB | 1 | GATK -Xmx 下限 g |
| TH-30 | GATK_MEM_MAX_GB | 8 | GATK -Xmx 上限 g |
| TH-31 | COHORT_MEM_MIN_GB | 8 | cohort -Xmx 下限 g |
| TH-32 | COHORT_MEM_MAX_GB | 32 | cohort -Xmx 上限 g |

> TH-17~20 编号已废弃不复用；TH-21 起编号保持不变以维持引用稳定。

钉钉里程碑模板（alerts.step_milestone）：

```
[GWAS][P1] 20260720批次 · Step 2 比对完成
样本: 4/4 成功 | Lane 合并 16/16
指标: mapped 98.7% | proper pair 94.2%
异常: [P1] L20260615001 mapped 91.3%（阈值 95%）← 需确认
产物: bam/*/*.sort.bam ×4  mtime 2026-09-14 15:22
日志: tail -f logs/sample_L20260615001.log
```

## 6. 输出与目录规约 <!-- MACHINE -->

- **`results/` 根下只允许 `<批次>_<执行日期>/` 目录，无任何散文件/子目录**（RUN-24 起）
- `results/<批次>_<执行日期>/`：samples.tsv、fastq_merged/fastq_clean/bam/gvcf/cohort/matrix/per_sample_vcf/qc、run_report.md、run_summary.json、logs/（执行日期启动时固定，跨 0 点不切换，见 DEC-01）
- `Output/<批次>_<执行日期>/`：`*.PASS.adjudicated.vcf.gz(+.tbi)`、`*multiqc_report.html`（v2.3.0 起）、md5sum.txt、MANIFEST.tsv、README.md；`Output/INDEX.md` 跨批次**累积**索引（扫描全部历史交付目录 ∪ 本次运行合并生成，RUN-29——曾只写本次运行批次，跨日运行会把历史交付挤出索引）
- 日志层级：`<批次>/logs/run_<ts>.log`（入口 tee，多批次写首个批次目录）→ `logs/pipeline_<ts>.log` → `sample_<样本>.log`；dry-run 全部不落盘
- run_summary.json：逐批写批次目录快照（DEC-10）；`--out` 显式覆盖时写该目录
## 7. 已知边界与风险 <!-- MACHINE -->

1. **cgroup 硬上限**：bwa-mem2 17G 索引 + 运行时/页缓存，硬上限 <24GiB 会被 OOM-kill
   （20GiB/24GiB-旧sort 实测失败；修正 sort 分配后 24G 通过，32G 实测峰值 25.59GB）。
   裸机 20G 属临界（依赖内核回收 mmap 干净页）。规划器对装不下的机器启动即快速失败。
2. **CPUQuota 在 WSL2 用户 slice 不生效**（cpu.max 缺失）：限核用 taskset（规划器已识别
   sched_getaffinity），迁移到 systemd 托管机时 cpu.max 可被识别。
3. 思源笔记为活文档：手册只是快照，分析学口径冲突时以笔记原文（MCP 回查）为准。
4. 跨日中断重跑进入新日期目录从头计算（同日才续跑）：长时间批次建议白天启动或
   次日同日补跑未完批次（DEC-01 设计行为）。

## 8. 变更日志（CHANGELOG）<!-- 人机共写：每方改动各记一行 -->

| 版本 | 日期 | 角色 | 变更 |
| --- | --- | --- | --- |
| 1.0.0 | 2026-09-14 | 人 | 依提示词确立目标与验收（REQ-01~12） |
| 1.1.0 | 2026-09-14 | 人 | 新增 REQ-13 交付目录 / REQ-14 多次运行目录隔离 / REQ-15 测试集 |
| 1.1.0 | 2026-09-14 | 机 | 完成实现与验证（RUN-01~20），建立本文档与 RUN_HISTORY、check_design 防漂移校验 |
| 2.0.0 | 2026-09-14 | 人 | 确认 REQ-14 执行日期启动时固定（跨 0 点不切换）；on-target 阈值 60%→8%（REQ-16/TH-14，DEC-15） |
| 2.0.0 | 2026-09-14 | 机 | 测试 75→69（含新增跨午夜锚）；tests+check_design 全绿；RUN-21 |
| 2.0.1 | 2026-09-14 | 机 | 健壮性修复（DEC-16）：Runner 自解析容器运行时绝对路径 + 子进程 PATH 兜底，PATH 受限环境可正常运行（复现验证 Step 1 成功 2/2）；测试 69→71 |
| 2.0.2 | 2026-09-14 | 机 | 根因修正（DEC-17）：用户提供关键事实——本机 `python` 为 singularity 容器别名，真实根因是**嵌套容器**（上轮 TZ=UTC 提示的来源）；新增启动卫兵拒绝容器内运行并指引宿主 python3；测试 71→73 |
| 2.1.0 | 2026-09-14 | 人 | 目录规约收紧（RUN-24）：results/ 根只允许批次_日期目录；运行日志归入批次目录；质疑全局 run_summary 的必要性 |
| 2.1.0 | 2026-09-14 | 机 | tee 改缓冲模式延迟落盘（单批次=本批次 logs/、多批次=首个批次、dry-run 零落盘）；run_summary 逐批写批次目录、不再写全局、dry-run 不写；历史散文件迁移清理；测试 73→76；实跑 --step 0 验证落盘位置 |
| 2.1.1 | 2026-09-15 | 机 | Q30 解析修复（RUN-25）：fastp json 键名按真实产物校准（整体 q30_rate 小数→百分数；R1/R2 由 read{1,2}_after_filtering 的 q30_bases/total_bases 自算），采集端统一百分数口径；test fixture 重建为真实结构；真实产物验证 85.76/86.82% |
| 2.2.0 | 2026-09-16 | 人 | 环境参数（钉钉通知地址等）不得硬编码进程序，改为 .env 文件实现；singularity 软件位置不变（pipeline 同级 singularity/ 目录） |
| 2.2.0 | 2026-09-16 | 机 | 实现 DEC-18：config.py 增 .env 解析加载（三源优先级），移除硬编码 webhook 默认值，密钥迁入 pipeline/.env（600，与代码同目录）+ .env.example 模板；dingtalk 未配置降级（静默跳过+配置指引）；启动日志报告 .env 加载状态；测试 75→86（test_envfile.py：解析语法/优先级/防回潮锚/降级）；`--notify-test` 实发验证成功（RUN-26） |
| 2.3.0 | 2026-09-16 | 人 | 交付目录 `delivery/` 改名 **`Output/`**；MultiQC 汇总报告随交付一并放入交付目录 |
| 2.3.0 | 2026-09-16 | 机 | 实现 DEC-02 更新：DELIVERY_DIR 默认 `$WORK/Output`（GWAS_DELIVERY_DIR 覆盖不变）；MultiQC 移至 Step 6 末尾——**全部流程（含裁决与每样本 VCF 重建）结束后生成最终报告，随后才交付导出**（用户澄清时序），export_delivery 增 extra_files（MultiQC 平铺进交付目录，幂等 mtime、入 md5sum.txt/MANIFEST/README）；旧 `Output/` 只读口径废止（改 0_raw_data/back 只读）；历史 delivery/ 迁移为 Output/（INDEX.md 路径同步）；顺带修复 3 处 dry-run 潜伏守卫缺失（gvcf.list 写入/merge 建目录/disk_usage，全流程 dry-run 曾仅同日实跑后可用）；测试 86→89（extra_files 布局/幂等/回退 + 全流程 dry-run 零落盘锚 + MultiQC 时序锚）；真实产物重导出验证 md5sum -c 全过（RUN-27） |
| 2.4.0 | 2026-09-16 | 人 | 文档清理：删除全部与 Illumina 比对相关的表述（比对不属于本流程），并清理已废弃内容（删除项经确认） |
| 2.4.0 | 2026-09-16 | 机 | DESIGN.md：删 §7 比对数据表（节号 8/9/10→7/8/9 重排）、REQ-09/DEC-06/07 行、DEC-15 收敛为 on-target 阈值口径、TH-17~20 注简化、CHANGELOG 2.0.0 行去比对字样、证据索引修正（删已废弃全局 run_summary 行与比对报告行、测试数→89、results/logs 路径）；RUN_HISTORY：总览金标准要点删除、各行一致率/比对报告字样清除（保留 GT 列互换等修复链）、指标表删 SNP一致率列；pipeline/README、tests/README、config.py 注释同步（-T→-R 修正）；PROMPT 原件经确认保留；比对产物文件实际已不存在（清 compare.pyc 残留）；tests+check_design 全绿（RUN-28） |
| 2.5.0 | 2026-09-16 | 人 | 依据 0_raw_data_test 实跑日志分析（RUN-29 前置）确认执行 4 项修复：norm 误报 WARN、对账日志文案 -T→-R、钉钉发送成功落日志、交付 INDEX 累积 |
| 2.5.0 | 2026-09-16 | 机 | ① norm_split 的 norm 命令 outputs 去 .tbi（该索引由后续 index 命令生成，混入检查必误报"产物缺失"）；② 对账 RESULT 文案 view -T→-R（命令本身已是 -R，DEC-11，文案过时误导排查）；③ send_markdown 成功时 logger.info 落一行"钉钉已发送: 标题"（此前 notify=on 但日志零通知痕迹）；④ INDEX.md 改累积合并（write_delivery_index：扫描 Output/ 全部历史交付目录 ∪ 本次运行；曾整体重写致历史交付 260422_20260914 从索引消失）并即时重生成找回；测试 89→92（norm outputs 口径/发送成功日志/索引累积锚）；tests+check_design 全绿；dry-run 零落盘复验（RUN-29） |
| 2.6.0 | 2026-09-16 | 人 | design_doc 移入 pipeline/ 随仓库管理（PROMPT 与内部版 README_pipeline 不入库）；仓库增加 GitHub Actions CI |
| 2.6.0 | 2026-09-16 | 机 | design_doc/ 迁至 pipeline/ 内（DESIGN/notes_code_reference/RUN_HISTORY 入库；PROMPT_GWAS_pipeline.md 含 webhook token、README_pipeline.md 自我声明不对外——两者进 .gitignore）；check_design 路径改为 pipeline 内相对（克隆仓库即可全绿跑 run_tests.sh）；新增 .github/workflows/ci.yml（py3.10/3.12 矩阵跑 run_tests.sh）；测试 CI 兼容化：rt 解析用例自包含（临时 PATH 植入假运行时，无 singularity 环境可验）、两个端到端 dry-run 用例加 --resource-profile low（CI runner 内存 ~16G < auto 档单样本峰值 19.8G 会快速失败）；92 tests+check_design 全绿（RUN-30） |

## 9. 证据索引 <!-- MACHINE -->

| 内容 | 路径 |
| --- | --- |
| 运行台账（每轮） | pipeline/design_doc/RUN_HISTORY.md |
| 最终验收运行（0_raw_data_test） | results/260422_20260914/、results/260529_20260914/（运行日志在各批次 logs/） |
| 交付 | Output/<批次>_20260914/ + INDEX.md（v2.3.0 前为 delivery/，已迁移并补 MultiQC） |
| 测试集 | pipeline/tests/（89 用例）+ pipeline/run_tests.sh |
| 使用说明/与笔记差异 | pipeline/README.md |
| 环境参数文件 | pipeline/.env（密钥，600）+ pipeline/.env.example（模板） |
