# pipeline 单元测试集（迁移/改动后的快速回归）

**纯 Python 标准库，无容器、无网络、无真实数据，约 5 秒跑完 109 个用例。**
每个用例对应本项目开发/运行中的真实踩坑或硬性口径，改动代码后先跑本测试再跑真实数据。

## 运行

```bash
cd pipeline
./run_tests.sh                    # 或 python3 -m unittest discover -s tests -v
```

## 用例文件 ↔ 踩坑映射

| 文件 | 覆盖 | 对应踩坑/要求 |
| --- | --- | --- |
| `test_runner.py` | tool() 拼装、binds、cpath 路径映射、**rt 绝对路径解析与 PATH 兜底**、幂等 SKIP、dry-run 不执行、超时 kill、capture、**统计类命令 timeout=600 传递（RUN-31）** | binds 变量遮蔽曾致 tool() 崩溃；产物存在时不执行命令（`exit 7` 也返回 0）；无输出长命令超时曾不生效（本测试集修复）；PATH 受限环境启动曾致 `singularity: not found`（RUN-22，本测试集修复） |
| `test_resource.py` | sort -m 格式、gatk/cohort -Xmx 格式、low/high 档精确值、快速失败、计划透明表 | `samtools sort -m 2.0G` 被解析为 2 字节——必须整数+单位（`128M`）；low 档必须 workers=1/T=4/128M/1g；内存装不下索引必须启动即 SystemExit |
| `test_scanner.py` | Illumina 平铺/外送子目录/**外送平铺（RUN-36）三种布局**、S 号与 Lane、R1/R2 不匹配与 0 字节无效标记、**布局互斥锚（Illumina 命名不得误认外送平铺）**、**同名冲突先认者优先**、**md5 失败样本名推导覆盖外送平铺**、合并（真实与 dry-run）、samples.tsv 列 | dry-run 下合并曾被误判失败致批次失败；无效样本不得影响其余样本；外送交付平铺数据曾两种识别器都不认 → 整批"无有效样本"静默跳过（RUN-36 事故）；平铺文件 md5 失败样本名曾误取整个文件名致失败集与有效集无交集（RUN-33 病根） |
| `test_dingtalk.py` | 表格降级为列表、单换行→`\n\n`、超长截断、列表可渲染性、**发送成功落 INFO 日志** | 钉钉 markdown 官方子集**不含表格**、换行必须 `\n\n`（曾整条消息渲染成竖线串）；notify=on 曾零通知痕迹无法事后确认（RUN-29） |
| `test_envfile.py` | .env 解析语法（注释/export/引号/行内注释/URL 含=?&）、三源优先级（环境变量>.env>默认）、webhook 默认空、**源码防回潮锚（不得出现 access_token=）**、未配置时 send/notify 降级不抛异常 | 钉钉 webhook 曾硬编码在 config.py 默认值里（密钥随代码泄露，RUN-26 迁 pipeline/.env）；reload 类用例须在 finally 中恢复真实 config 防污染 |
| `test_alerts.py` | **P0/P1/P2 三级判定**（v2.9.0 DEC-21：质量阈值全量 P2、NTC 污染 P1、mapped<90 P1 报错不中断）、最差级别、里程碑模板字段、产物摘要 | 分级体系：P0=阻断中断 / P1=严重不中断 / P2=质量提示；NTC 污染曾为 P0（RUN-33 降级锚）；RUN-34 新增 reads 不足/NTC reads 占比/深度 CV/recal 观测数检查 |
| `test_parsers.py` | fastqc zip（小写状态）、fastp json、markdup/hsmetrics 按表头名、flagstat、samtools stats、bcftools stats、norm 统计、mosdepth summary | fastqc 状态是小写 pass/fail（曾按大写比较误报）；bcftools stats SN 行带文件 ID 列；norm 统计行是 7 字段动态表头；mosdepth 是 6 列且靶区深度取 `total_region` 行；flagstat `in total` 行无百分比（本测试集修复的潜伏 bug）；recal 观测数按 RecalTable1 表头名取列只加 M 行（RUN-34） |
| `test_variant_post.py` | 矩阵 `./.` 裁决（DP≥20/cells_total）、rebuild 只换 GT、GT 列显式映射、**norm_split 产物检查口径（norm 不查 .tbi）** | **GT 列序互换**曾致两样本基因型对调；矩阵列序必须按 calling 顺序显式映射；norm 命令 outputs 曾混入 .tbi（由后续 index 命令生成），首跑必误报"产物缺失"（RUN-29 修复） |
| `test_misc.py` | 日志 tee 自动落盘、Logger(None) 控制台模式、`results/<批次>_<日期>` 命名、**run_date 启动时固定（跨 0 点锚）**、**全流程 dry-run success+零落盘锚**、**外送平铺布局端到端锚（RUN-36）**、**样本名白名单 P0 阻断锚（非常规字符→中断分析，DEC-19 v2.8.0）**、**--version 锚**、**相对 --input 尊重 GWAS_RAW_DATA 覆盖锚（RUN-32 事故）**、**MultiQC 时序锚（裁决/重建之后、交付导出之前，v2.7.0 起锚定 step6 方法）**、交付导出（`Output/<批次>_<日期>` 命名、md5sum/MANIFEST/README/幂等 mtime、**MultiQC extra_files 布局与回退**）、**INDEX.md 累积锚（历史交付∪本次运行）** | 手动 shell 重定向曾被认为是必需（现 stdout 自动 tee 进运行日志）；单批次运行不得覆盖全局汇总；交付源未更新不得重拷；执行日期须取自启动时间戳一次性派生，否则跨午夜运行会写进两个日期目录；全流程 dry-run 曾因 gvcf.list 写入/merge 建目录/disk_usage 三处无守卫，仅同日实跑后可用、跨日直接 FileNotFoundError（RUN-27 修复）；MultiQC 时序为用户口径——全部流程结束 QC 齐全后才出最终报告再交付（RUN-27）；交付目录 v2.3.0 由 delivery 改名 Output 且 MultiQC 随交付；INDEX.md 曾整体重写只含本次运行批次，历史交付被挤出索引（RUN-29 修复） |

## 新服务器迁移后的建议验证顺序

1. `./run_tests.sh` —— 109 用例全绿（验证 Python 版本兼容与全部纯逻辑）
2. `cp .env.example .env`（pipeline/ 下）并填入钉钉地址 → `python3 run_pipeline.py --notify-test` —— 钉钉连通性（webhook/关键词，需能出网）
3. `python3 run_pipeline.py --dry-run --batch <小批次>` —— 容器/参考文件/路径可用性与资源计划推导
4. `python3 run_pipeline.py --resource-profile low --dry-run` —— 低配档口径（workers=1、sort 128M、GATK 1g）
5. 单批次冒烟（如 2 样本批次），再放全量

## 约束

- 测试不得访问网络（钉钉仅测纯函数 normalize 与未配置降级分支，不发真实请求）
- 测试不得依赖容器/真实数据目录（fixture 全部在 tempfile 中构造）
- 若新增模块/解析器，请同步补用例并在本表登记对应踩坑
