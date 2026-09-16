# 思源笔记《下机数据流程》代码备查手册

> 导出时间：2026-09-14（由思源笔记 MCP 逐篇读取后自动提取）
> 笔记位置：/生信相关/项目笔记/GWAS-全基因组关联分析/下机数据流程/
> 用途：构建 pipeline 各软件模块时核对命令与参数。**代码块按笔记原文逐字导出（含笔记中的注释），
> 其中部分块是笔记对命令的讲解/伪代码而非可直接执行的脚本，使用时以 Step 语义为准。**
> 若与提示词冲突：以笔记原文（可经思源 MCP 回查）为准。


---

## 笔记《0-准备》  （共 1 个代码块，其中 1 个含流程命令）


### 0-准备 · 代码块 0

```bash
#!/bin/bash
# download_reference.sh
export http_proxy=127.0.0.1:7890
export https_proxy=127.0.0.1:7890

set -exuo pipefail
# -e 命令出错立即退出
# -x 打印执行过程
# -u 使用未定义变量立即报错
# -o pipefail 管道中任一环节失败即整体失败


# 以脚本所在目录为基准，保证从任何位置调用都一致
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 支持外部覆盖：REF_DIR=/data/ref bash download_reference.sh
REF_DIR="${REF_DIR:-$SCRIPT_DIR/reference}"
RAW_NAME="raw"
RAW_DIR="$REF_DIR/$RAW_NAME"

# 容器镜像下载
SINGULARITY_DIR="${SINGULARITY_DIR:-$SCRIPT_DIR/singularity}"
mkdir -p "${SINGULARITY_DIR}"

mkdir -p "${RAW_DIR}"
cd "${REF_DIR}"

# ===== GATK Resource Bundle (Google Cloud) =====
# 需要 gsutil 或通过 HTTPS 下载

# 参考基因组 FASTA
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
ln -sfn "${RAW_NAME}/Homo_sapiens_assembly38.fasta" "${REF_DIR}/hg38.fa"

# FASTA 索引
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
ln -sfn "${RAW_NAME}/Homo_sapiens_assembly38.fasta.fai" "${REF_DIR}/hg38.fa.fai"

# 序列字典
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict
ln -sfn "${RAW_NAME}/Homo_sapiens_assembly38.dict" "${REF_DIR}/hg38.dict"

# dbSNP 138
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
bgzip -f -@ "$(nproc)" "${RAW_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"
ln -sfn "${RAW_NAME}/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
tabix -p vcf "${RAW_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
ln -sfn "${RAW_NAME}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" "${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"

# Mills and 1000G gold standard indels
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
ln -sfn "${RAW_NAME}"/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* .

# af-only-gnomad (可选，用于 CNV/Somatic 流程，此处备用)
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -c -P "${RAW_DIR}" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

cd "${SINGULARITY_DIR}"

# 安装并行parallel
# 1) 下载（已存在则跳过，便于脚本重跑）
if [ ! -f "${PARALLEL_TGZ}" ]; then
    wget -q "https://mirrors.tuna.tsinghua.edu.cn/gnu/parallel/${PARALLEL_TGZ}" \
        || { echo "[ERROR] 下载失败: ${PARALLEL_TGZ}"; exit 1; }
fi

# 2) 解压（目录已存在则跳过，避免重复解压/中断重跑时报错）
if [ ! -d "${PARALLEL_SRC}" ]; then
    tar -xjf "${PARALLEL_TGZ}" \
        || { echo "[ERROR] 解压失败，宿主缺 bzip2？可先: sudo apt-get install -y bzip2"; exit 1; }
fi

# 3) 相对软链接：${SINGULARITY_DIR}/parallel -> parallel-20260822/src/parallel
chmod +x "${PARALLEL_SRC}/src/parallel"
ln -sfn "${PARALLEL_SRC}/src/parallel" parallel

# 4) 验证（用 ./ 前缀避免 PATH 里同名命令干扰）
./parallel --version

declare -A IMAGES=(
    ["fastqc_0.12.1"]="quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    ["fastp_1.3.6"]="quay.io/biocontainers/fastp:1.3.6--h43da1c4_0"
    ["bwa-mem2_2.3"]="quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0"
    ["samtools_1.24"]="quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    ["gatk_4.6.2.0"]="quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
    ["bcftools_1.24"]="quay.io/biocontainers/bcftools:1.24--h487d631_1"
    ["mosdepth_0.3.14"]="quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1"
    ["multiqc_1.35"]="quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1"
)
for name in "${!IMAGES[@]}"; do
    SIF="${SINGULARITY_DIR}/${name}.sif"
    echo "=========================================="
    echo "Pulling: ${name}"
    echo "Source:  ${IMAGES[$name]}"
    echo "=========================================="
    
    # 幂等：已存在则跳过
    if [ -f "${SIF}" ]; then
        echo "[SKIP] ${name}.sif 已存在"
        echo ""
        continue
    fi
    sudo docker pull "${IMAGES[$name]}"
    sudo singularity build "${SINGULARITY_DIR}/${name}.sif" docker-daemon://"${IMAGES[$name]}"
    #sudo singularity pull --name "${SINGULARITY_DIR}/${name}.sif" "${IMAGES[$name]}"
    echo "[OK] ${name}.sif save"
    echo ""
done

echo "All images pulled successfully!"

#ls -1 "${SINGULARITY_DIR}"/*.sif

# ===== 生成 BWA-MEM2 索引 =====
singularity exec \
    --bind ${REF_DIR}:/ref \
    ${SINGULARITY_DIR}/bwa-mem2_2.3.sif \
    bwa-mem2 index /ref/hg38.fa
# 内存不足可以直接下载

ln -sfn "${RAW_NAME}"/hg38.fa* .

# 补 .fai 和 .dict
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools faidx /data/genome/genome.fa
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
    gatk CreateSequenceDictionary -R /data/genome/genome.fa
# chr1 长度比对（GRCh38 标准值）
grep -P '^chr1\t' genome/genome.fa.fai
#   期望: chr1  248956422 ...

grep -P '^chr1\t' reference/raw/Homo_sapiens_assembly38.fasta.fai
#   期望: chr1  248956422 ...

# 确认 dbsnp VCF 的主染色体 contig（应列出 chr1 ... chrM）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
    bcftools view -h /data/reference/raw/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    | grep '^##contig=<ID=chr' | sed 's/.*ID=//;s/,.*//' | head -30

echo "Reference preparation complete!"
```

---

## 笔记《1-原始数据 QC + 修剪》  （共 13 个代码块，其中 5 个含流程命令）


### 1-原始数据 QC + 修剪 · 代码块 0

```bash
mkdir -p qc/fastqc_raw
# 执行fastqc
ls 20260720/data/raw_data/*/*_R{1,2}.fastq.gz | \
  xargs singularity exec --bind "$PWD:/data" singularity/fastqc_0.12.1.sif \
    fastqc -t 4 -o /data/qc/fastqc_raw
```

### 1-原始数据 QC + 修剪 · 代码块 1  <!-- 讲解/输出示例块 -->

```bash
# fastqc结果快速检查
for z in *.zip; do
  n="${z%.zip}"
  echo "── $n"
  unzip -p "$z" "$n/fastqc_data.txt" | \
    awk '/^>>/{m=$0; sub(/^>>/,"",m); sub(/\t.*/,"",m); s=$0; sub(/^.*\t/,"",s);
         if(m=="END_MODULE") next;
         if(m=="Basic Statistics"||m=="Per base sequence quality"||m=="Overrepresented sequences"||m=="Adapter Content"||s!="PASS") print "  "m" = "s}'
done
```

### 1-原始数据 QC + 修剪 · 代码块 2  <!-- 讲解/输出示例块 -->

```bash
遍历 26 个 *_fastqc.zip
  └→ 从 zip 里提取 fastqc_data.txt（纯文本摘要）
      └→ awk 过滤：只显示 PASS 以外的行 + 几个重点模块

awk命令详解：
awk '/^>>/{                           # 只处理 ">>" 开头的行
    m=$0; sub(/^>>/,"",m);            # 砍掉 ">>" 前缀 → m = "Per base sequence quality PASS"
    sub(/\t.*/,"",m);                 # 砍掉制表符及之后 → m = "Per base sequence quality"
    s=$0; sub(/^.*\t/,"",s);          # 砍掉制表符之前 → s = "PASS"

    if(m=="END_MODULE") next;         # 跳过终止标记行

    if(m=="Basic Statistics"          # 这几项始终显示（最有价值）
       || m=="Per base sequence quality"
       || m=="Overrepresented sequences"
       || m=="Adapter Content"
       || s!="PASS")                  # 其他项只在非 PASS 时显示
        print "  "m" = "s
}'

    是 Basic Statistics？ ──yes──→ 打印
    是 Per base sequence quality？ ──yes──→ 打印
    是 Overrepresented sequences？ ──yes──→ 打印
    是 Adapter Content？ ──yes──→ 打印
    是 END_MODULE？ ──yes──→ 跳过
    s != "PASS"？ ──yes──→ 打印（自动捕获所有 WARN/FAIL）
    s == "PASS" ──→ 不打印（忽略，不污染输出）

没有 WARN/FAIL 的模块不显示，扫一眼就能定位问题
结果如下：
── NA12878_R1_fastqc
  Basic Statistics = pass
  Per base sequence quality = pass
  Per sequence quality scores = pass
  Per base sequence content = pass
  Per sequence GC content = warn
  Per base N content = pass
  Sequence Length Distribution = pass
  Sequence Duplication Levels = pass
  Overrepresented sequences = pass
  Adapter Content = warn
── NA12878_R2_fastqc
  Basic Statistics = pass
  Per base sequence quality = pass
  Per sequence quality scores = pass
  Per base sequence content = pass
  Per sequence GC content = warn
  Per base N content = pass
  Sequence Length Distribution = pass
  Sequence Duplication Levels = pass
  Overrepresented sequences = pass
  Adapter Content = warn
```

### 1-原始数据 QC + 修剪 · 代码块 3  <!-- 讲解/输出示例块 -->

```bash
# 污染检查
## 成对（同一文件 R1+R2 一起查），覆盖三类样本
for n in CB261007720-OSW_R1_fastqc CB261007720-OSW_R2_fastqc \
         L200100435_R1_fastqc L200100435_R2_fastqc \
         NA12878_R1_fastqc NA12878_R2_fastqc; do
  echo "── $n"
  unzip -p "$n.zip" "$n/fastqc_data.txt" | \
    sed -n '/>>Adapter Content/,/>>END_MODULE/p' | \
    grep -v "^#" | grep -v "^>>" | \
    tail -6 | \
    awk -F'\t' '{printf "  %s\t%.2f%%\n", $1, $2}'
done
```

### 1-原始数据 QC + 修剪 · 代码块 4  <!-- 讲解/输出示例块 -->

```bash
── CB261007720-OSW_R1_fastqc
  128-129	3.31%
  130-131	3.60%
  132-133	3.92%
  134-135	4.26%
  136-137	4.61%
  138-139	4.96%
── CB261007720-OSW_R2_fastqc
  128-129	3.40%
  130-131	3.70%
  132-133	4.02%
  134-135	4.36%
  136-137	4.72%
  138-139	5.08%
── L200100435_R1_fastqc
  128-129	3.97%
  130-131	4.32%
  132-133	4.69%
  134-135	5.08%
  136-137	5.49%
  138-139	5.91%
── L200100435_R2_fastqc
  128-129	4.08%
  130-131	4.44%
  132-133	4.82%
  134-135	5.21%
  136-137	5.63%
  138-139	6.06%
── NA12878_R1_fastqc
  128-129	3.90%
  130-131	4.23%
  132-133	4.58%
  134-135	4.93%
  136-137	5.31%
  138-139	5.71%
── NA12878_R2_fastqc
  128-129	3.97%
  130-131	4.30%
  132-133	4.65%
  134-135	5.01%
  136-137	5.39%
  138-139	5.79%
```

### 1-原始数据 QC + 修剪 · 代码块 5  <!-- 讲解/输出示例块 -->

```bash
检查逻辑速览：

全量快速扫描 → 发现 末端 Adapter Content WARN
                  ↓
          查具体数值（本条命令）→ 确认末端 5.8~6.1% 超标
                  ↓
           决策：开启 fastp 修剪（实际上没有污染也推荐fastp修剪）
```

### 1-原始数据 QC + 修剪 · 代码块 6

```bash
fastqc -a 自定义接头文件.txt 输入.fastq.gz
```

### 1-原始数据 QC + 修剪 · 代码块 7  <!-- 讲解/输出示例块 -->

```bash
# 自定义接头文件（adapter.txt）
MyCustomAdapter_Forward	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
MyCustomAdapter_Reverse	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

### 1-原始数据 QC + 修剪 · 代码块 8  <!-- 讲解/输出示例块 -->

```bash
不修剪的 read：  [真实序列 140bp] + [接头序列 10bp]
                                                    ↓
BWA-MEM2 比对到参考基因组：这 10bp 接头在基因组上找不到匹配
  → 末尾被 soft-clip（切掉）→ CIGAR 变成 140M10S
  → 虽不影响 mapped 率，但增加比对模糊度
```

### 1-原始数据 QC + 修剪 · 代码块 9  <!-- 讲解/输出示例块 -->

```bash
接头序列在基因组上无同源区域，通常不会导致假 SNP。
但 soft-clip 末端 reads 质量偏低，HC 会降低该区域的 genotyping 置信度。
极端情况：接头序列碰巧和基因组某处部分匹配 → 假比对 → 假 SNP
```

### 1-原始数据 QC + 修剪 · 代码块 10

```bash
mkdir -p fastq_clean qc/fastp qc/fastqc_trim

while read -r sm r1 r2; do
  echo ">>> fastp $sm"
  singularity exec --bind "$PWD:/data" singularity/fastp_1.3.6.sif \
    fastp -i "/data/$r1" -I "/data/$r2" \
    -o "/data/fastq_clean/${sm}_R1.fastq.gz" -O "/data/fastq_clean/${sm}_R2.fastq.gz" \
    -h "/data/qc/fastp/${sm}.html" -j "/data/qc/fastp/${sm}.json" \
    --length_required 36 --thread 4
done < samples.tsv
```

### 1-原始数据 QC + 修剪 · 代码块 11

```bash
while read -r sm r1 r2; do          # ① 循环读清单，每行取 3 列
  echo ">>> fastp $sm"              # ② 打印进度提示
  singularity exec ... fastp ...    # ③ 容器里跑 fastp
done < samples.tsv                  # ④ 从清单文件读入

读取 R1/R2
  ↓
自动检测接头（从 overlap 推断）
  ↓
切除 read 末端接头序列
  ↓
质量过滤（去掉低质量碱基/reads）
  ↓
长度过滤（<36bp 丢弃）
  ↓
输出 clean fastq + html + json 报告
```

### 1-原始数据 QC + 修剪 · 代码块 12

```bash
# post-trim FastQC 复检（抽 4 个文件验证 adapter 已清除）
singularity exec --bind "$PWD:/data" singularity/fastqc_0.12.1.sif \
  fastqc -t 4 -o /data/qc/fastqc_trim \
  /data/fastq_clean/NA12878_R1.fastq.gz /data/fastq_clean/NA12878_R2.fastq.gz \
  /data/fastq_clean/CB261007723-OSW_R1.fastq.gz /data/fastq_clean/CB261007723-OSW_R2.fastq.gz

cd qc/fastqc_trim
for z in *.zip; do
  n="${z%.zip}"
  echo "── $n"
  unzip -p "$z" "$n/fastqc_data.txt" | \
    awk '/^>>/{m=$0; sub(/^>>/,"",m); sub(/\t.*/,"",m); s=$0; sub(/^.*\t/,"",s);
         if(m=="END_MODULE") next;
         if(m=="Basic Statistics"||m=="Per base sequence quality"||m=="Adapter Content"||s!="PASS") print "  "m" = "s}'
done
```

---

## 笔记《2-比对》  （共 17 个代码块，其中 8 个含流程命令）


### 2-比对 · 代码块 0  <!-- 讲解/输出示例块 -->

```bash
# 生成samples.tsv
awk -F, 'NR>1{print $3"\t20260720/data/raw_data/"$3"/"$3"_R1.fastq.gz\t20260720/data/raw_data/"$3"/"$3"_R2.fastq.gz"}' \
  20260720/data/samplesheet.csv > samples.tsv
```

### 2-比对 · 代码块 1  <!-- 讲解/输出示例块 -->

```bash
2S 20M 1D 2I 20M
```

### 2-比对 · 代码块 2  <!-- 讲解/输出示例块 -->

```bash
M + I + S 的总和 = SEQ 字段的碱基数
```

### 2-比对 · 代码块 3

```bash
# 单独流程验证（单样本测试）

## ① bwa-mem2 比对 → 管道直接 samtools sort
mkdir -p bam/NA12878 qc/flagstat qc/stats

singularity exec --bind "$PWD:/data" singularity/bwa-mem2_2.3.sif \
  bwa-mem2 mem -t 8 -K 100000000 -Y \
  -R "@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA" \
  /data/reference/genome/genome.fa \
  /data/fastq_clean/NA12878_R1.fastq.gz /data/fastq_clean/NA12878_R2.fastq.gz \
  | singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools sort -@ 8 -m 4G -O BAM -o /data/bam/NA12878/NA12878.sort.bam -
```

### 2-比对 · 代码块 4  <!-- 讲解/输出示例块 -->

```sam
@RG	ID:NA12878	SM:NA12878	PL:ILLUMINA
```

### 2-比对 · 代码块 5

```bash
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools index /data/bam/NA12878/NA12878.sort.bam
```

### 2-比对 · 代码块 6

```bash
# 随机访问一下，能秒回就说明索引生效
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools view -c /data/bam/NA12878/NA12878.sort.bam chr1:1000000-9000000

# 26500
# 索引目标区域的reads数
```

### 2-比对 · 代码块 7

```bash
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools flagstat /data/bam/NA12878/NA12878.sort.bam > /data/qc/flagstat/NA12878.sorted.flagstat

cat qc/flagstat/NA12878.sorted.flagstat

# 结果如下
# 格式说明：每行是"QC-passed + QC-failed"两个数，第二个几乎总是 0，直接看第一个数就行
6813859 + 0 in total (QC-passed reads + QC-failed reads)
6733928 + 0 primary
0 + 0 secondary
79931 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
6813270 + 0 mapped (99.99% : N/A)            <- 观察这个
6733339 + 0 primary mapped (99.99% : N/A)
6733928 + 0 paired in sequencing
3366964 + 0 read1
3366964 + 0 read2
6450400 + 0 properly paired (95.79% : N/A)    <- 观察这个
6732850 + 0 with itself and mate mapped
489 + 0 singletons (0.01% : N/A)              <- 观察这个
242774 + 0 with mate mapped to a different chr
217338 + 0 with mate mapped to a different chr (mapQ>=5)
```

### 2-比对 · 代码块 8  <!-- 讲解/输出示例块 -->

```bash
6450400 + 0 properly paired (95.79% : N/A)
...
242774 + 0 with mate mapped to a different chr   ← 占 paired 的 ~3.6%
```

### 2-比对 · 代码块 9

```bash
# 跨染色体匹配验证
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools view /data/bam/NA12878/NA12878.sort.bam | \
  awk '$7 != "=" && $7 != "*" {print $3, $7}' | \
  sort | uniq -c | sort -rn | head -20
```

### 2-比对 · 代码块 10  <!-- 讲解/输出示例块 -->

```bash
   2137 chr2 chr1
   2115 chr1 chr2
   1669 chr3 chr1
   1603 chr1 chr6
   1591 chr6 chr1
   1586 chr1 chr3
   1581 chr4 chr1
   1572 chr1 chr4
   1517 chr7 chr1
   1488 chr1 chr7
   1484 chr5 chr1
   1484 chr2 chr6
   1478 chr6 chr2
   1477 chr3 chr2
   1467 chr5 chr2
   1461 chr2 chr3
   1461 chr1 chr5
   1460 chr2 chr5
   1406 chr4 chr2
   1391 chr7 chr2
```

### 2-比对 · 代码块 11

```bash
# 质量分布检查
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools view -F 2304 /data/bam/NA12878/NA12878.sort.bam | \
  awk '$7 != "=" && $7 != "*" {print $5}' | sort -n | uniq -c

# 分布检查增加占比可视化
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
  samtools view -F 2304 /data/bam/NA12878/NA12878.sort.bam | \
  awk '$7 != "=" && $7 != "*" { c[$5]++; t++ }
       END {
         max = 0
         for (q = 0; q <= 60; q++) if (c[q] > max) max = c[q]
         for (q = 0; q <= 60; q++) {
           if (c[q] == 0) continue
           n = int(c[q] * 50 / max)
           bar = ""
           for (i = 0; i < n; i++) bar = bar "#"
           printf "MAPQ=%2d  %7d  %5.1f%%  %s\n", q, c[q], c[q] * 100 / t, bar
         }
       }'
```

### 2-比对 · 代码块 12  <!-- 讲解/输出示例块 -->

```bash
# 结果
MAPQ= 0    22841    9.4%  #####
MAPQ= 1      658    0.3%  
MAPQ= 2      547    0.2%  
MAPQ= 3      725    0.3%  
...(省略)
MAPQ=57      291    0.1%  
MAPQ=58      411    0.2%  
MAPQ=59      280    0.1%  
MAPQ=60   193582   79.7%  ##################################################
```

### 2-比对 · 代码块 13  <!-- 讲解/输出示例块 -->

```bash
chr1 拷贝 A:  ...ACGTACGT [G] ACGTACGT ACGT...
chr2 拷贝 B:  ...ACGTACGT [T] ACGTACGT ACGT...
                    ↑ 完全一致区段    ↑ 差异位点（唯一的"指纹"）
```

### 2-比对 · 代码块 14

```bash
singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools stats  /data/bam/NA12878/NA12878.sort.bam > /data/qc/stats/NA12878.sorted.samtools.stats
```

### 2-比对 · 代码块 15  <!-- 讲解/输出示例块 -->

```bash
SN	raw total sequences:	6813859
SN	reads mapped:	6813270
SN	reads properly paired:	6450400
SN	reads unmapped:	589
SN	average length:	150
SN	average quality:	36.1
SN	insert size average:	290.5
SN	insert size standard deviation:	30.2
SN	reads duplicated:	0
...
```

### 2-比对 · 代码块 16

```bash
# 批量处理
mkdir -p qc/flagstat qc/stats

while read -r sm r1 r2; do
  mkdir -p "bam/$sm"
  if [ -f "bam/$sm/$sm.sort.bam" ] && [ -f "bam/$sm/$sm.sort.bam.bai" ]; then
    echo "[SKIP] $sm 已存在（含索引）"
  else
    echo ">>> 比对 $sm"
    singularity exec --bind "$PWD:/data" singularity/bwa-mem2_2.3.sif \
      bwa-mem2 mem -t 8 -K 100000000 -Y \
      -R "@RG\tID:${sm}\tSM:${sm}\tPL:ILLUMINA" \
      /data/reference/genome/genome.fa "/data/$r1" "/data/$r2" \
      | singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
        samtools sort -@ 8 -m 4G -O BAM -o "/data/bam/$sm/$sm.sort.bam" -
    singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
      samtools index "/data/bam/$sm/$sm.sort.bam"
  fi


  # 统计分析
  singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools flagstat "/data/bam/$sm/$sm.sort.bam" > "/data/qc/flagstat/${sm}.sorted.flagstat"
  singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools stats  "/data/bam/$sm/$sm.sort.bam" > "/data/qc/stats/${sm}.sorted.samtools.stats"
done < samples.tsv

# 批量汇总检查（全样本 mapped% / properly paired% 一览），字段匹配
for f in qc/flagstat/*.sorted.flagstat; do
  sm=$(basename "$f" .sorted.flagstat)
  awk '
    /in total/     { total=$1 }
    $4=="mapped"   { mapped=$1 }
    $4=="paired"   { paired=$1 }
    $4=="properly" { pp=$1 }
    END { printf "%s: mapped %.2f%%  properly_paired %.2f%%\n",
          sm, mapped/total*100, pp/paired*100 }
  ' sm="$sm" "$f"
done
```

---

## 笔记《3-去重+校准》  （共 16 个代码块，其中 6 个含流程命令）


### 3-去重+校准 · 代码块 0

```bash
while read -r sm _r1 _r2; do
  IN_BAM="/data/bam/$sm/$sm.sort.bam"
  OUT_BAM="/data/bam/$sm/$sm.markdup.bam"
  METRICS="/data/bam/$sm/$sm.markdup.metrics"
  FLAGSTAT="/data/qc/flagstat/${sm}.markdup.flagstat"
  STATS="/data/qc/stats/${sm}.markdup.samtools.stats"

  # ① 断点续跑：BAM 与 METRICS 都齐全才算"完成"
  if [[ ! -f "$BAM" || ! -f "$METRICS" ]]; then
    echo ">>> MarkDuplicates $sm"
    singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
      gatk --java-options -Xmx2g MarkDuplicates \
      -I "$IN_BAM" \
      -O "$OUT_BAM" \
      -M "$METRICS" \
      --CREATE_INDEX true
  else
    echo "[SKIP] $sm 已存在（bam + metrics 均齐全）"
  fi

  # ② 统计阶段也加幂等
  if [[ ! -f "$FLAGSTAT" ]]; then
    singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
      samtools flagstat "$BAM" > "$FLAGSTAT"
  fi
  if [[ ! -f "$STATS" ]]; then
    singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
      samtools stats "$BAM" > "$STATS"
  fi

done < samples.tsv
```

### 3-去重+校准 · 代码块 1  <!-- 讲解/输出示例块 -->

```bash
samples.tsv ──循环──> sort.bam ──MarkDuplicates──> markdup.bam + markdup.metrics
                                                       │
                              ┌────────────────────────┴──────────┐
                              ▼                                    ▼
                     qc/flagstat/*.markdup.flagstat    qc/stats/*.markdup.samtools.stats
```

### 3-去重+校准 · 代码块 2  <!-- 讲解/输出示例块 -->

```bash
BIN          CoverageMult    all_sets    non_optical_sets
覆盖倍数     相对平均倍数    分子数       排除光学重复后的分子数
```

### 3-去重+校准 · 代码块 3  <!-- 讲解/输出示例块 -->

```bash
重复率高（>30%）？
├── ESTIMATED_LIBRARY_SIZE 足够大（≫ 测序对数）？
│   └── ✅ 只是测太深，降采样或接受，无需重建
└── ESTIMATED_LIBRARY_SIZE 偏小？
    └── ❌ 文库复杂度不足（原始分子太少）
        └── 重建文库 → 重新建库 + 重测序（唯一有效出路）
```

### 3-去重+校准 · 代码块 4  <!-- 讲解/输出示例块 -->

```bash
{
  printf "SAMPLE\tREAD_PAIRS\tDUP%%\tELS\tELS/Pairs\tOptical\n"
  for m in bam/*/*.markdup.metrics; do
    sm=$(basename "$(dirname "$m")")
    awk -F'\t' '
      /^##/ {next}
      /^#/  {next}
      $1=="LIBRARY" {hdr=1; next}
      hdr==1 {
        printf "%s\t%s\t%.1f%%\t%.0f\t%.2f\t%s\n",
               "'"$sm"'", $3, $9*100, $10, $10/$3, $8;
        exit
      }' "$m"
  done
} | sort -t$'\t' -k3,3n | column -t
```

### 3-去重+校准 · 代码块 5  <!-- 讲解/输出示例块 -->

```bed
chr11   35067210        35067370        chr11:35067210:35067370:chr11_35067210_35067370 160     +
chr11   36314632        36314792        chr11:36314632:36314792:chr11_36314632_36314792 160     +
chr11   36365270        36365430        chr11:36365270:36365430:chr11_36365270_36365430 160     +
chr11   36410393        36410553        chr11:36410393:36410553:chr11_36410393_36410553 160     +
chr1    1238150 1238310 chr1:1238150:1238310:chr1_1238150_1238310       160     +
chr1    1243464 1243624 chr1:1243464:1243624:chr1_1243464_1243624       160     +
chr1    1312033 1312193 chr1:1312033:1312193:chr1_1312033_1312193       160     +
chr1    2552441 2552601 chr1:2552441:2552601:chr1_2552441_2552601       160     +
```

### 3-去重+校准 · 代码块 6  <!-- 讲解/输出示例块 -->

```bash
# 生成排序 BED
awk 'NF>=3{print $1"\t"$2"\t"$3}' targets.bed | sort -k1,1V -k2,2n -u > ./reference/targets.sorted.bed
echo "targets.sorted.bed 行数: $(wc -l < ./reference/targets.sorted.bed)"
```

### 3-去重+校准 · 代码块 7

```bash
# 把 3 列的 BED 转成带 @HD/@SQ header 的 interval_list
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx2g BedToIntervalList \
  -I /data/reference/targets.sorted.bed \
  -O /data/reference/targets.sorted.interval_list \
  -SD /data/reference/genome/genome.dict
echo "interval_list 头部:"; head -4 reference/targets.sorted.interval_list
```

### 3-去重+校准 · 代码块 8  <!-- 讲解/输出示例块 -->

```bash
targets.bed（6列）
   │  awk 只取前3列
   ▼
(3列 BED)
   │  sort -k1,1V -k2,2n -u（版本排序+去重）
   ▼
./reference/targets.sorted.bed
   │  BedToIntervalList -SD genome.dict（加 header）
   ▼
targets.sorted.interval_list ──> CollectHsMetrics --BAIT/--TARGET
                                     samtools/mosdepth 吃 targets.bed 也可
```

### 3-去重+校准 · 代码块 9

```bash
mkdir -p qc/hsmetrics

singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx4g CollectHsMetrics \
  -I /data/bam/NA12878/NA12878.markdup.bam \
  -O /data/qc/hsmetrics/NA12878.markdup.hsmetrics.txt \
  --BAIT_INTERVALS /data/reference/targets.sorted.interval_list \
  --TARGET_INTERVALS /data/reference/targets.sorted.interval_list \
  -R /data/reference/genome/genome.fa
```

### 3-去重+校准 · 代码块 10  <!-- 讲解/输出示例块 -->

```bash
awk -F'\t' '
  /^## METRICS/ {next}
  /^BAIT_SET/ {
    for (i=1;i<=NF;i++) h[$i]=i
    next
  }
  $1=="targets" && NF>=60 {
    pq=$(h["PF_UQ_BASES_ALIGNED"]) # 去重后口径，非去重口径可用$(h["PF_BASES_ALIGNED"])
    printf "PCT_SELECTED %.4f\n",  $(h["PCT_SELECTED_BASES"])
    printf "MEAN_CVG      %.1f\n", $(h["MEAN_TARGET_COVERAGE"])
    printf "MED_CVG       %.1f\n", $(h["MEDIAN_TARGET_COVERAGE"])
    printf "FOLD          %.1f\n", $(h["FOLD_ENRICHMENT"])
    printf "20X           %.4f\n", $(h["PCT_TARGET_BASES_20X"])
    printf "on-bait       %.4f\n", $(h["ON_BAIT_BASES"])/pq
    printf "on-target     %.4f\n", $(h["ON_TARGET_BASES"])/pq
    printf "off-target    %.4f\n", $(h["OFF_BAIT_BASES"])/pq
  }' qc/hsmetrics/NA12878.markdup.hsmetrics.txt
```

### 3-去重+校准 · 代码块 11  <!-- 讲解/输出示例块 -->

```bash
PCT_SELECTED 0.2590
MEAN_CVG      82.8
MED_CVG       84.0
FOLD          433.7
20X           0.9447
on-bait       0.1766
on-target     0.1057
off-target    0.9298
```

### 3-去重+校准 · 代码块 12

```bash
mkdir -p qc/hsmetrics

while read -r sm _r1 _r2; do
  if [[ -f qc/hsmetrics/$sm.markdup.hsmetrics.txt ]]; then
    echo "[SKIP] $sm 已有 hsmetrics"
  else
    echo ">>> HsMetrics $sm"
    singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
      gatk --java-options -Xmx4g CollectHsMetrics \
      -I /data/bam/$sm/$sm.markdup.bam \
      -O /data/qc/hsmetrics/$sm.markdup.hsmetrics.txt \
      --BAIT_INTERVALS /data/reference/targets.sorted.interval_list \
      --TARGET_INTERVALS /data/reference/targets.sorted.interval_list \
      -R /data/reference/genome/genome.fa
  fi
done < samples.tsv

{
  printf "SAMPLE\tPCT_SELECTED\tMEAN_CVG\tMED_CVG\tFOLD\t20X\ton-bait\ton-target\toff-target\n"
  for m in qc/hsmetrics/*.markdup.hsmetrics.txt; do
    sm=$(basename "$m" .markdup.hsmetrics.txt)
    awk -F'\t' -v sm="$sm" '
      /^## METRICS/ {next}
      /^BAIT_SET/ {
        for (i=1;i<=NF;i++) h[$i]=i
        next
      }
      $1=="targets" && NF>=60 && !seen {
        pq=$(h["PF_UQ_BASES_ALIGNED"])
        printf "%s\t%.4f\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\n",
               sm,
               $(h["PCT_SELECTED_BASES"]), $(h["MEAN_TARGET_COVERAGE"]),
               $(h["MEDIAN_TARGET_COVERAGE"]), $(h["FOLD_ENRICHMENT"]),
               $(h["PCT_TARGET_BASES_20X"]),
               $(h["ON_BAIT_BASES"])/pq, $(h["ON_TARGET_BASES"])/pq,
               $(h["OFF_BAIT_BASES"])/pq
        seen=1
      }' "$m"
  done
} | sort -t$'\t' -k8,8n | column -t
```

### 3-去重+校准 · 代码块 13

```bash
# BQSR（碱基质量校准）
while read -r sm _r1 _r2; do
  # ① 拟合误差模型（用已知位点）
  if [ ! -f "bam/$sm/$sm.recal.table" ]; then
    echo ">>> BaseRecalibrator $sm"
    singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
      gatk --java-options -Xmx2g BaseRecalibrator \
      -R /data/reference/genome/genome.fa \  # 参考基因组
      -I "/data/bam/$sm/$sm.markdup.bam" \  # 输入=去重后的 BAM
      --known-sites /data/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz \  # 标准答案册①：遮罩已知 SNP
      --known-sites /data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \  # 标准答案册②：遮罩已知 indel
      -O "/data/bam/$sm/$sm.recal.table"  # 产物=修正表
  fi
  # ② 按上一步输出的修正表模型重写 BAM 里的碱基质量
  if [ ! -f "bam/$sm/$sm.markdup.BQSR.bam" ]; then
    echo ">>> ApplyBQSR $sm"
    singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
      gatk --java-options -Xmx2g ApplyBQSR \
      -R /data/reference/genome/genome.fa \
      -I "bam/$sm/$sm.markdup.bam" \  # 输入=去重后的 BAM（和1相同）
      --bqsr-recal-file "bam/$sm/$sm.recal.table" \  # 读①产出的修正表
      -O "bam/$sm/$sm.markdup.BQSR.bam"  # 产物=校准后 BAM
  fi
  # recal 阶段统计（三次 stats 之第三套，进 MultiQC）
  singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools flagstat "bam/$sm/$sm.markdup.BQSR.bam" > "qc/flagstat/${sm}.recal.flagstat"
  singularity exec --bind "$PWD:/data" singularity/samtools_1.24.sif \
    samtools stats  "bam/$sm/$sm.markdup.BQSR.bam" > "qc/stats/${sm}.recal.samtools.stats"
done < samples.tsv
```

### 3-去重+校准 · 代码块 14  <!-- 讲解/输出示例块 -->

```bash
# flagstat 逐行 diff
cd qc/flagstat
for f in *.markdup.flagstat; do
  sm=${f%.markdup.flagstat}                      # 去掉后缀拿样本名
  if diff -q "$f" "${sm}.recal.flagstat" >/dev/null 2>&1; then
    echo "PASS  $sm"
  else
    echo "DIFF  $sm"
    diff "$f" "${sm}.recal.flagstat"
  fi
done
```

### 3-去重+校准 · 代码块 15

```bash
cd ~/nvme2/work/pipeline/gwas
singularity exec --bind $PWD:/data singularity/samtools_1.24.sif \
  samtools quickcheck -v bam/*/*.markdup*.bam && echo 'ALL_BAM_OK'
```

---

## 笔记《4-变异检测》  （共 23 个代码块，其中 21 个含流程命令）


### 4-变异检测 · 代码块 0  <!-- 讲解/输出示例块 -->

```json
BQSR.bam → HaplotypeCaller（-ERC GVCF）→ 每个样本一个 .g.vcf.gz
                                          ↓ 小样本不用单样本 VCF
                            4 个样本 .g.vcf.gz → CombineGVCFs → GenotypeGVCFs
                                          ↓
                                联合分型 VCF → 硬过滤 → 最终 VCF
```

### 4-变异检测 · 代码块 1

```bash
mkdir -p gvcf

i=0
while read -r sm _r1 _r2; do
  [ -f "gvcf/$sm.g.vcf.gz" ] && echo "[SKIP] $sm" && continue
  echo ">>> HC $sm"
  ( singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
      gatk --java-options -Xmx2g HaplotypeCaller \  # JVM 堆内存 2G（panel 数据够了，默认 4G 是浪费）
      -R /data/reference/genome/genome.fa \  # 	参考基因组
      -I "/data/bam/$sm/$sm.markdup.BQSR.bam" \  # BQSR 后的干净 BAM
      -O "/data/gvcf/$sm.g.vcf.gz" \  # 输出 gVCF（gzip 压缩）
      -ERC GVCF  \  # gVCF 模式，存证据不拍板
      -L /data/reference/targets.sorted.bed \
      --interval-padding 100 \  # 	只 call 目标区间 ± 100bp
      --native-pair-hmm-threads 4 \  # "配对 HMM 打分"用 4 线程
      -ploidy 2  # 人类二倍体（默认就是 2，显式写出防误用）
  ) &  # 把 HC 放进子 shell 后台——让多个样本同时跑
  i=$((i+1))  # 计数器到 3 就 wait 同步一次：一批 3 个跑完再放下一批，控制并发数
  [ $((i % 3)) = 0 ] && wait && echo "--- 一批完成 ---"
done < samples.tsv
wait
echo "完成: $(ls gvcf/*.g.vcf.gz 2>/dev/null | wc -l) / 12"
```

### 4-变异检测 · 代码块 2  <!-- 讲解/输出示例块 -->

```bash
12 个 gVCF ──CombineGVCFs──▶ 1 个合并 gVCF ──GenotypeGVCFs──▶ cohort.raw.vcf.gz
  gvcf/          "证据拼桌"      cohort/           "统一拍板"       (12样本全位点 GT)
```

### 4-变异检测 · 代码块 3

```bash
mkdir -p cohort
ls gvcf/*.g.vcf.gz > gvcf.list
echo "gVCF 数量: $(wc -l < gvcf.list)"

singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx4g CombineGVCFs \
  -R /data/reference/genome/genome.fa \
  -V gvcf.list \
  -O cohort/cohort.g.vcf.gz
```

### 4-变异检测 · 代码块 4

```bash
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx4g GenotypeGVCFs \
  -R /data/reference/genome/genome.fa \
  -V cohort/cohort.g.vcf.gz \
  -O cohort/cohort.raw.vcf.gz
```

### 4-变异检测 · 代码块 5

```bash
## ① 样本数（应符合样本总数，此次运行示例应为12）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools query -l cohort/cohort.raw.vcf.gz | wc -l


## ② 变异总数（1MB 区域预期几百到几千）：列变异记录 → 数条数
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -H cohort/cohort.raw.vcf.gz | wc -l
# 结果：17090

## ③ SNP vs INDEL 比例粗看（REF 长度 1=SNP）：量 REF 长度 → 分类计数
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -H cohort/cohort.raw.vcf.gz \
  | awk '{print length($4)}' \
  | sort -n \
  | uniq -c \
  | sort -k2,2n \
  | awk 'BEGIN{print "REF长度\t记录数\t粗分类型"}
         {if($2==1) t="SNP(或插入)";
          else if($2<=9) t="短indel";
          else t="长REF/复杂位点";
          print $2"\t"$1"\t"t}'
# 结果
REF长度 记录数  粗分类型
1       15752   SNP(或插入)
2       581     短indel
3       276     短indel
4       133     短indel
5       127     短indel
6       28      短indel
7       24      短indel
8       14      短indel
9       28      短indel
10      10      长REF/复杂位点
11      18      长REF/复杂位点
12      7       长REF/复杂位点
13      13      长REF/复杂位点
14      3       长REF/复杂位点
15      9       长REF/复杂位点
16      6       长REF/复杂位点
17      4       长REF/复杂位点
18      6       长REF/复杂位点
19      3       长REF/复杂位点
20      1       长REF/复杂位点
21      6       长REF/复杂位点
22      1       长REF/复杂位点
23      5       长REF/复杂位点
24      1       长REF/复杂位点
25      3       长REF/复杂位点
27      3       长REF/复杂位点
28      1       长REF/复杂位点
29      5       长REF/复杂位点
30      1       长REF/复杂位点
31      5       长REF/复杂位点
32      1       长REF/复杂位点
33      1       长REF/复杂位点
37      1       长REF/复杂位点
39      1       长REF/复杂位点
44      1       长REF/复杂位点
45      1       长REF/复杂位点
49      1       长REF/复杂位点
55      1       长REF/复杂位点
58      1       长REF/复杂位点
64      1       长REF/复杂位点
65      3       长REF/复杂位点
67      1       长REF/复杂位点
90      1       长REF/复杂位点
96      1       长REF/复杂位点
```

### 4-变异检测 · 代码块 6

```bash
## 长尾区间落点判断，整个 raw VCF 先切成"targets 区间内"，再数长 REF，若条数较低（占比较低）就不管了
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -H -T reference/targets.sorted.bed cohort/cohort.raw.vcf.gz \
  | wc -l
# 结果：10072


## 落点检查
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -T reference/targets.sorted.bed cohort/cohort.raw.vcf.gz \
  | awk 'length($4)>=30{print $1"\t"$2"\tREF_len="length($4)"\tALT="$5}'
# 结果
chr1    7829912 REF_len=55      ALT=G
chr1    163758726       REF_len=67      ALT=G
chr1    247416614       REF_len=64      ALT=G
chr6    32624621        REF_len=33      ALT=*,G
chr6    32624661        REF_len=96      ALT=A
chr7    1845497 REF_len=65      ALT=T
chr7    1845542 REF_len=65      ALT=*,C
chr7    1996826 REF_len=90      ALT=C
chr17   30237266        REF_len=44      ALT=A
```

### 4-变异检测 · 代码块 7

```bash
> # 生成外扩100bp的目标文件
> awk '{print $1"\t"($2-100)"\t"($3+100)}' reference/targets.sorted.bed > /tmp/targets.pad100.bed
> singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
>   bcftools view -H -T /tmp/targets.pad100.bed cohort/cohort.raw.vcf.gz | wc -l
> # 输出数量和变异总数一致或近似
>
```

### 4-变异检测 · 代码块 8

```bash
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -T reference/targets.sorted.bed cohort/cohort.raw.vcf.gz \
  | awk 'length($4)>=30'
```

### 4-变异检测 · 代码块 9

```bash
# 前置一步：把多等位/MIXED 记录拆成一行一个，否则会丢失多等位位点 + 左对齐 + REF 校验（-f 与数据同源）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools norm -f reference/genome/genome.fa -m -any cohort/cohort.raw.vcf.gz -Oz -o cohort/cohort.raw.split.vcf.gz
## bcftools norm -f <参考基因组> = 给规范化过程加一个"参考坐标系"

# 建索引（后续 SelectVariants / 多文件统计都要用）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
bcftools index -t cohort/cohort.raw.split.vcf.gz
```

### 4-变异检测 · 代码块 10

```bash
printf '\n%-26s %8s %6s %5s %7s %7s %6s %6s %8s\n' 文件 records SNPs MNPs indels others noALT multi 差额
for f in cohort/cohort.raw.vcf.gz cohort/cohort.raw.split.vcf.gz; do
  singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif bcftools stats "$f" | \
  awk -v name="${f##*/}" '/^SN/{
      if($0~/number of records:/)rec=$NF;
      else if($0~/number of SNPs:/)snp=$NF;
      else if($0~/number of MNPs:/)mnp=$NF;
      else if($0~/number of indels:/)ind=$NF;
      else if($0~/number of others:/)oth=$NF;
      else if($0~/number of no-ALTs:/)noal=$NF;
      else if($0~/number of multiallelic sites:/)ma=$NF}
    END{d=snp+mnp+ind+oth+noal-rec;
      printf "%-26s %8d %6d %5d %7d %7d %6d %6d %+8d\n", name, rec, snp, mnp, ind, oth, noal, ma, d}'
done
echo "差额 = 各类之和 − records；差额>0 = 同行多类型被重复计数(MIXED)；差额<0 = 有记录未被归入任何类别"
```

### 4-变异检测 · 代码块 11

```bash
# 差额>0 或 <0 时，按 ALT 形态分类，找"未归类"的记录类型
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
bcftools view -H cohort/cohort.raw.split.vcf.gz | \
awk '{ref=$4; alt=$5;
      if(alt==".") noalt++;
      else if(ref~/^[ACGT]$/ && alt~/^[ACGT]$/) snp++;
      else if(alt~/^\*$/) star++;
      else if(ref~/^[ACGT]+$/ && alt~/^[ACGT]+$/ && length(ref)!=length(alt)) indel++;
      else other++}
     END{print "SNP="snp, "indel="indel, "star(*)="star, "noALT="noalt, "other="other, "total="NR}'
```

### 4-变异检测 · 代码块 12

```bash
# SNP / INDEL 分开（阈值不同，必须分开处理）
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx2g SelectVariants -V cohort/cohort.raw.split.vcf.gz--select-type SNP -O cohort/cohort.snp.vcf.gz

singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx2g SelectVariants -V cohort/cohort.raw.split.vcf.gz --select-type INDEL -O cohort/cohort.indel.vcf.gz
```

### 4-变异检测 · 代码块 13

```bash
# SNP 硬过滤（GATK best practices 阈值）
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx2g VariantFiltration -V cohort/cohort.snp.vcf.gz -O cohort/cohort.snp.hardfiltered.vcf.gz \
  -filter "QD < 2.0" --filter-name QD2 \
  -filter "QUAL < 30.0" --filter-name QUAL30 \
  -filter "SOR > 3.0" --filter-name SOR3 \
  -filter "FS > 60.0" --filter-name FS60 \
  -filter "MQ < 40.0" --filter-name MQ40 \
  -filter "MQRankSum < -12.5" --filter-name MQRankSum \
  -filter "ReadPosRankSum < -8.0" --filter-name ReadPosRankSum
```

### 4-变异检测 · 代码块 14

```bash
# INDEL 硬过滤
singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
  gatk --java-options -Xmx2g VariantFiltration -V cohort/cohort.indel.vcf.gz -O cohort/cohort.indel.hardfiltered.vcf.gz \
  -filter "QD < 2.0" --filter-name QD2 \
  -filter "QUAL < 30.0" --filter-name QUAL30 \
  -filter "FS > 200.0" --filter-name FS200 \
  -filter "SOR > 10.0" --filter-name SOR10 \
  -filter "ReadPosRankSum < -20.0" --filter-name ReadPosRankSum
```

### 4-变异检测 · 代码块 15

```bash
# ④ 合并回单一 VCF + 建索引
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools concat -a cohort/cohort.snp.hardfiltered.vcf.gz cohort/cohort.indel.hardfiltered.vcf.gz \
  -Oz -o cohort/cohort.hardfiltered.vcf.gz
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools index -t cohort/cohort.hardfiltered.vcf.gz
```

### 4-变异检测 · 代码块 16

```bash
mkdir -p results qc/bcftools_stats

# 全口径：构建位点×样本基因型矩阵（全部位点；./. = 缺失 ≠ 0/0，分析时务必区分）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools query -R /data/reference/targets.sorted.bed \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  /data/cohort/cohort.hardfiltered.vcf.gz > results/genotype_matrix.tsv
echo "矩阵行数: $(wc -l < results/genotype_matrix.tsv)"
head -3 results/genotype_matrix.tsv
```

### 4-变异检测 · 代码块 17

```bash
# 前置检查，判断标签完整性
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools query -f '%FILTER\n' /data/cohort/cohort.hardfiltered.vcf.gz | sort | uniq -c | sort -rn
## 正常跑完 VariantFiltration 的记录，FILTER 列只会是 PASS 或某个标签名，不会出现 .。`. 的出现等于一条铁证——某条过滤命令漏跑或漏生效了

# PASS 子集：PASS 位点 + 质控详情（GT:DP:AD——下游关联分析的正式输入）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -f PASS -Oz -o /data/cohort/cohort.PASS.vcf.gz \
  /data/cohort/cohort.hardfiltered.vcf.gz  # 这个裸参数是输入文件
## 建立索引
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools index -t /data/cohort/cohort.PASS.vcf.gz

singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n' \
  -o /data/results/genotype_detail_PASS.tsv \
  /data/cohort/cohort.PASS.vcf.gz 

# 行数对账：detail 行数 应等于 PASS VCF 记录数
echo "detail 行数 : $(wc -l < results/genotype_detail_PASS.tsv)"
echo "PASS 记录数 : $(singularity exec --bind $PWD:/data singularity/bcftools_1.24.sif bcftools view -H cohort/cohort.PASS.vcf.gz | wc -l)"

# 验"模板展开对不对",列数对账：应为 5 + 3N
# 5 来自括号外的固定字段（CHROM / POS / ID / REF / ALT），3N 来自括号内 3 个字段 × N 个样本。
N=$(singularity exec --bind $PWD:/data singularity/bcftools_1.24.sif bcftools query -l cohort/cohort.PASS.vcf.gz | wc -l)
echo "期望列数 : $(( 5 + 3 * N ))"
head -1 results/genotype_detail_PASS.tsv | awk -F'\t' '{print "实际列数 :", NF}'

# 字段形态：DP 在第 7、10、13… 列，不该是字面量
awk -F'\t' '{for(i=7;i<=NF;i+=3) if($i=="DP") bad++}
            END{print "字面量 DP 格子数:", bad+0}' results/genotype_detail_PASS.tsv
```

### 4-变异检测 · 代码块 18

```bash
## VCF 正式统计（进 MultiQC）
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools stats /data/cohort/cohort.hardfiltered.vcf.gz > qc/bcftools_stats/cohort.hardfiltered.stats

singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools stats /data/cohort/cohort.PASS.vcf.gz > qc/bcftools_stats/cohort.PASS.stats

ls qc/bcftools_stats/
```

### 4-变异检测 · 代码块 19

```bash
# 格式检查-看 stats 到底输出了什么，单纯看看格式长什么样方便后续写命令
echo "=== stats 前 30 行（诊断）==="
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools stats /data/cohort/cohort.raw.vcf.gz | head -30
```

### 4-变异检测 · 代码块 20

```bash
# FILTER 标签分布：PASS vs FAIL 分布 + raw/PASS 的 Ti/Tv
echo "--- FILTER 标签分布 ---"
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
  bcftools view -H cohort/cohort.PASS.vcf.gz | awk '{print $7}' | sort | uniq -c | sort -rn
```

### 4-变异检测 · 代码块 21

```bash
# 对比Raw和硬过滤的Ti/Tv(转换数 / 颠换数)
## 全口径分析
echo '===== 全量 Ti/Tv ====='
for f in cohort/cohort.raw.vcf.gz cohort/cohort.PASS.vcf.gz; do
  printf '%-38s' "$f";
  singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif bcftools stats "$f" | awk '/^TSTV\t/{printf "ts=%s tv=%s Ti/Tv=%s (1stALT=%s)\n", $3, $4, $5, $8}';
done
## 单口径-纯双等位SNP位点
echo "===== 纯双等位 SNP Ti/Tv（单等位口径）====="
for spec in "raw 全量|cohort/cohort.raw.vcf.gz|" "PASS 子集|cohort/cohort.PASS.vcf.gz|"; do
  IFS='|' read -r name f extra <<< "$spec"
  printf '%-12s %-36s ' "$name" "$f"
  singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif \
    bcftools view $extra -H "$f" | \
  awk 'length($4)==1 && length($5)==1{a=$4;b=$5; if((a=="A"&&b=="G")||(a=="G"&&b=="A")||(a=="C"&&b=="T")||(a=="T"&&b=="C"))ts++; else tv++} END{printf "ts=%d  tv=%d  Ti/Tv=%.2f\n",ts,tv,ts/tv}'
done

# 诊断位点集合差异（有没有静默丢位点）
echo "===== ② 位点集合差异：raw vs hardfiltered（查丢位点/新增）====="
singularity exec --bind "$PWD:/data" singularity/bcftools_1.24.sif bcftools stats cohort/cohort.raw.vcf.gz cohort/cohort.PASS.vcf.gz | \
  awk '/^ID\t/{f[$2]=$3; if(NF>3) f[$2]=f[$2]" + "$4}
     /^TSTV\t/{printf "id %s (%s): ts=%s tv=%s Ti/Tv=%s\n", $2, f[$2], $3, $4, $5}'
```

### 4-变异检测 · 代码块 22

```bash
mkdir -p per_sample_vcf
BC="singularity exec --bind $PWD:/data $PWD/singularity/bcftools_1.24.sif bcftools"

while read -r sm _r1 _r2; do
  $BC view -s "$sm" -Oz -o "/data/per_sample_vcf/${sm}.hardfiltered.vcf.gz" /data/cohort/cohort.hardfiltered.vcf.gz
  $BC index -t "/data/per_sample_vcf/${sm}.hardfiltered.vcf.gz"
  $BC view -s "$sm" -Oz -o "/data/per_sample_vcf/${sm}.PASS.vcf.gz" /data/cohort/cohort.PASS.vcf.gz
  $BC index -t "/data/per_sample_vcf/${sm}.PASS.vcf.gz"
done < samples.tsv

# 验收
src=$($BC view -H cohort/cohort.PASS.vcf.gz | wc -l)
for f in per_sample_vcf/*.PASS.vcf.gz; do
  printf "%s  records=%s (src=%s)  samples=%s\n" \
    "$(basename "$f")" "$($BC view -H "$f" | wc -l)" "$src" "$($BC query -l "$f" | wc -l)"
done
```

---

## 笔记《5-测序质量》  （共 8 个代码块，其中 5 个含流程命令）


### 5-测序质量 · 代码块 0

```bash
mkdir -p qc/mosdepth

# mosdepth ×2（markdup 后 + BQSR 后各一次），3 任务一批并行
i=0
while read -r sm _r1 _r2; do
  for tag in md bqsr; do
    if [ "$tag" = "md" ]; then BAM="bam/$sm/$sm.markdup.bam"; else BAM="bam/$sm/$sm.markdup.BQSR.bam"; fi
    if [ -f "qc/mosdepth/${sm}.${tag}.mosdepth.summary.txt" ]; then echo "[SKIP] $sm.$tag"; continue; fi
    echo ">>> mosdepth $sm.$tag (后台)"
    (
      singularity exec --bind "$PWD:/data" singularity/mosdepth_0.3.14.sif \
        mosdepth --by /data/reference/targets.sorted.bed --fast-mode -t 4 \
        "/data/qc/mosdepth/${sm}.${tag}" "$BAM"
    ) &
    i=$((i+1))
    if [ $((i % 3)) -eq 0 ]; then wait; echo "--- 一批完成 ---"; fi
  done
done < samples.tsv
wait
echo "mosdepth 完成: $(ls qc/mosdepth/*.mosdepth.summary.txt | wc -l) 个 summary"
ls qc/mosdepth/*.regions.bed.gz | wc -l
ls qc/mosdepth/*.per-base.bed.gz | wc -l
```

### 5-测序质量 · 代码块 1  <!-- 讲解/输出示例块 -->

```python
import sys, gzip, os
matrix = sys.argv[1]
regions_dir = sys.argv[2]
samples = sys.argv[3].split()
out = sys.argv[4]
DP_MIN = 20

depth = {}
for sm in samples:
    f = f"{regions_dir}/{sm}.bqsr.regions.bed.gz"   # ★ 无 .mosdepth. 中缀！
    if not os.path.exists(f):
        print(f"[WARN] 缺 regions: {f}", file=sys.stderr)
        continue
    d = {}
    with gzip.open(f, 'rt') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            c, s, e = parts[0], int(parts[1]), int(parts[2])
            mean = float(parts[3])
            d[(c, s, e)] = mean
    depth[sm] = d
    print(f"[INFO] {sm}: {len(d)} regions")

n_total = 0; n_filled = 0; n_kept = 0
with open(matrix) as f, open(out, 'w') as fo:
    for line in f:
        p = line.rstrip('\n').split('\t')
        chrom, pos, ref, alt = p[:4]
        gts = p[4:]
        posi = int(pos)
        new_gts = []
        for sm, gt in zip(samples, gts):
            if gt == './.':
                n_total += 1
                mean = None
                if sm in depth:
                    for (c, s, e), m in depth[sm].items():
                        if c == chrom and s <= posi < e:
                            mean = m
                            break
                if mean is not None and mean >= DP_MIN:
                    new_gts.append('0/0')
                    n_filled += 1
                else:
                    new_gts.append('./.')
                    n_kept += 1
            else:
                new_gts.append(gt)
        fo.write('\t'.join(p[:4] + new_gts) + '\n')
print(f"done: {out}")
print(f"裁决: ./ 总计 {n_total} -> 填 0/0: {n_filled} · 保留 ./.: {n_kept}")
```

### 5-测序质量 · 代码块 2

```bash
# 需要 python3（服务器应有）；用 bcftools query 拿样本名顺序
SAMPLE_LIST=$(singularity exec --bind $PWD:/data singularity/bcftools_1.24.sif \
  bcftools query -l /data/cohort/cohort.PASS.vcf.gz | tr '\n' ' ')

python3 adjudicate_nocall.py \
  results/genotype_matrix.tsv qc/mosdepth "$SAMPLE_LIST" \
  results/genotype_matrix.adjudicated.tsv

# 验收：./. 数量应下降（达标部分转成 0/0）
echo "裁决前 ./.: $(awk -F'\t' '{for(i=5;i<=NF;i++) if($i=="./.") c++} END{print c}' results/genotype_matrix.tsv)"
echo "裁决后 ./.: $(awk -F'\t' '{for(i=5;i<=NF;i++) if($i=="./.") c++} END{print c}' results/genotype_matrix.adjudicated.tsv)"
awk -F'\t' '{print NF}' results/genotype_matrix.adjudicated.tsv | sort | uniq -c

# 分样本检查
for f in per_sample_vcf/*.GT.adjudicated.tsv; do
  printf "%s  ./.: %s\n" "$(basename "$f")" "$(grep -c $'\t\./\.' "$f")"
done
```

### 5-测序质量 · 代码块 3  <!-- 讲解/输出示例块 -->

```python
#!/usr/bin/env python3
"""
方案一（最终版）：裁决后矩阵 → 每样本合规 VCF（保留全 FORMAT，只替换 GT）
流程：
  1. bcftools view -s <sm> 从 cohort.PASS.vcf.gz 拆出该样本原始 VCF（保留全字段）
  2. 读取裁决后矩阵，取该样本列的新 GT
  3. 按 (CHROM,POS) 匹配，替换 VCF 每行的 GT 字段（FORMAT 第一个字段）
  4. 输出 + 压缩 + 索引
用法：
  python3 rebuild_vcf_python.py <sm> <gt_tsv> <in.vcf> <out.vcf>
"""
import sys, gzip

sm = sys.argv[1]
gt_tsv = sys.argv[2]
in_vcf = sys.argv[3]
out_vcf = sys.argv[4]

new_gt = {}
with open(gt_tsv) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 5:
            new_gt[(p[0], int(p[1]))] = p[4]
        elif len(p) == 3:
            new_gt[(p[0], int(p[1]))] = p[2]

n_total = 0
n_changed = 0
with open(in_vcf) as f, open(out_vcf, 'w') as fo:
    for line in f:
        if line.startswith('#'):
            fo.write(line)
            continue
        p = line.rstrip('\n').split('\t')
        chrom, pos = p[0], int(p[1])
        n_total += 1
        if len(p) >= 10:
            fmt = p[8].split(':')
            sm_cols = p[9].split(':')
            if fmt[0] == 'GT' and (chrom, pos) in new_gt:
                sm_cols[0] = new_gt[(chrom, pos)]
                p[9] = ':'.join(sm_cols)
                n_changed += 1
        fo.write('\t'.join(p) + '\n')

print(f"[{sm}] 总记录={n_total}  GT已改={n_changed}")
```

### 5-测序质量 · 代码块 4

```bash
#!/usr/bin/env bash
set -euo pipefail

BC="singularity exec --bind $PWD:/data singularity/bcftools_1.24.sif bcftools"
OUTDIR=per_sample_vcf
mkdir -p "$OUTDIR"

# 样本列映射：samples.tsv 顺序 = 矩阵列顺序（第1样本=矩阵列5）
declare -A COL=(
  [CB261007720-OSW]=5   [CB261007720-SAL]=6
  [CB261007723-OSW]=7   [CB261007723-SAL]=8
  [CB261007729-OSW]=9   [CB261007729-SAL]=10
  [CB261007733-OSW]=11  [CB261007733-SAL]=12
  [L200100435]=13       [L20260615001]=14
  [NA12878]=15          [NA18544]=16
)

for sm in "${!COL[@]}"; do
  col=${COL[$sm]}
  echo "===== $sm (矩阵列 $col) ====="
  $BC view -s "$sm" --min-ac 0 /data/cohort/cohort.PASS.vcf.gz > "/tmp/${sm}.raw.vcf"
  awk -v c="$col" -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$c}' \
    results/genotype_matrix.adjudicated.tsv > "/tmp/${sm}.gt.tsv"
  python3 rebuild_vcf_python.py "$sm" "/tmp/${sm}.gt.tsv" "/tmp/${sm}.raw.vcf" "/tmp/${sm}.adj.vcf"
  $BC view "/tmp/${sm}.adj.vcf" -Oz -o "$OUTDIR/${sm}.PASS.adjudicated.vcf.gz"
  $BC index -t "$OUTDIR/${sm}.PASS.adjudicated.vcf.gz"
  echo "[OK] $sm 完成"
done

echo '===== 验收：全部样本 ====='
for f in "$OUTDIR"/*.PASS.adjudicated.vcf.gz; do
  sm=$(basename "$f" .PASS.adjudicated.vcf.gz)
  rec=$($BC view -H "$f" | wc -l)
  nsm=$($BC query -l "$f" | wc -l)
  echo "$sm  records=$rec  samples=$nsm"
done
```

### 5-测序质量 · 代码块 5

```bash
while read -r sm _r1 _r2; do
  echo ">>> CollectHsMetrics $sm"
  singularity exec --bind "$PWD:/data" singularity/gatk_4.6.2.0.sif \
    gatk --java-options -Xmx2g CollectHsMetrics \
    -I "bam/$sm/$sm.markdup.BQSR.bam" \
    -O "qc/hsmetrics/${sm}.hs_metrics.txt" \
    -R /data/reference/genome/genome.fa \
    -BAIT_INTERVALS /data/reference/targets.sorted.interval_list \
    -TARGET_INTERVALS /data/reference/targets.sorted.interval_list
done < samples.tsv
```

### 5-测序质量 · 代码块 6  <!-- 讲解/输出示例块 -->

```bash
for f in qc/hsmetrics/*.hs_metrics.txt; do
  echo -n "$(basename "$f" .hs_metrics.txt): "
  awk -F'\t' '
    /^BAIT_SET/ {for(i=1;i<=NF;i++) h[$i]=i; next}
    $1=="targets" && NF>=60 {
      pq=$(h["PF_UQ_BASES_ALIGNED"])
      printf "on-bait=%.1f%%  on-target=%.1f%%  MEAN=%.0fx  median=%.0fx  >=20x=%.1f%%  FOLD80=%.2f\n",
        $(h["ON_BAIT_BASES"])/pq*100, $(h["ON_TARGET_BASES"])/pq*100,
        $(h["MEAN_TARGET_COVERAGE"]), $(h["MEDIAN_TARGET_COVERAGE"]),
        $(h["PCT_TARGET_BASES_20X"])*100, $(h["FOLD_80_BASE_PENALTY"])
    }' "$f"
done
```

### 5-测序质量 · 代码块 7

```bash
mkdir -p qc/multiqc

singularity exec --bind "$PWD:/data" singularity/multiqc_1.35.sif \
  multiqc /data/qc \
  -o /data/qc/multiqc \
  --title "GWAS Panel 下机数据 QC"
```

---

## 笔记《一些零散知识点》  （共 3 个代码块，其中 1 个含流程命令）


### 一些零散知识点 · 代码块 0  <!-- 讲解/输出示例块 -->

```bash
sudo docker pull broadinstitute/picard:3.4.0
sudo singularity build picard_3.4.0.sif docker-daemon://broadinstitute/picard:3.4.0
```

### 一些零散知识点 · 代码块 1

```bash
singularity exec picard_3.4.0.sif java -Xmx8G -jar /usr/picard/picard.jar FastqToSam -help
```

### 一些零散知识点 · 代码块 2  <!-- 讲解/输出示例块 -->

```perl
#!/usr/bin/env perl
#
# fastq_to_ubam.pl — 批量将 FASTQ 文件通过 Picard FastqToSam 转换为 uBAM
#
# 依赖: Parallel::ForkManager (可通过 cpanm Parallel::ForkManager 安装)
#       若未安装，脚本会自动降级为串行模式
#
# 用法:
#   perl fastq_to_ubam.pl [选项]
#
# 选项:
#   --raw-dir   <dir>    原始 FASTQ 数据根目录            [默认: raw_data]
#   --output    <dir>    输出 uBAM 目录                   [默认: ubam_output]
#   --sif       <file>   Picard Singularity 镜像路径       [默认: picard_3.4.0.sif]
#   --java-mem  <str>    Java 堆内存 (per-job)            [默认: 8G]
#   --threads   <int>    并行任务数                       [默认: 20]
#   --platform  <str>    测序平台                         [默认: illumina]
#   --center    <str>    测序中心                         [默认: BI]
#   --run-date  <str>    RUN_DATE (ISO 8601)              [默认: 自动使用当天]
#   --sample-name <map>  自定义 SAMPLE_NAME 映射          [格式: dir=sample,dir=sample]
#   --dry-run             仅打印命令，不实际执行
#   --verbose             详细输出
#   --help                显示帮助信息
#
# 目录结构要求 (raw_data 下):
#   raw_data/
#   ├── <样本目录>/
#   │   ├── *_R1.fastq.gz
#   │   └── *_R2.fastq.gz
#
# 输出结构:
#   ubam_output/
#   ├── <样本名称>.bam
#   ├── <样本名称>.log
#   └── fastqtosam_<timestamp>.summary.tsv
#

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use POSIX qw(strftime);
use Scalar::Util qw(looks_like_number);
use Fcntl qw(:flock);

# ---------------------------------------------------------------------------
# 命令行参数
# ---------------------------------------------------------------------------
my $raw_data_dir  = "raw_data";
my $output_dir    = "ubam_output";
my $picard_sif    = "picard_3.4.0.sif";
my $java_mem      = "8G";
my $threads       = 20;
my $platform      = "illumina";
my $center        = "BI";
my $run_date      = "";          # 空则自动使用当天日期
my $sample_map_str = "";         # dir1=sample1,dir2=sample2
my $dry_run       = 0;
my $verbose       = 0;
my $show_help     = 0;

GetOptions(
    "raw-dir=s"    => \$raw_data_dir,
    "output=s"     => \$output_dir,
    "sif=s"        => \$picard_sif,
    "java-mem=s"   => \$java_mem,
    "threads=i"    => \$threads,
    "platform=s"   => \$platform,
    "center=s"     => \$center,
    "run-date=s"   => \$run_date,
    "sample-name=s"=> \$sample_map_str,
    "dry-run"      => \$dry_run,
    "verbose"      => \$verbose,
    "help"         => \$show_help,
) or die "选项解析失败，使用 --help 查看帮助\n";

if ($show_help) {
    print_help();
    exit 0;
}

# ---------------------------------------------------------------------------
# 解析用户自定义的 SAMPLE_NAME 映射 (dir1=sample1,dir2=sample2)
# ---------------------------------------------------------------------------
my %custom_sample_names;
if ($sample_map_str) {
    my @pairs = split /,/, $sample_map_str;
    for my $pair (@pairs) {
        my ($dir, $name) = split /=/, $pair, 2;
        $custom_sample_names{$dir} = $name if defined $name;
    }
}

# ---------------------------------------------------------------------------
# RUN_DATE 默认值
# ---------------------------------------------------------------------------
$run_date = strftime("%Y-%m-%dT%H:%M:%S%z", localtime) unless $run_date;

# ---------------------------------------------------------------------------
# 验证输入
# ---------------------------------------------------------------------------
die "错误: raw-data 目录不存在: $raw_data_dir\n" unless -d $raw_data_dir;
$raw_data_dir = abs_path($raw_data_dir);

make_path($output_dir) unless $dry_run;
$output_dir = abs_path($output_dir);

die "错误: Picard SIF 文件不存在: $picard_sif\n"
    unless $dry_run || -f $picard_sif;

# ---------------------------------------------------------------------------
# 发现样本 — 扫描 raw_data 下所有包含 FASTQ 的子目录
# ---------------------------------------------------------------------------
sub discover_samples {
    my ($raw_dir) = @_;
    my %samples;

    opendir(my $dh, $raw_dir) or die "无法打开目录 $raw_dir: $!\n";
    my @entries = sort readdir($dh);
    closedir($dh);

    for my $entry (@entries) {
        next if $entry =~ /^\./;
        my $subdir = File::Spec->catdir($raw_dir, $entry);
        next unless -d $subdir;

        # 搜索 R1/R2 FASTQ 文件
        my ($r1, $r2);
        opendir(my $sdh, $subdir) or do {
            warn "警告: 无法读取 $subdir: $!\n";
            next;
        };
        my @files = sort readdir($sdh);
        closedir($sdh);

        for my $f (@files) {
            next if $f =~ /^\./;
            if ($f =~ /_R1\.fastq\.gz$/i || $f =~ /_R1\.fq\.gz$/i || $f =~ /_1\.fastq\.gz$/i || $f =~ /_1\.fq\.gz$/i) {
                $r1 = File::Spec->catfile($subdir, $f);
            } elsif ($f =~ /_R2\.fastq\.gz$/i || $f =~ /_R2\.fq\.gz$/i || $f =~ /_2\.fastq\.gz$/i || $f =~ /_2\.fq\.gz$/i) {
                $r2 = File::Spec->catfile($subdir, $f);
            }
        }

        if ($r1 && $r2) {
            $samples{$entry} = {
                dirname  => $entry,
                r1       => $r1,
                r2       => $r2,
            };
        } elsif ($r1) {
            warn "警告: $entry 仅有 R1，跳过（需要双端数据）\n";
        } elsif ($r2) {
            warn "警告: $entry 仅有 R2，跳过（需要双端数据）\n";
        } else {
            warn "警告: $entry 中未找到 FASTQ 文件，跳过\n";
        }
    }

    return %samples;
}

# ---------------------------------------------------------------------------
# 从目录名推导读组字段
#
# 规则:
#   SAMPLE_NAME    — 优先使用用户自定义映射，否则尝试提取前缀
#                    CB261007720-OSW → CB261007720 (SAL/OSW 后缀剥离)
#                    NA12878         → NA12878
#   READ_GROUP_NAME — 目录名本身，如 CB261007720-OSW
#   LIBRARY_NAME    — 目录名 + .lib
#   PLATFORM_UNIT   — 目录名 + .unit (理想情况下应为 run_barcode.lane)
# ---------------------------------------------------------------------------
sub derive_read_group {
    my ($info) = @_;
    my $dirname = $info->{dirname};

    # SAMPLE_NAME
    my $sample_name;
    if (exists $custom_sample_names{$dirname}) {
        $sample_name = $custom_sample_names{$dirname};
    } elsif ($dirname =~ /^(.+)-(?:OSW|SAL)$/i) {
        # CB261007720-OSW → CB261007720
        $sample_name = $1;
    } else {
        $sample_name = $dirname;
    }

    # READ_GROUP_NAME — 使用目录名
    my $read_group_name = $dirname;

    # LIBRARY_NAME — 简单派生
    my $library_name = $dirname . ".lib";

    # PLATFORM_UNIT — 理想格式: run_barcode.lane，此处用目录名做占位
    my $platform_unit = $dirname . ".unit";

    return {
        sample_name    => $sample_name,
        read_group_name=> $read_group_name,
        library_name   => $library_name,
        platform_unit  => $platform_unit,
    };
}

# ---------------------------------------------------------------------------
# 构建单个 FastqToSam 命令
# ---------------------------------------------------------------------------
sub build_command {
    my ($info, $rg, $output_bam) = @_;

    my $fastq1 = $info->{r1};
    my $fastq2 = $info->{r2};

    # 构建 Java 命令
    my @cmd = (
        "singularity", "exec", $picard_sif,
        "java", "-Xmx$java_mem",
        "-jar", "/usr/picard/picard.jar",
        "FastqToSam",
        "FASTQ=$fastq1",
        "FASTQ2=$fastq2",
        "OUTPUT=$output_bam",
        "READ_GROUP_NAME=$rg->{read_group_name}",
        "SAMPLE_NAME=$rg->{sample_name}",
        "LIBRARY_NAME=$rg->{library_name}",
        "PLATFORM_UNIT=$rg->{platform_unit}",
        "PLATFORM=$platform",
        "SEQUENCING_CENTER=$center",
        "RUN_DATE=$run_date",
        "SORT_ORDER=queryname",
    );

    return @cmd;
}

# ---------------------------------------------------------------------------
# 打印帮助
# ---------------------------------------------------------------------------
sub print_help {
    print <<"HELP";
用法: perl fastq_to_ubam.pl [选项]

选项:
  --raw-dir   <dir>     原始 FASTQ 数据根目录                [默认: raw_data]
  --output    <dir>     输出 uBAM 目录                       [默认: ubam_output]
  --sif       <file>    Picard Singularity 镜像路径          [默认: picard_3.4.0.sif]
  --java-mem  <str>     Java 堆内存，每个任务独立分配        [默认: 8G]
  --threads   <int>     并行任务数                           [默认: 20]
  --platform  <str>     测序平台                             [默认: illumina]
  --center    <str>     测序中心                             [默认: BI]
  --run-date  <str>     RUN_DATE (ISO 8601 格式)             [默认: 当天日期]
  --sample-name <map>   自定义 SAMPLE_NAME 映射
                        格式: dir1=s1,dir2=s2
                        例: --sample-name "CB261007720-OSW=SAMPLE1,NA12878=REF1"
  --dry-run             仅打印命令，不实际执行
  --verbose             详细输出
  --help                显示该帮助

目录结构:
  raw_data/
  ├── <样本目录>/
  │   ├── *_R1.fastq.gz
  │   └── *_R2.fastq.gz

输出:
  ubam_output/
  ├── <read_group_name>.bam      — 输出的 uBAM 文件
  ├── <read_group_name>.log      — stderr/stdout 日志
  └── fastqtosam_<时间戳>.tsv   — 运行汇总

HELP
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
sub main {
    print "========================================\n";
    print "  FastqToSam 批量转换 - 运行配置\n";
    print "========================================\n";
    print "  原始数据目录:  $raw_data_dir\n";
    print "  输出目录:      $output_dir\n";
    print "  Picard SIF:    $picard_sif\n";
    print "  Java 堆内存:   $java_mem (per-job)\n";
    print "  并行任务数:    $threads\n";
    print "  测序平台:      $platform\n";
    print "  测序中心:      $center\n";
    print "  RUN_DATE:      $run_date\n";
    print "  模式:          " . ($dry_run ? "干运行 (不执行命令)" : "正式运行") . "\n";
    if (%custom_sample_names) {
        print "  自定义 SAMPLE_NAME:\n";
        for (keys %custom_sample_names) {
            print "    $_ => $custom_sample_names{$_}\n";
        }
    }
    print "========================================\n\n";

    # 发现样本
    my %samples = discover_samples($raw_data_dir);

    if (!%samples) {
        print "未发现任何样本，退出。\n";
        exit 0;
    }

    print "发现 " . scalar(keys %samples) . " 个样本:\n";
    for my $key (sort keys %samples) {
        my $rg = derive_read_group($samples{$key});
        printf "  %-30s  SAMPLE_NAME=%s  RG=%s\n",
            $key, $rg->{sample_name}, $rg->{read_group_name};
    }
    print "\n";

    # 构建任务列表
    my @tasks;
    for my $key (sort keys %samples) {
        my $info = $samples{$key};
        my $rg   = derive_read_group($info);

        my $output_bam = File::Spec->catfile($output_dir, $rg->{read_group_name} . ".bam");
        my $log_file   = File::Spec->catfile($output_dir, $rg->{read_group_name} . ".log");
        my @cmd        = build_command($info, $rg, $output_bam);

        push @tasks, {
            info       => $info,
            rg         => $rg,
            command    => \@cmd,
            output_bam => $output_bam,
            log_file   => $log_file,
        };
    }

    # 干运行模式：仅打印命令
    if ($dry_run) {
        my $idx = 0;
        for my $task (@tasks) {
            $idx++;
            print "--- 任务 $idx: $task->{rg}{read_group_name} ---\n";
            print "  " . join(" ", @{$task->{command}}) . "\n";
            print "  日志: $task->{log_file}\n\n";
        }
        print "干运行完成。以上命令未实际执行。\n";
        exit 0;
    }

    # ---------------------------------------------------------------------------
    # 并行执行 (Parallel::ForkManager)
    # ---------------------------------------------------------------------------
    my $has_pfm = eval {
        require Parallel::ForkManager;
        Parallel::ForkManager->import;
        1;
    };

    if ($has_pfm) {
        print "使用 Parallel::ForkManager 并行执行 (max workers = $threads)\n\n";
    } else {
        print "Parallel::ForkManager 未安装，降级为串行执行。\n";
        print "安装方法: cpanm Parallel::ForkManager\n\n";
    }

    my $pm;
    if ($has_pfm) {
        $pm = Parallel::ForkManager->new($threads);
    }

    my $total     = scalar @tasks;
    my $completed = 0;
    my $failed    = 0;
    my @failed_tasks;

    my $start_time = time;

    # 汇总文件
    my $timestamp  = strftime("%Y%m%d_%H%M%S", localtime);
    my $summary_tsv = File::Spec->catfile($output_dir, "fastqtosam_${timestamp}.tsv");
    open(my $summary_fh, ">", $summary_tsv) or die "无法写入汇总文件: $!\n";
    print $summary_fh join("\t", qw(
        ReadGroupName SampleName LibraryName PlatformUnit
        OutputBam Status ExitCode DurationSec LogFile
    )) . "\n";

    for my $task (@tasks) {
        my $rg   = $task->{rg};
        my @cmd  = @{$task->{command}};

        if ($has_pfm) {
            # 子进程开始
            $pm->start and next;

            my $task_start = time;
            my ($exit_code, $status);

            # 重定向输出到日志文件
            my $cmd_str = join(" ", map { /[\s\*?{}]/ ? "'$_'" : $_ } @cmd);
            system("$cmd_str > '$task->{log_file}' 2>&1");
            $exit_code = $?;

            my $duration = time - $task_start;

            if ($exit_code == 0) {
                $status = "SUCCESS";
            } else {
                $status = "FAILED";
            }

            # 每个子进程独立写汇总行
            my $summary_line = join("\t",
                $rg->{read_group_name},
                $rg->{sample_name},
                $rg->{library_name},
                $rg->{platform_unit},
                $task->{output_bam},
                $status,
                ($exit_code == 0 ? 0 : ($exit_code >> 8)),
                $duration,
                $task->{log_file},
            ) . "\n";

            # 加文件锁写入
            flock($summary_fh, LOCK_EX) or warn "无法加锁汇总文件: $!\n";
            seek($summary_fh, 0, 2);  # 追加到末尾
            print $summary_fh $summary_line;
            flock($summary_fh, LOCK_UN);

            print "[$status] $rg->{read_group_name}  (${duration}s)\n";

            $pm->finish($exit_code);
        } else {
            # ---- 串行执行 ----
            my $task_start = time;

            print "执行: $rg->{read_group_name} ...\n";
            if ($verbose) {
                print "  命令: " . join(" ", @cmd) . "\n";
            }

            my $cmd_str = join(" ", map { /[\s\*?{}]/ ? "'$_'" : $_ } @cmd);
            system("$cmd_str > '$task->{log_file}' 2>&1");
            my $exit_code = $?;
            my $duration  = time - $task_start;

            my $status;
            if ($exit_code == 0) {
                $completed++;
                $status = "SUCCESS";
                print "  [OK]  ${duration}s\n";
            } else {
                $failed++;
                $status = "FAILED";
                my $rc = $exit_code >> 8;
                print "  [FAIL] 退出码=$rc  ${duration}s  日志: $task->{log_file}\n";
                push @failed_tasks, $task;
            }

            print $summary_fh join("\t",
                $rg->{read_group_name},
                $rg->{sample_name},
                $rg->{library_name},
                $rg->{platform_unit},
                $task->{output_bam},
                $status,
                ($exit_code == 0 ? 0 : ($exit_code >> 8)),
                $duration,
                $task->{log_file},
            ) . "\n";
        }
    }

    if ($has_pfm) {
        # 回收所有子进程并统计
        $pm->wait_all_children;

        # 从汇总文件中重新统计
        seek($summary_fh, 0, 0);
        while (<$summary_fh>) {
            chomp;
            next if /^ReadGroupName/;  # 跳过标题行
            my @f = split /\t/;
            if ($f[5] && $f[5] eq "SUCCESS") {
                $completed++;
            } elsif ($f[5] && $f[5] eq "FAILED") {
                $failed++;
            }
        }
    }

    close($summary_fh);

    my $total_duration = time - $start_time;
    my $wall_clock = sprintf("%d:%02d:%02d",
        int($total_duration / 3600),
        int(($total_duration % 3600) / 60),
        $total_duration % 60);

    print "\n========================================\n";
    print "  运行汇总\n";
    print "========================================\n";
    print "  总任务数:    $total\n";
    print "  成功:        $completed\n";
    print "  失败:        $failed\n";
    print "  总耗时:      $wall_clock\n";
    print "  汇总文件:    $summary_tsv\n";
    print "========================================\n";

    if ($failed > 0) {
        print "\n有 $failed 个任务失败，请检查对应日志文件。\n";
        exit 1;
    }
}

main();
__END__
```

---

## 附录：正文散文中出现的独立命令行（不在代码块内，供交叉核对）

