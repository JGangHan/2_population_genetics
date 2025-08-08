预处理流程：
- vcftools 拆分染色体
- gatk SelectVariants 拆分 SNP 和 Indel

## 1. SNP 质控信息
**（1）常见 SNP 质控软件、参数和阈值**
| 参数                                   | 中文释义                | 推荐默认值         | 可调范围      | 主要影响                     | 常见适用场景                        | 推荐质控工具                                          |
| ------------------------------------ | ------------------- | ------------- | --------- | ------------------------ | ----------------------------- | ----------------------------------------------- |
| **QD (Quality by Depth)**            | 质量深度比               | `< 2.0` 剔除    | 1.5–3.0   | 排除低质量/低深度支持的变异           | 默认 2.0 适合大多数群体重测序；低深度可降至 1.5  | **GATK VariantFiltration**                      |
| **MQ (Mapping Quality)**             | 比对质量                | `< 40.0` 剔除   | 35–60     | 反映测序比对的可信度               | 高复杂度基因组可调高至 50–60，低复杂度可放宽至 35 | **GATK VariantFiltration**                      |
| **FS (Fisher Strand Bias)**          | 链偏性检验（Fisher 检验）    | `> 60.0` 剔除   | 40–80     | 检测正反链测序支持比例是否异常          | 默认 60；高置信度需求可用 40             | **GATK VariantFiltration**                      |
| **SOR (Symmetric Odds Ratio)**       | 对称比值比（链偏性补充指标）      | `> 3.0` 剔除    | 2.5–4.0   | 辅助判断链偏性                  | 高深度数据可降至 2.5 过滤更严格            | **GATK VariantFiltration**                      |
| **MQRankSum**                        | 比对质量秩和检验            | `< -12.5` 剔除  | -15 – -8  | 比较 REF 与 ALT 等位基因的比对质量差异 | 默认 -12.5；高质量数据可用 -10 提高灵敏度    | **GATK VariantFiltration**                      |
| **ReadPosRankSum**                   | 读长位置秩和检验            | `< -8.0` 剔除   | -10 – -6  | 判断变异位点是否偏向读长两端           | 默认 -8 平衡灵敏度与特异性               | **GATK VariantFiltration**                      |
| **MAF (Minor Allele Frequency)**     | 最小等位基因频率            | `≥ 0.05` 保留   | 0.01–0.10 | 过滤低频变异                   | 群体结构分析常用 0.05；GWAS 可用 0.01    | **VCFtools** 或 **PLINK**                        |
| **等位基因数 (biallelic)**                | 双等位要求               | 2（双等位）        | 固定        | 仅保留纯 SNP                 | 适用于只分析 SNP 的项目                | **VCFtools**（`--min-alleles 2 --max-alleles 2`） |
| **缺失率 (max-missing)**                | 样本基因型完整率            | ≥ 0.90        | 0.70–0.95 | 控制缺失数据比例                 | 群体分析建议 ≥0.90；低深度可放宽至 0.80     | **VCFtools**（`--max-missing`）或 **PLINK**        |
| **平均测序深度 (min-meanDP)**              | 每个位点平均测序深度          | ≥ 3           | 2–8       | 剔除低覆盖位点                  | 高深度测序可调至 8–10；低深度可用 2–3       | **VCFtools**（`--min-meanDP`）                    |
| **HWE (Hardy-Weinberg Equilibrium)** | Hardy–Weinberg 平衡检验 | `p > 1e-6` 保留 | 1e-6–1e-4 | 检测群体遗传异常                 | GWAS 或纯种群体建议加此步；杂交群体可不使用      | **VCFtools**（`--hwe`）或 **PLINK**                |

**（2）常用 SNP 质控组合策略**
- GATK → 位点质量控制（上游）
  - 针对测序/比对质量层面的过滤，比如：QD、MQ、FS、SOR、MQRankSum、ReadPosRankSum
  - 依赖 BAM 中的比对信息，只有在变异检测阶段（GATK/DeepVariant 等）才能计算
  - 剔除测序或比对产生的低可信变异，确保留下的位点都是技术上可靠的
- vcftools → 基于群体统计过滤（中游）
  - 针对群体水平 SNP 位点统计学结果进行过滤，比如：MAF（小等位基因频率）、max-missing（缺失率）、min-meanDP（平均深度）、biallelic（双等位）
  - 去掉在群体中罕见或缺失太多的位点，以及低深度等潜在假阳性位点
- PLINK → 下游补充质控
  - 常用于 GWAS、PCA、结构分析前的最终 QC：HWE（Hardy–Weinberg 平衡）、亲缘关系过滤（IBD、PI_HAT）、样本级缺失率过滤、LD 剔除（--indep-pairwise）
  - 确保群体学分析符合遗传学假设，避免因近交、群体分化等造成假阳性
**
- GATK：过滤掉技术上“不可信”的变异位点
- vcftools：过滤掉群体学分析中“不合格”的位点/个体
- PLINK：过滤掉统计分析中“不合理”的样本和位点


## 2. GATK 质检 + bcftools 质控
**gatk 质控特征**
- GATK VariantFiltration 针对 VCF 文件的 INFO 字段(QD、MQ、FS、SOR、MQRankSum、ReadPosRankSum)进行质检
- GATK VariantFiltration 不能直接删除位点，而是添加是否质量合格标签
- 在 VCF 文件中给不满足条件的 SNP 加上 FILTER 标签（Filter）
- 后续可用 bcftools 或 gatk SelectVariants 删除不合格位点
- 硬过滤





## 3. vcftools 质控特征
- 样本/群体级过滤（sample/population-level filter），针对 VCF 的基因型信息（FORMAT + GT 字段） 进行质控
- 包括：过滤低 MAF 变异（去掉低频等位基因）、限制缺失率（max-missing）、保留双等位 SNP（min-alleles / max-alleles）、控制平均测序深度（min-meanDP）。
- 关注变异在整个群体中的表现，确保最终保留的位点在统计分析中可用
- 会直接生成新的 VCF（不只是打标签），不合格位点直接被剔除。
- 硬过滤


## 2. 
SNP 质控（GATK 和 vcftools）
**常用质控软件有 gatk、vcftools、plink等**  
**首先用 GATK过滤**
```
# 1. 加载gatk
module load /workspace/public/x86/software/modules/tool/GATK-4.0.0.0


# 2. 指定变量
VCF_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP"
OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter"
mkdir -p $OUTPUT_DIR


# 3. gatk 过滤
# --java-options 指定线程
# --filter-expression  多项过滤指标阈值（标准）
# --filter-name "Filter" 质量不合格变异添加 Filter 标签
# 对 26 条染色体逐条过滤
for CHROM in {1..26}; do
  (
    gatk --java-options "-Xmx4G -Djava.io.tmpdir=./" VariantFiltration \
      -V "${VCF_DIR}/chr${CHROM}.recode_snp.vcf.gz" \
      -O "${OUTPUT_DIR}/chr${CHROM}.filtered.vcf.gz" \
      --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
      --filter-name "Filter" &
  ) &
done
wait


for CHROM in {1..26}; do
  (
    # 打标签
    gatk --java-options "-Xmx4G -Djava.io.tmpdir=./" VariantFiltration \
      -V "${VCF_DIR}/chr${CHROM}.recode_snp.vcf.gz" \
      -O "${OUTPUT_DIR}/chr${CHROM}.filtered.vcf.gz" \
      --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
      --filter-name "Filter"

    # 只保留过滤通过的SNP（FILTER字段为PASS）
    bcftools view -f PASS "${OUTPUT_DIR}/chr${CHROM}.filtered.vcf.gz" -Oz -o "${OUTPUT_DIR}/chr${CHROM}.filtered.pass.vcf.gz"

    # 重新索引
    tabix -p vcf "${OUTPUT_DIR}/chr${CHROM}.filtered.pass.vcf.gz"
  ) &
done
wait


```
  
**再用 vcftools 过滤**
```
# 1. 加载模块
module load /workspace/public/x86/software/modules/tool/vcftools-0.1.16 


# 2. 指定变量和文件路径
VCF_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter"
OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter"
mkdir -p $OUTPUT_DIR


# 3. 26 条染色体进行并行过滤
# --gzvcf：输入文件（bgzip 压缩的 VCF 文件）；
# --maf 0.05：去除小于5%最小等位基因频率（MAF）的 SNP；
# --min-alleles 2 / --max-alleles 2：只保留二等位变异（即真正的 SNP，去除多等位或复杂变异）；
# --max-missing 0.9：最少90%的样本有基因型（允许10%缺失）；
# --min-meanDP 3：平均测序深度小于3的位点将被剔除；
# --recode：输出新的 VCF；
# --recode-INFO-all：保留所有 INFO 字段；
# --out：输出文件前缀（生成 .recode.vcf 文件）；

for CHROM in {1..26}; do
  (
    echo "Processing chr${CHROM}..."

    vcftools --gzvcf "${VCF_DIR}/chr${CHROM}.filtered.vcf.gz" \
      --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.9 --min-meanDP 3 \
      --recode --recode-INFO-all --out "${OUTPUT_DIR}/chr${CHROM}.filtered.snp"

    # 删除特定 contig 信息
    grep -v "##contig=<ID=LWLT*" "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf" > "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf"

    # 修复缺失的基因型格式
    perl -pe 's/\s\.:/\t.\/.:/g' "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf" > "${OUTPUT_DIR}/chr${CHROM}.snp.fix.recode.vcf"

    # 删除中间文件
    rm "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf"
    rm "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf"

    echo "chr${CHROM} done."
  ) &
done
wait
```

### 6. 质控后处理
- 合并各染色体 .vcf.gz 文件
- 染色体编号重命名
- 各 vcf 文件一致性检查
```
#!/bin/bash

# 理论上每条染色体 .vcf.gz 文件的头部信息一致，所以仅需提取一个染色体的头部信息，之后将每条染色体的数据信息部分逐个添加

# 输出文件路径
OUT_VCF="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf"

# 获取 header 并添加
grep "^#" /workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr1.snp.fix.recode.vcf > "$OUT_VCF"

# 合并各染色体 body
# grep -v 反向检索
for CHR in {1..26}; do
  FILE="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr${CHR}.snp.fix.recode.vcf"
  echo "Merging: $FILE"
  grep -v "^#" "$FILE" >> "$OUT_VCF"
done
```

#### 6.2 查看 vcf.gz 文件相关信息
1. 个体名称和个体数量
```
# 直接输出样品名称
bcftools query -l qmrb_50sheep_snp.vcf
# 样品数量
bcftools query -l qmrb_50sheep_snp.vcf | wc -l
```
2. 头部信息染色体编号
```
grep "^##contig=" ./qmrb_50sheep_snp.vcf
```
3. 数据部分染色体编号
**最好在这一步就提前检查后续涉及的 vcf 文件染色体名称和长度是否一致，如果发现问题最好提前解决**
```
# cut -f1  提取第一个字段

# 拆分质控后的染色体文件
cut -f1 /workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr26.snp.fix.recode.vcf | grep -v "^#" | sort | uniq
# 合并文件
cut -f1 /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf | grep -v "^#" | sort | uniq
# 第三方处理文件
cut -f1 /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf | grep -v "^#" | sort | uniq
```


#### 6.3 vcf 文件染色体重命名
将染色体编号命名为标准格式：最好是直接用数字代表染色体（1，2，3...），而不是字符串表示染色体（chr1，chr2，chr3......）  
**感觉染色体重命名应该在最先合并的vcf文件中处理，但刚合并vcf文件中除了常染色体，还有很多其他的 contig，可能不太好处理**  
**需要注意：因为染色体名称同时出现在头部信息和数据部分，为了严谨，最好一起修改**
1. 染色体编号重命名
```
# 染色体重命名
input_file="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf"
output_file="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep.vcf"

# CM001582.2    chr1
# ........
# CM001607.2    chr26

# 使用sed命令替换染色体名称
# sed -e：表示后面跟的是一个 sed 的表达式（script）。
# 's/CM001582.2/1/g'：s/旧内容/新内容/g
# s：表示“substitute”（替换）    /CM001582.2/：要被替换的原始字符串     /1/：替换成的内容      g：表示“全局替换”，即一行中出现多次也全部替换。
sed -e 's/CM001582.2/1/g' \
    -e 's/CM001583.2/2/g' \
    -e 's/CM001584.2/3/g' \
    -e 's/CM001585.2/4/g' \
    -e 's/CM001586.2/5/g' \
    -e 's/CM001587.2/6/g' \
    -e 's/CM001588.2/7/g' \
    -e 's/CM001589.2/8/g' \
    -e 's/CM001590.2/9/g' \
    -e 's/CM001591.2/10/g' \
    -e 's/CM001592.2/11/g' \
    -e 's/CM001593.2/12/g' \
    -e 's/CM001594.2/13/g' \
    -e 's/CM001595.2/14/g' \
    -e 's/CM001596.2/15/g' \
    -e 's/CM001597.2/16/g' \
    -e 's/CM001598.2/17/g' \
    -e 's/CM001599.2/18/g' \
    -e 's/CM001600.2/19/g' \
    -e 's/CM001601.2/20/g' \
    -e 's/CM001602.2/21/g' \
    -e 's/CM001603.2/22/g' \
    -e 's/CM001604.2/23/g' \
    -e 's/CM001605.2/24/g' \
    -e 's/CM001606.2/25/g' \
    -e 's/CM001607.2/26/g' \
    $input_file > $output_file
```

2. 头部信息中仅保留常染色体
```
#!/bin/bash
input_vcf="qmrb_50sheep.vcf.gz"
output_vcf="qmrb_50sheep_chr1to26.vcf.gz"
# 常染色体 1-26
zcat "$input_vcf" | \
awk '
BEGIN {
    for (i=1; i<=26; i++) keep[i]=1
}
{
    if ($0 ~ /^##contig=<ID=/) {
        match($0, /ID=([^,>]+)/, arr)
        if (arr[1] in keep)
            print
    } else if ($0 ~ /^##/) {
        print
    } else if ($0 ~ /^#CHROM/) {
        print
        in_data=1
    } else if (in_data) {
        print
    }
}' | bgzip -c > "$output_vcf"

# 索引
tabix -p vcf "$output_vcf"

grep "^##contig=" ./qmrb_50sheep_snp.vcf | head   # 重命名前染色体
##contig=<ID=CM001582.2,length=275406953,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001583.2,length=248966461,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001584.2,length=223996068,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001585.2,length=119216639,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001586.2,length=107836144,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001587.2,length=116888256,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001588.2,length=100009711,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001589.2,length=90615088,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001590.2,length=94583238,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=CM001591.2,length=86377204,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>

grep "^##contig=" ./qmrb_50sheep.vcf | head   # 重命名后染色体
##contig=<ID=1,length=275406953,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=2,length=248966461,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=3,length=223996068,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=4,length=119216639,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=5,length=107836144,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=6,length=116888256,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=7,length=100009711,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=8,length=90615088,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=9,length=94583238,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=10,length=86377204,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>

bcftools view -h ./qmrb_50sheep_chr1to26.vcf.gz | grep "^##contig=" | tail  # 头部信息末尾多余染色体信息被移除
##contig=<ID=17,length=72251135,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=18,length=68494538,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=19,length=60445663,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=20,length=51049468,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=21,length=49987992,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=22,length=50780147,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=23,length=62282865,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=24,length=41976827,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=25,length=45223504,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
##contig=<ID=26,length=44047080,assembly=GCA_000298735.2_Oar_v4.0_genomic.fna>
```
### 7. 遇见问题
#### 7.1 gatk 过滤后文件损毁
1. vcftools 对 chr9 和 chr15.filtered.vcf.gz 进行质控过程命令卡住（注意不是报错或停止）  
若 vcf.gz 文件损毁，通过 gunzip -c 解压会直接报错
3. 经逐项检查，发现是 gatk 过滤步骤后，chr9 和 chr15 染色体 vcf 文件损毁（或不完善），用 gatk 重新过滤即可
#### 7.2 














