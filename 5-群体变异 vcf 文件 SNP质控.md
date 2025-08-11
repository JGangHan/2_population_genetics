- **经 vcftools 拆分染色体，与 gatk SelectVariants 拆分 SNP 和 Indel 后的 vcf 文件**

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

- **GATK：过滤掉技术上“不可信”的变异位点**
- **vcftools：过滤掉群体学分析中“不合格”的位点/个体**
- **PLINK：过滤掉统计分析中“不合理”的样本和位点**


## 2. GATK 质检 + bcftools 质控
- **gatk 质检特征**
  - GATK VariantFiltration 针对 VCF 文件的 INFO 字段(QD、MQ、FS、SOR、MQRankSum、ReadPosRankSum)进行质检
  - GATK VariantFiltration 不能直接删除位点，而是添加是否质量合格标签
  - 在 VCF 文件中给不满足条件的 SNP 加上 Filter 标签
  - 后续可用 bcftools 或 gatk SelectVariants 删除不合格位点
  - 硬过滤  
**逐条染色体 vcf 文件质控**
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

## 3. vcftools 质控特征
- 样本/群体级过滤（sample/population-level filter），针对 VCF 的基因型信息（FORMAT + GT 字段） 进行质控
- 包括：过滤低 MAF 变异（去掉低频等位基因）、限制缺失率（max-missing）、保留双等位 SNP（min-alleles / max-alleles）、控制平均测序深度（min-meanDP）。
- 关注变异在整个群体中的表现，确保最终保留的位点在统计分析中可用
- 会直接生成新的 VCF（不只是打标签），不合格位点直接被剔除。
- 硬过滤  
**继续逐条染色体 vcf 文件质控**
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

    vcftools --gzvcf "${VCF_DIR}/chr${CHROM}.filtered.pass.vcf.gz" \
      --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.9 --min-meanDP 3 \
      --recode --recode-INFO-all --out "${OUTPUT_DIR}/chr${CHROM}.filtered.snp"

    # 修复缺失的基因型格式，在 VCF 中，样本的基因型缺失值应写成 ./. 而不是单独一个 .
    # 将  “ .:” 替换为 "  ./.:"
    perl -pe 's/\s\.:/\t.\/.:/g' "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf" > "${OUTPUT_DIR}/chr${CHROM}.snp.fix.recode.vcf"

    # 删除中间文件
    rm "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf"
    rm "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf"

    echo "chr${CHROM} done."
  ) &
done
wait
```

## 4. 可能遇见的报错
**gatk 过滤后文件损毁**
- 内容：vcftools 对 chr9 和 chr15.filtered.vcf.gz 进行质控过程命令卡住（注意不是报错或停止，染色体随机）
- **检查方式：gunzip -c 解压检查，若 vcf.gz 文件损毁，则直接报错**
- 结果：gatk snp 质检步骤后，chr9 和 chr15 染色体 vcf 文件损毁（或不完善），用 gatk 重新过滤即可














