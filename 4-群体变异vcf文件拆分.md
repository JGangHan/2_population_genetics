**流程：按染色体拆分 → 提取变异类型（SNP 和 Indel） → 逐染色体质控**  
**Why？**  
1. 按染色体拆分-处理效率更高，利于并行化计算
- 群体变异 VCF 文件体积大，若同时处理，耗时长，占用内存大，容易报错
2. 按染色体拆分-有助于后续分析
- 下游分析，如GWAS、群体结构分析、连锁不平衡LD计算等，以染色体为单位。
- 后续可快速定位目标染色体变异、避免每次都读入冗余信息。
3. 按变异类型拆分-SNP 和 Indel 的质量控制标准是不同的：
- SNP	QD < 2.0、FS > 60.0、SOR > 3.0、MQ < 40.0 等
- Indel	通常使用 FS > 200.0、SOR > 10.0 等更严格的阈值
- 混合质控使用统一参数：对 SNP 过松、对 Indel 过严，后续分析难度增大


## 1. 群体 vcf.gz 文件 → 按染色体拆分
**使用 vcftools**
```
# 根据物种染色体数量修改：绵羊26对常染色体，1对性染色体**

# 1. 加载软件（视具体情况而定）
module load /workspace/public/x86/software/modules/tool/vcftools-0.1.16 


# 2. 检查是否包含 contig 片段
# （1）尽量在参考基因组比对之前去除常染色体、性染色体和线粒体之外的contig片段
# （2）NCBI 参考基因组还要检查染色体编码（CM001582.2）和染色体号（chr1）对应关系
# 编码          编号
# CM001582.2    chr1
# CM001607.2    chr26
# CM001608.2    chrX

zcat /workspace/public/zongzhan/qmrb_sheep/jointcall/qmrb_50sheep.vcf.gz | awk '/^#/ {print} !/^#/ {exit}'
# 输出结果为染色体编码，而非编号


# 3. 变量
# 输入路径
INPUT_VCF="/workspace/public/zongzhan/qmrb_sheep/jointcall/qmrb_50sheep.vcf.gz"
# 输出路径
OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split"
mkdir -p ${OUTPUT_DIR}


# 4. 按染色体拆分，并针对不同染色体 vcf.gz 文件重新命名
# CM001582.2    chr1
# .........
# CM001607.2    chr26
# CM001608.2    chrX

# vcftools
# --gzvcf 输入文件为 .vcf.gz 格式
# --chr 提取指定染色体
# --recode 输出为新的 vcf 文件
# --recode-INFO-all：保留所有 INFO 字段

# bgzip VCF 文件压缩为 .vcf.gz 格式 -f：若存在则强制覆盖
# tabix 索引压缩后 .vcf.gz 文件；-p vcf 指定格式类型为 VCF

# 常染色体 CM001582.2 至 CM001607.2
for i in $(seq 582 607)
do
    OLD_CHR="CM001${i}.2"
    NEW_NUM=$((i - 581))  # 将染色体编号转换为1到26的范围
    (
        vcftools --gzvcf ${INPUT_VCF} --chr ${OLD_CHR} --recode --recode-INFO-all --out ${OUTPUT_DIR}/chr${NEW_NUM} &&
        bgzip -f ${OUTPUT_DIR}/chr${NEW_NUM}.recode.vcf &&
        tabix -p vcf ${OUTPUT_DIR}/chr${NEW_NUM}.recode.vcf.gz
    ) & 
done
# 等待所有后台任务完成
wait


# 5. 重新检查一遍
# 查看头部信息
bcftools view -h ./jointcall/split/chr2.recode.vcf.gz | head
# 查看数据部分
bcftools view -H ./jointcall/split/chr2.recode.vcf.gz | head
```



## 2. 单个染色体 vcf.gz 文件 → 按变异类型拆分
**使用 gatk SelectVariants 拆分SNP 和 Indel 两种变异类型**
```
# 1. 加载软件
module load /workspace/public/x86/software/modules/tool/GATK-4.0.0.0


# 2. 变量
INPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split"
SNP_OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP"
INDEL_OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/indel"
# 注意这里的参考基因组和比对用参考基因组必须完全一致
REFERENCE="/workspace/public/zongzhan/Reference/Sheep_Oar_4.0/ncbi_dataset/data/GCA_000298735.2/GCA_000298735.2_Oar_v4.0_genomic.fna"

# 创建输出目录（如果不存在）
mkdir -p ${SNP_OUTPUT_DIR}
mkdir -p ${INDEL_OUTPUT_DIR}


# 3. gatk SelectVariants
# 输入目录下所有 .vcf 文件列表
VCF_FILES=(${INPUT_DIR}/*.vcf.gz)

# 遍历每个VCF文件，从中分别拆分出 SNP 和 Indel 两种变异类型
for vcf_file in "${VCF_FILES[@]}"; do
    base_name=$(basename "${vcf_file}" .vcf.gz)  # basename 命令可去除去除路径与扩展名
    
    # GATK SelectVariants 对单个染色体 vcf 文件拆分 SNP
    gatk SelectVariants \
        -R "${REFERENCE}" \
        -V "${vcf_file}" \
        --select-type-to-include SNP \
        -O "${SNP_OUTPUT_DIR}/${base_name}_snp.vcf.gz"

    # GATK SelectVariants 对单个染色体 vcf 文件拆分 Indel
    gatk SelectVariants \
        -R "${REFERENCE}" \
        -V "${vcf_file}" \
        --select-type-to-include INDEL \
        -O "${INDEL_OUTPUT_DIR}/${base_name}_indel.vcf.gz"
done
wait


# 4. 成功拆分结果
# chr1.recode_indel.vcf.gz
# chr1.recode_snp.vcf.gz
```




