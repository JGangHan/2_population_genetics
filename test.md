






```


#################################### 双湖羊 ##################################
# 代码不全，仅挑选关键部分测试运行，具体参考 github
/data/hyr/hjg/raw.SNP.vcf.gz
/data/hyr/hjg/snp448.vcf.gz

###### 按染色拆分snp文件
# 输入路径
INPUT_VCF="/data/hyr/hjg/raw.SNP.vcf.gz"
# 输出路径
OUTPUT_DIR="/data/hyr/hjg/split"
mkdir -p ${OUTPUT_DIR}

### 将vcftools改为bcftools
for i in $(seq 458 483)
do
    OLD_CHR="NC_019${i}.2"
    NEW_NUM=$((i - 457))  # 将染色体编号转换为1到26的范围
    (
        bcftools view -r ${OLD_CHR} ${INPUT_VCF} -O z -o ${OUTPUT_DIR}/chr${NEW_NUM}.recode.vcf.gz &&
        tabix -p vcf ${OUTPUT_DIR}/chr${NEW_NUM}.recode.vcf.gz
    ) & 
done
# 等待所有后台任务完成
wait


####### 使用gatk 或 bcftools 对每条染色体进行过滤，代码不完整
VCF_DIR="/data/hyr/hjg/split"
OUTPUT_DIR="/data/hyr/hjg/split/GATK-filter"
mkdir -p $OUTPUT_DIR

# 1. 使用 GATK 进行过滤
# 只对染色体 10 进行过滤
CHROM=10
(
   gatk --java-options "-Xmx16G" VariantFiltration \
   -V "${VCF_DIR}/chr${CHROM}.recode.vcf.gz" \
   -O "${OUTPUT_DIR}/chr${CHROM}.filtered.vcf.gz" \
   --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" \
   --filter-name "Filter"
)

bcftools view -f PASS "${OUTPUT_DIR}/chr10.filtered.vcf.gz" -Oz -o "${OUTPUT_DIR}/chr10.filtered.pass.vcf.gz"

# 关键代码环节正常，直接跳到合并步骤



##### 合并质控后的每条染色体snp文件
OUT_VCF="/data/hyr/hjg/shuanghu_sheep.vcf"

# 获取 header 并添加（处理 .vcf.gz 文件）
zgrep "^#" /data/hyr/hjg/split/chr1.recode.vcf.gz > "$OUT_VCF"

# 合并各染色体 body
for CHR in {1..26}; do
  FILE="/data/hyr/hjg/split/chr${CHR}.recode.vcf.gz"
  echo "Merging: $FILE"
  # 使用 zgrep 反向检索（去掉 header 部分）并追加到输出文件
  zgrep -v "^#" "$FILE" >> "$OUT_VCF"
done

# 压缩输出文件（确保文件是 .vcf.gz 格式）
bgzip "$OUT_VCF"

# 重新索引（生成 .tbi 索引文件）
tabix -p vcf "${OUT_VCF}.gz"


# 查看样品名称
bcftools query -l /data/hyr/hjg/shuanghu_sheep.vcf.gz
bcftools query -l /data/hyr/hjg/shuanghu_sheep.vcf.gz | wc -l
bcftools query -l /data/hyr/hjg/snp_448.vcf.gz
bcftools query -l /data/hyr/hjg/snp_448.vcf.gz | wc -l








################################# vcf 文件合并前预处理
shuanghu_sheep.vcf
snp448.vcf
# 下方部分步骤跳过，仅测试关键环节代码

# 检查和修改样品名
# 对齐vcf文件染色体名称

# 1. 样本名前添加前缀（可选，避免样本名称重复）
# 两个或多个合并文件都需要检查
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.vcf.gz > samples.txt
# 添加前缀
awk '{print $0 "\tQM" $0}' samples.txt > rename.txt
# 执行样本名替换
bcftools reheader -s rename.txt -o /workspace/qmrb_50sheep_chr1to26.QMprefix.vcf.gz /workspace/mrb_50sheep_chr1to26.vcf.gz
tabix -p vcf /workspace/qmrb_50sheep_chr1to26.QMprefix.vcf.gz
# 重新检查样本名称
bcftools query -l /workspace/qmrb_50sheep_chr1to26.QMprefix.vcf.gz

# 2. 去掉特殊字符-和_(两个连接线)
# 如果存在特殊字符，后续plink软件进行处理和分析会反复复制样本名
grep -E '[-_]' samples.txt
# 若有，去掉
sed -i 's/[-_]//g' samples.txt
# 二次检查
grep -E '[-_]' samples.txt
# 替换 
bcftools reheader -s rename.txt -o /workspace/qmrb_50sheep_chr1to26.QMprefix.vcf.gz /workspace/mrb_50sheep_chr1to26.vcf.gz
# 检查
bcftools query -l /workspace/qmrb_50sheep_chr1to26.QMprefix.vcf.gz


# 3. 提取特定群体（若需要，可选）
# 压缩
bgzip -c /workspace/snp448.rename.vcf > /workspace/snp448.rename.vcf.gz
tabix -p vcf /workspace/snp448.rename.vcf.gz

# 提取前缀为 QL、GG、SG、AW、HB、DM、GB 样本
# 指定样本名
bcftools query -l /workspace/snp448.rename.vcf.gz | grep -E "^(QL|GG|SG|AW|HB|DM|GB|SGTS)" > selected_samples.txt
# 提取
bcftools view -S selected_samples.txt -Oz -o /workspace/snp448_selected.vcf.gz /workspace/snp448.rename.vcf.gz
# 检查
bcftools query -l /workspace/snp448_selected.vcf.gz | wc -l
tabix -p vcf /workspace/snp448_selected.vcf.gz


# 4. 检查样本名是否有重叠
bcftools query -l qmrb_50sheep.sorted.vcf.gz > qm_samples.txt
bcftools query -l snp448.sorted.vcf.gz > chip_samples_renamed.txt
comm -12 <(sort qm_samples.txt) <(sort chip_samples_renamed.txt)  # 若输出为空，则样本无重名


# 5. 对齐头部染色体名称（包含头部和数据部分）
# 格式一：1 2 3 4 5 6 7
# 格式二：chr1 chr2 chr3 chr4 chr5 chr6 chr7
# 格式三（genebank）：CM001600
# 格式四（refseq）：NC_ 或 NW_
bcftools view -h File_A.vcf.gz | grep "^##contig"
bcftools view -h File_B.vcf.gz  | grep "^##contig"

# 例子，将 NC_019${i}.2 格式修改为 12 3 ...26 格式
INPUT_VCF="/data/hyr/hjg/shuanghu_sheep.vcf.gz"
OUTPUT_VCF="/data/hyr/hjg/shuanghu_sheep_modified.vcf.gz"
OUTPUT_DIR="/data/hyr/hjg"
bcftools view -h "$INPUT_VCF" > header.txt
# 使用 sed 修改头部中的染色体名称
sed -i 's/NC_019458.2/1/g' header.txt
sed -i 's/NC_019459.2/2/g' header.txt
sed -i 's/NC_019460.2/3/g' header.txt
sed -i 's/NC_019461.2/4/g' header.txt
sed -i 's/NC_019462.2/5/g' header.txt
sed -i 's/NC_019463.2/6/g' header.txt
sed -i 's/NC_019464.2/7/g' header.txt
sed -i 's/NC_019465.2/8/g' header.txt
sed -i 's/NC_019466.2/9/g' header.txt
sed -i 's/NC_019467.2/10/g' header.txt
sed -i 's/NC_019468.2/11/g' header.txt
sed -i 's/NC_019469.2/12/g' header.txt
sed -i 's/NC_019470.2/13/g' header.txt
sed -i 's/NC_019471.2/14/g' header.txt
sed -i 's/NC_019472.2/15/g' header.txt
sed -i 's/NC_019473.2/16/g' header.txt
sed -i 's/NC_019474.2/17/g' header.txt
sed -i 's/NC_019475.2/18/g' header.txt
sed -i 's/NC_019476.2/19/g' header.txt
sed -i 's/NC_019477.2/20/g' header.txt
sed -i 's/NC_019478.2/21/g' header.txt
sed -i 's/NC_019479.2/22/g' header.txt
sed -i 's/NC_019480.2/23/g' header.txt
sed -i 's/NC_019481.2/24/g' header.txt
sed -i 's/NC_019482.2/25/g' header.txt
sed -i 's/NC_019483.2/26/g' header.txt
# 重新应用修改后的头部
bcftools reheader -h header.txt "$INPUT_VCF" -o "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF"
bcftools view -h "$OUTPUT_VCF" | grep "^##contig" | less 




# 6. 对齐数据部分染色体名称（包含头部和数据部分）
# 因需要合并的不同 vcf 文件染色体明明规则可能不一致，若不一致则肯定会报错，无法合并
# 时间较长
bcftools view -H File_A.vcf.gz | cut -f1 | sort | uniq
bcftools view -H File_B.vcf.gz | cut -f1 | sort | uniq

# 使用 zcat 和 sed 修改数据部分中的染色体名称
zcat shuanghu_sheep_modified.vcf.gz | \
sed 's/NC_019458.2/1/g' | \
sed 's/NC_019459.2/2/g' | \
sed 's/NC_019460.2/3/g' | \
sed 's/NC_019461.2/4/g' | \
sed 's/NC_019462.2/5/g' | \
sed 's/NC_019463.2/6/g' | \
sed 's/NC_019464.2/7/g' | \
sed 's/NC_019465.2/8/g' | \
sed 's/NC_019466.2/9/g' | \
sed 's/NC_019467.2/10/g' | \
sed 's/NC_019468.2/11/g' | \
sed 's/NC_019469.2/12/g' | \
sed 's/NC_019470.2/13/g' | \
sed 's/NC_019471.2/14/g' | \
sed 's/NC_019472.2/15/g' | \
sed 's/NC_019473.2/16/g' | \
sed 's/NC_019474.2/17/g' | \
sed 's/NC_019475.2/18/g' | \
sed 's/NC_019476.2/19/g' | \
sed 's/NC_019477.2/20/g' | \
sed 's/NC_019478.2/21/g' | \
sed 's/NC_019479.2/22/g' | \
sed 's/NC_019480.2/23/g' | \
sed 's/NC_019481.2/24/g' | \
sed 's/NC_019482.2/25/g' | \
sed 's/NC_019483.2/26/g' | \
bgzip >  shuanghu_sheep_modified_body.vcf.gz

# 为输出文件生成索引
tabix -p vcf shuanghu_sheep_modified_body.vcf.gz
bcftools view -H shuanghu_sheep_modified_body.vcf.gz | cut -f1 | sort | uniq



# 7. 提取合并文件共同SNP位点
# 文件A：位点位置
bcftools query -f '%CHROM\t%POS\n' snp448.vcf.gz > snp448_position.txt
# 位置文件可能存在格式错误，必须检查一遍
head snp448_position.txt

# 文件B：提取SNP子集
bcftools query -f '%CHROM\t%POS\n' shuanghu_sheep_modified_body.vcf.gz | wc -l  # 36930565 
gunzip -c shuanghu_sheep_modified_body.vcf.gz > shuanghu_sheep_modified_body.vcf
# 提取子集
vcftools --vcf shuanghu_sheep_modified_body.vcf --positions snp448_position.txt --recode --out shuanghu_common  #475528

# 文件A-B共同位点位置
bcftools query -f '%CHROM\t%POS\n' shuanghu_common.recode.vcf > shuanghu_position.txt
head shuanghu_position.txt
wc -l shuanghu_position.txt

# 文件A：提取子集
gunzip -c snp448.vcf.gz > snp448.vcf
# 提取交集
vcftools --vcf snp448.vcf --positions shuanghu_position.txt --recode --out snp448_common
bcftools query -f '%CHROM\t%POS\n' snp448_common.recode.vcf | wc -l   #475528
# 头部信息
bcftools view -h snp448_common.recode.vcf | less




# 8. SNP 排序
# 默认按照染色体名称（CHROM 字段）和位点位置（POS 字段）进行双重排序
# 如果染色体命名是 "1", "2", ..., "26"，排序将是按照自然数逻辑排序，没有问题。
# 但如果是 "chr1", "chr10", "chr2"，默认是字符串排序，结果可能是：chr1 chr10 chr2
# 所以最好直接用数字形式表示染色体
# 文件A
bcftools sort shuanghu_common.recode.vcf -O v -o shuanghu_common.sorted.vcf
head shuanghu_common.sorted.vcf
# 文件B
bcftools sort snp448_common.recode.vcf -O v -o snp448_common.sorted.vcf
head snp448_common.sorted.vcf


# 9. 转为 ped map 格式，通过 plink 合并
# 格式转换
vcftools --vcf shuanghu_common.sorted.vcf --plink --out test
vcftools --vcf snp448_common.sorted.vcf --plink --out source

# plink 合并
plink --chr-set 26 --file test --merge source.ped source.map --recode --out merged_population
# 输出 merged_population.ped map bed bim fam 等五种文件

# 查看个体和SNP数量
plink --bfile merged_population --freq --out summary






############## 后续 pca 分析和发育树分析代码
# 1. 转换为 VCF 文件格式（供 Beagle 使用）
plink --chr-set 26 --file merged_population --recode vcf --out merged_population

# 2. Beagle 进行缺失基因型填充
# merged_population_bg.vcf 是经过 Beagle imputation/校正后生成的文件，Beagle会重新编码基因型、可能去掉某些注释或合并部分信息，从而减少文件大小。
bgzip merged_population.vcf
# 时间较长
java -Xmx8g -jar /PATH/TO/beagle.01Mar24.d36.jar gt=merged_population.vcf.gz out=merged_population_bg nthreads=8
# 输出 merged_population_bg.vcf.gz

# 3. 检查填充后的 VCF 文件是否存在重复变异位点（可选但推荐）
# plink 可以读入压缩的 vcf.gz 文件，但仅能输出未经压缩 vcf 文件
plink --chr-set 26 --vcf merged_population_bg.vcf.gz --list-duplicate-vars ids-only suppress-first  --out merged_population_remove
head merged_population_remove.dupvar

# 4. 去除重复 SNP 位点，并输出为 vcf 格式
plink --chr-set 26 \
      --vcf merged_population_bg.vcf.gz \
      --exclude merged_population_remove.dupvar \
      --recode vcf \
      --out merged_population_unique
bcftools view -H merged_population_unique.vcf | wc -l
bcftools query -l merged_population_unique.vcf 


# 5. 提取子集（若需）
bcftools query -l merged_population_unique.vcf | grep -E "^(SGTS|QM)" > keep_samples.txt
wc -l keep_samples.txt
bcftools view -S keep_samples.txt merged_population_bg.fixed.vcf -Ov -o subset_SGTS_QM.vcf
bcftools query -l subset_SGTS_QM.vcf 


# 不管是否提取子集，都需要质控
# 6. snp 质控
# 转为 bed bim fam 格式进行过滤
plink --vcf subset_SGTS_QM.vcf --make-bed --out subset_SGTS_QM --chr-set 26

plink --bfile subset_SGTS_QM --chr-set 26 \
  --mind 0.10 \
  --geno 0.10 \
  --make-bed --out subset_qc1_missing

plink --bfile subset_qc1_missing --chr-set 26 \
  --maf 0.05 \
  --make-bed --out subset_qc2_maf

# 检查 snp 过滤情况
wc -l subset_SGTS_QM.bim        # 初始 SNP 数
wc -l subset_qc1_missing.bim                  # 缺失率过滤后 SNP 数
wc -l subset_qc2_maf.bim    


# 7. pca 分析
# 修剪
plink --bfile subset_qc2_maf --double-id --indep-pairwise 50 10 0.2 --out pruned --chr-set 26
plink --bfile subset_qc2_maf \
      --extract pruned.prune.in \
      --make-bed \
      --out subset_qc2_maf_pruned \
      --chr-set 26

# PCA
gcta64 --bfile subset_qc2_maf_pruned --make-grm --autosome-num 26 --out subset_grm_matrix
# 输出文件 grm_matrix.grm.bin grm_matrix.grm.N.bin grm_matrix.grm.id

gcta64 --grm subset_grm_matrix --pca 3 --out subset_pca_matrix
# .eigenval文件，用于计算贡献率，第一个主成分的贡献率=第一个主成分的特征值/总的特征值之和。即文件第一行除以所有行之和，依次类推
# 另一个文件包含样品空间坐标

# 8. 进化树（全部个体或子集）
# 经质控修剪后的文件
plink --bfile merged_qc2_maf_pruned \
      --recode \
      --out merged_qc2_maf_pruned

awk '{print $1, $2}' merged_qc2_maf_pruned.fam | wc -l

# ped2fasta.py 代码见下方
python ped2fasta.py merged_qc2_maf_pruned.ped merged_qc2_maf_pruned.fasta

# 该步骤耗时较长，40W 位点 219 个体 1h
/PATH/TO/rapidNJ-master/bin/rapidnj -i fa mergedt_qc2_maf.fasta -b 1000 > mergedt_qc2_maf.tre
head mergedt_qc2_maf.tre

# 之后用 iTOL 工具可视化


ped2fasta.py

dos2unix ped2fasta.py
chmod +x ped2fasta.py

# 代码内容如下，保存为单独的 ped2fasta.py 文件，最好在工作目录下
#!/usr/bin/env python3

import sys

# 定义将.ped格式转换为.fasta格式的函数
def ped_to_fasta(ped_file, output_file):
    # 打开输入的.ped文件和输出的.fasta文件
    with open(ped_file, 'r') as ped, open(output_file, 'w') as fasta:
        # 逐行读取.ped文件
        for line in ped:
            # 将每一行按空格分割成多个列
            columns = line.strip().split()

            # 第一列是样本ID
            sample_id = columns[0]

            # 从第7列开始是基因型数据
            genotypes = columns[6:]

            # 将基因型数据转换为单一字符串
            # 为简化，我们只取每个SNP的一个等位基因
            # 这将给我们一个单倍体序列
            sequence = ''.join([geno[0] for geno in genotypes])

            # 将转换后的数据写入.fasta文件
            fasta.write(f">{sample_id}\n{sequence}\n")

# 主函数
if __name__ == "__main__":
    # 检查命令行参数的数量
    if len(sys.argv) != 3:
        print("使用方法: python3 脚本名称.py <输入.ped> <输出.fasta>")
        sys.exit(1)

    # 获取输入和输出文件的路径
    ped_file = sys.argv[1]
    output_file = sys.argv[2]

    # 调用函数进行转换
    ped_to_fasta(ped_file, output_file)



#################################### 双湖羊 ##################################
```


















```







# 如何查看 chr1_snp.vcf.gz 文件和 chr1_indel.vcf.gz，这两者有什么区别
# 如何确定质检前SNP数量
# 是否涉及群体相关的标准，是否基于单个染色体进行过滤，如果增减个体是否需要重新过滤
# vcf 文件包含哪些变异类型，是否仅有SNP和Indel两种
# 是否每一次都需要过滤
# gatk和vcftool两种过滤有什么异同
# 为什么对染色体逐条过滤，这是标准分析流程吗
常用质控软件有 gatk、vcftools、plink等 正确吗，过滤过程有什么区别
vcftools 和bcftools 区别




























/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.vcf.gz 中包含 50 个样本，在所有样本名称前添加前缀“QM_”
提取原始样本名列表
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.vcf.gz > samples.txt

生成重命名表（带前缀的样本名）
awk '{print $0 "\tQM_" $0}' samples.txt > rename.txt

执行样本名替换
bcftools reheader -s rename.txt \
  -o /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz \
  /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.vcf.gz
  
tabix -p vcf /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz

# 检查
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz

  


# snp448.vcf 中选择指定样本
cp /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf /workspace/public/zongzhan/qmrb_sheep/snp448.vcf
bgzip -c /workspace/public/zongzhan/qmrb_sheep/snp448.vcf > /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz
tabix -p vcf /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz
# 根据前缀QL、GG、SG、AW、HB、DM、GB, 在 /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf.gz 文件中提取包含上述前缀的样本，并生成新的 .vcf.gz 文件

bcftools query -l /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz > 1.txt

#SG 变为 SG0 
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz \
  | grep -E "^(QL|GG|SG0|AW|HB|DM|GB)" > selected_samples.txt

bcftools view -S selected_samples.txt \
  -Oz -o /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz \
  /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz

wc -l selected_samples.txt  # 148
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | wc -l # 148
tabix -p vcf /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz



###合并前准备工作
# 检查样本名是否有重叠
bcftools query -l qmrb_50sheep.sorted.vcf.gz > qm_samples.txt
bcftools query -l snp448.sorted.vcf.gz > chip_samples_renamed.txt
comm -12 <(sort qm_samples.txt) <(sort chip_samples_renamed.txt)  # 若输出为空，则样本无重名

# 查看信息部分染色体长度是否一致
bcftools view -h /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | grep "^##contig"
bcftools view -h /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz  | grep "^##contig"

# 查看数据部分染色体名称是否一致
bcftools view -H /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | cut -f1 | sort | uniq
bcftools view -H /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz  | cut -f1 | sort | uniq


# 排序
# 按照标准顺序排列 SNP 位点，减少合并过程产生报错的可能

# 默认是按照染色体名称（CHROM 字段）和位点位置（POS 字段）进行双重排序
# 如果染色体命名是 "1", "2", ..., "26"，排序将是按照自然数逻辑排序，没有问题。
# 但如果是 "chr1", "chr10", "chr2"，默认是字符串排序，结果可能是：chr1 chr10 chr2
# 所以最好直接用数字形式表示染色体

# 对芯片数据文件排序
bcftools sort /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz -Oz -o snp448.sorted.vcf.gz
tabix -p vcf snp448.sorted.vcf.gz

# 对重测序文件（你可能已经是排序过的）
bcftools sort /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz -Oz -o qmrb_50sheep.sorted.vcf.gz
tabix -p vcf qmrb_50sheep.sorted.vcf.gz



# 正式合并
bcftools merge -m all -Oz -o merged.vcf.gz snp448.sorted.vcf.gz qmrb_50sheep.sorted.vcf.gz 
tabix -p vcf merged.vcf.gz
报错
# The REF prefixes differ: C vs T (1,1)
# Failed to merge alleles at 1:41555 in snp448.sorted.vcf.gz




######## 另一种合并方式，基于某一个 vcf 文件中的位点位置信息，用 plink 软件提取另一个 vcf 文件中的对应位点相关数据
snp448.sorted.vcf.gz qmrb_50sheep.sorted.vcf.gz 

用vcftools根据list文件提取指定vcf文件，代码
vcftools --gzvcf /data/liyf/00jointcall/data/1052.vcf.gz --keep 30.txt --recode --out 3goat
提取3goatvcf文件数据的snp点，然后根据这个位点信息去提取3个山羊的待鉴定的vcf文件的位点
提取第一列和第二列


# 提取 source 群体位点位置
# 参考群体 snp 数量
bcftools query -f '%CHROM\t%POS\n' snp448.sorted.vcf.gz > snp448_position.txt
wc -l snp448_position.txt # 482709

# 目标群体 snp 数量
bcftools query -f '%CHROM\t%POS\n' qmrb_50sheep.sorted.vcf.gz | wc -l   # 20202525 


# 移除头部的无用信息
# grep -v "^##" 3goat > new_positions_file
# 将两个文件放到一个目录下便于后续操作，后面是目标文件目录，之后解压缩文件
# cp All.raw.snp.vcf.gz /data/wangzq/
# gunzip All.raw.snp.vcf.gz
# 根据这个位置信息提取待鉴定资源的位点


# 先解压 再提取
gunzip -c qmrb_50sheep.sorted.vcf.gz > qmrb_50sheep.sorted.vcf

# 提取 test snp 子集
# 404796
vcftools --vcf qmrb_50sheep.sorted.vcf --positions snp448_position.txt --recode --out qmrb_50sheep_common

# 提取 source snp 子集
bcftools query -f '%CHROM\t%POS\n' qmrb_50sheep_common.recode.vcf > qmrb_position.txt
wc -l qmrb_position.txt
gunzip -c snp448.sorted.vcf.gz > snp448.sorted.vcf
vcftools --vcf snp448.sorted.vcf --positions qmrb_position.txt --recode --out snp448_common

# 404796
bcftools query -f '%CHROM\t%POS\n' snp448_common.recode.vcf | wc -l
bcftools query -l snp448_common.recode.vcf > 1.txt


# 转为 ped map 格式，通过 plink 合并
module load /workspace/public/x86/software/modules/tool/plink-1.90
# 转为 ped map 格式
vcftools --vcf qmrb_50sheep_common.recode.vcf --plink --out test  # 50 个体；404796 snp
vcftools --vcf snp448_common.recode.vcf --plink --out ref   # 148 个体；404796 snp
# 合并
plink --chr-set 26 --file test --merge ref.ped ref.map --recode --out merged_population
# 输出 merged_population.ped map bed bim fam 等五种文件
# 查看个体和SNP数量
plink --bfile merged_population --freq --out summary  # 198 个体，404796 snps


# 
之后将得到的merged_population文件使用begale软件进行填充
java -jar /opt/software/beagle.28Sep18.793.jar gt=all_3goat.vcf out=./allmerge_bg nthreads=8
将重复的位点信息输入到merge-remove.dupvar文件中
plink --chr-set 29 --file all_3goat --list-duplicate-vars ids-only suppress-first -out merge-remove
根据文件去重
plink -chr-set 29 --file all_3goat --exclude merge-remove.dupvar --recode --out all_3goat_unique
将去重后的文件转换成vcf文件
plink --chr-set 29 --file all_3goat_unique --recode-vcf --out all_3goat_unique



# Step 3: 将合并后的PED文件转换成VCF格式（供 Beagle 使用）
plink --chr-set 26 --file merged_population --recode vcf --out merged_population
# Step 4: 使用 Beagle 进行缺失基因型填充
java -Xmx16g -jar /workspace/public/x86/software/tool/beagle-5.4/beagle.01Mar24.d36.jar gt=merged_population.vcf out=merged_population_bg nthreads=16
bcftools query -f '%CHROM\t%POS\n' merged_population_bg.vcf.gz | wc -l

# Step 5: 检查填充后的 VCF 文件是否存在重复变异位点（可选但推荐）
plink --chr-set 26 --vcf merged_population_bg.vcf.gz \
  --list-duplicate-vars ids-only suppress-first \
  --double-id \
  --out merged_population_remove
bcftools query -f '%CHROM\t%POS\n' merged_population_remove.vcf.gz | wc -l


# Step 6: 去除重复 SNP 位点，并输出为 vcf 格式
plink --chr-set 26 \
      --vcf merged_population_bg.vcf.gz \
      --exclude merged_population_remove.dupvar \
      --double-id \
      --recode vcf \
      --out merged_population_unique
bcftools query -f '%CHROM\t%POS\n' merged_population_unique.vcf | wc -l




# 数据修剪 prune
# 转为 bed bim fam 格式，通常要去掉一对性染色体
plink --vcf merged_population_unique.vcf --double-id --make-bed --out merged_population_unique --chr-set 26

# 为什么要 prune
plink --bfile merged_population_unique --double-id --indep-pairwise 50 10 0.2 --out pruned --chr-set 26

--bfile 表示使用 .bed + .bim + .fam 三个文件作为输入；
--file 表示使用 .ped + .map 文件作为输入；
# --indep 50 10 0.2是一个常用的设置，其中：
50：表示窗口大小为50个SNP。
10：表示每次移动窗口时，窗口移动10个SNP。
0.2：表示窗口内的方差膨胀因子（VIF）阈值为0.2。VIF值大于这个阈值的SNP将被去除。
得到两个文件一个out，一个in
?两个文件格式是什么

plink --bfile merged_population_unique \
      --extract pruned.prune.in \
      --make-bed \
      --out merged_population_pruned \
      --chr-set 26

plink --bfile merged_population --freq --out summary





#gcta 和 gatk 傻傻分不清楚
gcta64 --bfile merged_population_pruned --make-grm --autosome-num 26 --out grm_matrix
# 输出文件 grm_matrix.grm.bin grm_matrix.grm.N.bin grm_matrix.grm.id
gcta64 --grm grm_matrix --pca 3 --out pca_matrix
# 输出文件 pca_matrix.eigenval  pca_matrix.eigenvec




计算贡献率
我们将得到的.eigenval文件，用于计算贡献率，第一个主成分的贡献率=第一个主成分的特征值/总的特征值之和。即文件第一行除以所有行之和，依次类推


我们在作图的时候，可能需要对指定地区的品种进行绘图，这时可以提取指定品种作图
awk '$1 == "TGB" ||      $1 == "XiDo" || $1 == "ALS" || $1 == "ABS"      || $1 == "ELS" || $1 == "LiNi" || $1 ==      "TRT" || $1 == "TBG"' 119snp-f-p.ped > Chinese.ped























# 用你原始使用的参考基因组 fasta（例如 Oar_v4.0.fa）
REF="/workspace/public/zongzhan/Reference/Sheep_Oar_4.0/ncbi_dataset/data/GCA_000298735.2/GCA_000298735.2_Oar_v4.0_genomic.fna"
bcftools norm -f $REF -Oz -o snp448.norm.vcf.gz snp448.sorted.vcf.gz
bcftools norm -f $REF -Oz -o qmrb_50sheep.norm.vcf.gz qmrb_50sheep.sorted.vcf.gz

# 索引
tabix -p vcf snp448.norm.vcf.gz
tabix -p vcf qmrb_50sheep.norm.vcf.gz
# 重新合并
bcftools merge -m all -Oz -o merged.vcf.gz qmrb_50sheep.norm.vcf.gz snp448.norm.vcf.gz

























```


