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








# 对 chr9
gatk VariantFiltration \
  -V "${VCF_DIR}/chr9.recode_snp.vcf.gz" \
  -O "${OUTPUT_DIR}/chr9.filtered.vcf.gz" \
  --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "Filter"

# 对 chr15
gatk --java-options "-Xmx4G -Djava.io.tmpdir=./" VariantFiltration \
  -V "${VCF_DIR}/chr15.recode_snp.vcf.gz" \
  -O "${OUTPUT_DIR}/chr15.filtered.vcf.gz" \
  --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "Filter"


# 解压原始 gzip 压缩的 VCF 文件
gunzip -c chr15.filtered.vcf.gz > chr15.vcf











# 加载 vcftools
module load /workspace/public/x86/software/modules/tool/vcftools-0.1.16

# 设置路径
VCF_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter"
OUTPUT_DIR="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter"
mkdir -p $OUTPUT_DIR

# 只处理 chr9 和 chr15
for CHROM in 9 15; do
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






for CHROM in {1..26}; do
  (
    vcftools --gzvcf "${VCF_DIR}/chr${CHROM}.filtered.vcf.gz" \
      --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.9 --min-meanDP 3 \
      --recode --recode-INFO-all --out "${OUTPUT_DIR}/chr${CHROM}.filtered.snp"

    # 删除特定的 contig 信息
    grep -v "##contig=<ID=LWLT*" "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf" > "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf"

    # 修复缺失的基因型格式
    cat "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf" | perl -pe "s/\s\.:/\t.\/.:/g" > "${OUTPUT_DIR}/chr${CHROM}.snp.fix.recode.vcf"

    # 删除中间文件
    rm "${OUTPUT_DIR}/chr${CHROM}.rm.snp.recode.vcf"
    rm "${OUTPUT_DIR}/chr${CHROM}.filtered.snp.recode.vcf" &
  ) &
done
wait









python merge.py list.txt merge.vcf.gz
#去除头信息
bcftools view -H /workspace/public/zongzhan/SLresult/sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/merge.vcf.gz > data_content.txt
#合并新头部和数据部分：
cat /workspace/public/zongzhan/SLresult/sheep/jointcall/head_no_LWLT.txt data_content.txt | bgzip > 98merge.vcf.gz
#为合并后的 VCF 文件生成索引：
bcftools index 98merge.vcf.gz



#!/bin/bash

# 输出文件路径
OUT_VCF="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf"

# 获取 header
grep "^#" /workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr1.snp.fix.recode.vcf > "$OUT_VCF"

# 合并各染色体 body
for CHR in {1..26}; do
  FILE="/workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr${CHR}.snp.fix.recode.vcf"
  echo "Merging: $FILE"
  grep -v "^#" "$FILE" >> "$OUT_VCF"
done
grep "^#CHROM" "$OUT_VCF"
bcftools view -H "$OUT_VCF" | head
bcftools view -h "$OUT_VCF" | head


bcftools view -h ./qmrb_50sheep_snp.vcf | head


bcftools view -h /workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr26.snp.fix.recode.vcf | head

# 样本数量和名称
bcftools query -l /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf
bcftools query -l /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf | wc -l

bcftools query -l /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf
bcftools query -l /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf | wc -l


# 染色体名称

# 统计数据部分染色体名称
cut -f1 /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf | grep -v "^#" | sort | uniq
cut -f1 /workspace/public/zongzhan/qmrb_sheep/jointcall/split/SNP/GATK-filter/vcftools-filter/chr26.snp.fix.recode.vcf | grep -v "^#" | sort | uniq
cut -f1 /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf | grep -v "^#" | sort | uniq

# 统计头部信息中染色体名称
grep "^##contig=" /workspace/public/zongzhan/448Tibetan_sheep/snp448.vcf
grep "^##contig=" ./qmrb_50sheep_snp.vcf




# 染色体重命名
input_file="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf"
output_file="/workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep.vcf"

# CM001582.2    chr1
# CM001607.2    chr26

# 使用sed命令替换染色体名称
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
	
echo "染色体名称替换完成，并保存到 $output_file"



echo "替换前的染色体列表："
cut -f1 "$input_file" | grep -v "^#" | sort | uniq

echo "替换后的染色体列表："
cut -f1 "$output_file" | grep -v "^#" | sort | uniq

# 头部信息
grep "^##contig=" ./qmrb_50sheep.vcf | head



# qmrb_50sheep.vcf.gz

#!/bin/bash
input_vcf="qmrb_50sheep.vcf.gz"
output_vcf="qmrb_50sheep_chr1to26.vcf.gz"

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

tabix -p vcf "$output_vcf"

bcftools view -h ./qmrb_50sheep_chr1to26.vcf.gz | head
zgrep "^##contig=" ./qmrb_50sheep_chr1to26.vcf.gz | head



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

bcftools query -l /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz \
  | grep -E "^(QL|GG|SG|AW|HB|DM|GB)" > selected_samples.txt

bcftools view -S selected_samples.txt \
  -Oz -o /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz \
  /workspace/public/zongzhan/qmrb_sheep/snp448.vcf.gz

bcftools query -l /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | wc -l
tabix -p vcf /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz






合并 /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz 和 /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz 准备步骤
bcftools view -h /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | grep "^##contig"
bcftools view -h /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz  | grep "^##contig"

bcftools view -H /workspace/public/zongzhan/qmrb_sheep/snp448_selected.vcf.gz | cut -f1 | sort | uniq
bcftools view -H /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_chr1to26.QMprefix.vcf.gz  | cut -f1 | sort | uniq



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

# 检查样本名是否有重叠
bcftools query -l qmrb_50sheep.sorted.vcf.gz > qm_samples.txt
bcftools query -l snp448.sorted.vcf.gz > chip_samples_renamed.txt
comm -12 <(sort qm_samples.txt) <(sort chip_samples_renamed.txt)  # 若输出为空，则样本无重名

# 正式合并
bcftools merge -m all -Oz -o merged.vcf.gz snp448.sorted.vcf.gz qmrb_50sheep.sorted.vcf.gz 
tabix -p vcf merged.vcf.gz
报错
# The REF prefixes differ: C vs T (1,1)
# Failed to merge alleles at 1:41555 in snp448.sorted.vcf.gz




snp448.sorted.vcf.gz qmrb_50sheep.sorted.vcf.gz 

用vcftools根据list文件提取指定vcf文件，代码
vcftools --gzvcf /data/liyf/00jointcall/data/1052.vcf.gz --keep 30.txt --recode --out 3goat
提取3goatvcf文件数据的snp点，然后根据这个位点信息去提取3个山羊的待鉴定的vcf文件的位点
提取第一列和第二列


# 提取 source 群体位点位置
bcftools query -f '%CHROM\t%POS\n' snp448.sorted.vcf.gz > snp448_position.txt
# snp 数量
wc -l snp448_position.txt # 482709

bcftools query -f '%CHROM\t%POS\n' qmrb_50sheep.sorted.vcf.gz > 1.txt
wc -l 1.txt # 20202525 

# 移除头部的无用信息
# grep -v "^##" 3goat > new_positions_file
# 将两个文件放到一个目录下便于后续操作，后面是目标文件目录，之后解压缩文件
# cp All.raw.snp.vcf.gz /data/wangzq/
# gunzip All.raw.snp.vcf.gz
# 根据这个位置信息提取待鉴定资源的位点

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


# 转为 ped map 格式，通过 plink 合并
vcftools --vcf qmrb_50sheep_common.recode.vcf --plink --out test
vcftools --vcf snp448_common.recode.vcf --plink --out source


module load /workspace/public/x86/software/modules/tool/plink-1.90 
plink --chr-set 26 --file test --merge source.ped source.map --recode --out merged_population
# 输出 merged_population.ped map bed bim fam 等五种文件
# 查看个体和SNP数量
plink --bfile merged_population --freq --out summary


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
java -Xmx8g -jar /workspace/public/x86/software/tool/beagle-5.4/beagle.01Mar24.d36.jar gt=merged_population.vcf out=merged_population_bg nthreads=8

# Step 5: 检查填充后的 VCF 文件是否存在重复变异位点（可选但推荐）
plink --chr-set 26 --vcf merged_population_bg.vcf.gz \
  --list-duplicate-vars ids-only suppress-first \
  --double-id \
  --out merged_population_remove


# Step 6: 去除重复 SNP 位点，并输出为 vcf 格式
plink --chr-set 26 \
      --vcf merged_population_bg.vcf.gz \
      --exclude merged_population_remove.dupvar \
      --double-id \
      --recode vcf \
      --out merged_population_unique





# 数据修剪 prune
# 转为 bed bim fam 格式，通常要去掉一对性染色体

plink --vcf merged_population_unique.vcf --double-id --make-bed --out merged_population_unique --chr-set 26

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


