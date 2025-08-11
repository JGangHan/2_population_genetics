
- **各 vcf 文件一致性检查**
- **合并质控后各染色体 .vcf.gz 文件**
- **染色体编号重命名**



## 1. 检查 vcf文件相关信息
- **个体名称和数量**
```
# 直接输出样品名称
bcftools query -l qmrb_50sheep_snp.vcf
# 样品数量
bcftools query -l qmrb_50sheep_snp.vcf | wc -l
```

- **头部信息染色体编号**
```
grep "^##contig=" ./qmrb_50sheep_snp.vcf
```

- **数据部分染色体编号**(若处理数个群体 vcf 文件，最好都检查一遍)  
```
# cut -f1  提取第一个字段
cut -f1 /workspace/public/zongzhan/qmrb_sheep/qmrb_50sheep_snp.vcf | grep -v "^#" | sort | uniq
```

## 2. 合并各染色体 vcf.gz 文件（理论上是）
- 经染色体拆分、变异类型拆分、gatk SNP 质控、vcftools SNP 质控后的逐条染色体 vcf.gz 文件
- **前期对染色体逐条处理，理论上 vcf 文件数据部分不存在 contig 信息，而头部信息中存在**
- **理论上每条染色体 vcf 文件头部信息一致，所以仅基于 chr1 染色体的头部信息和数据信息，将其他染色体的数据信息部分逐条添加**
```
#!/bin/bash

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


## 3. 再次检查合并后的 vcf 文件


## 4. vcf 文件染色体重命名
- **将染色体编号修改为标准格式：最好是直接用数字代表染色体（1，2，3...），而不是字符串表示染色体（chr1，chr2，chr3......）** 
- **最好刚获得 vcf 文件就进行处理**  
- **同时修改头部信息和数据部分**


**4.1 染色体编号重命名**  
- **同时修改 vcf 文件头部信息和数据部分染色体格式**
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


**4.2 头部信息中仅保留常染色体**
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
```


**4.3 结果示例** 
```
# 重命名前染色体信息
grep "^##contig=" ./qmrb_50sheep_snp.vcf | head
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

# 重命名后染色体信息
grep "^##contig=" ./qmrb_50sheep.vcf | head
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

# 头部信息末尾多余染色体信息被移除
bcftools view -h ./qmrb_50sheep_chr1to26.vcf.gz | grep "^##contig=" | tail
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
