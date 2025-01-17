






### 测序深度/覆盖度
```
file="
B22704_DSW63001-V
B24517_DSW63007-V
B24560_DSW63009-V
"
for d in $file;
do
#BWA mem generate sam file
genomeCoverageBed \
        -ibam /home/jianglin/ljiang/xiaohong/10xMongolianSheep/''$d'.sorted.uniqe.rg.dedup.realn.bam' \
        >''$d'.bam.hist.txt'
done
```
.bam.hist.txt 输出结果如下：列1：染色体；列2：覆盖度（从低到高）；列3：计数具有该覆盖度的基因组位置数量；列4：染色体长度；列5：具有该覆盖度的基因组长度占比（0.0458617 表示 chr1 覆盖度为 0 的基因组区域占该染色体的4.6%）
```
chr1	0	12630621	275406953	0.0458617
chr1	1	2184731	275406953	0.00793274
chr1	2	2785809	275406953	0.0101152
              .....
chrX	1836	3	135185801	2.21917e-08
chrX	1837	4	135185801	2.95889e-08
chrX	1838	3	135185801	2.21917e-08
```


### 构建索引
#### 1. tabix 对 .gz 压缩文件构建索引
```
# tabix 要求文件必须是通过 bgzip 压缩的，并对其进行索引。
tabix -p vcf sample.vcf.gz
```
#### 2. samtools index 创建 .bam 二进制文件索引
```
samtools index sample.sorted.rmdup.bam
```

















































