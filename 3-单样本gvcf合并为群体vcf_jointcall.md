
## 1. gvcf 文件 和 vcf 文件区别
- 单样本 - 仅存在 gvcf 
- 多样本 - 可同时存在 gvcf 和 vcf
- gvcf 文件：含基因组变异和非变异信息，用于多样本联合分析
- vcf 文件：仅含变异位点，适用于结果展示和下游分析
- vcf 文件基于多个 gvcf 文件产生  
```
# 查看命令
bcftools view -h ./jointcall/split/chr2.recode.vcf.gz | head  # 查看头部信息
bcftools view -H ./jointcall/split/chr2.recode.vcf.gz | head  # 查看数据部分，染色体编号未发生变化
```
  
### 1.1 gVCF（genomic VCF）格式
- 单样本分析（多样本联合变异检测_joint genotyping 输入文件
- 格式与 vcf 文件类似，用于同时记录变异位点和非变异位点信息
- 非变异区域会被压缩成一个“区间”表示，称为 reference blocks；
- 每个位点（或区间）都有一个质量分数，说明其为“非变异”的可信度。
- 使用 <NON_REF> 占位符描述参考等位基因以外的所有潜在变异；
- 多样本、大群体变异信息文件也可存在 gvcf 格式

```
#1. 头部文件内容

##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the reference block">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the reference block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM  POS     ID      REF     ALT         QUAL    FILTER  INFO            FORMAT                  Sample1               Sample2               Sample3

# 2. 数据部分内容

chr1    123456  .       A       G,<NON_REF>    50      PASS    .               GT:DP:GQ:PL           0/1:100:50:100,0,200    1/1:90:60:300,60,0     0/0:80:99:0,120,1800
chr1    123457  .       C       <NON_REF>      .       .       END=123788      GT:DP:GQ:MIN_DP:PL    0/0:95:99:85:0,120,1800 0/0:90:99:90:0,120,1800 0/0:80:99:80:0,120,1800
chr1    123789  .       T       C,<NON_REF>    99      PASS    .               GT:DP:GQ:PL           0/0:85:99:0,120,1800    0/1:80:90:100,0,100     0/1:75:85:100,0,200
chr2    234000  .       A       <NON_REF>      .       .       END=234566      GT:DP:GQ:MIN_DP:PL    0/0:70:99:65:0,120,1800 0/0:65:99:60:0,120,1800 0/0:75:99:70:0,120,1800
chr2    234567  .       G       A,<NON_REF>    70      PASS    .               GT:DP:GQ:PL           1/1:90:80:300,80,0      0/1:85:85:100,0,100     0/1:95:87:100,0,100
```
- ALT = G,<NON_REF>：这个位置有替代等位基因 G，并记录 <NON_REF>（为 joint genotyping 准备）；
- <NON_REF> 表示“我们看过这里，没有发现任何变异”；
- END=123788 是这个区间的结束位置；
- GT=0/0 表示基因型是纯合参考
- 0/1:100:50:100,0,200  1/1:90:60:300,60,0  0/0:80:99:0,120,1800 表示三个样本对应的三种基因型**基因型似然得分**，值越小可能性越高
- 0/1:100:50:100,0,200：“100,0,200”从左到右依次为100 0 200，分别代表纯合参考0/0、杂合0/1、纯合突变1/1，中间为0，所以该样本为杂合0/1
- 1/1:90:60:300,60,0：最右侧为0，该样本为纯合突变1/1

### 1.2 vcf（Variant Call Format）文件
**用于存储基因组变异信息的标准文件格式，由头部（Header）和变异数据（Body）两部分组成**    

```
# 1. 头部信息

##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=1,length=248956422>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1 sample2
```
- CHROM：染色体名称或序列名称（例如 1, X, chr2 等）。
- POS：变异的位置，指示变异在染色体上的起始位置。
- ID：变异的标识符，通常是 dbSNP 的 rsID，如 rs12345。如果没有标识符，则用 . 表示。
- REF：参考基因组中该位置的碱基或碱基序列。
- ALT：变异的碱基或碱基序列，可能是单一的替代碱基，也可以是多种碱基的组合。
- QUAL：质量值，表示变异的置信度。较高的值通常表示较高的可信度。
- FILTER：该变异是否通过过滤。常见的值是 PASS，表示该变异通过了质量控制；其他值（如 LowQual）表示未通过过滤。
- INFO：变异的其他信息，通常是以键值对形式表示，例如变异深度（DP）、变异类型（Type）等。
- FORMAT：指定样本基因型信息的格式。通常包含诸如基因型（GT）、深度（DP）、基因型质量（GQ）等字段。
- 样本列（sample1, sample2, ...）：每个样本的基因型信息。基因型通常以 0/0（纯合参考）、0/1（杂合）、1/1（纯合变异）等表示，同时冒号后显示基因频率（1/1:20,80），多个样本之间由制表符分隔。
```
# 2. 数据部分:列出每个变异的具体信息。每一行对应一个变异，列按照头部的列标题顺序排列。

##fileformat=VCFv4.2
##source=MyVariantCaller
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO        FORMAT          Sample1         Sample2         Sample3
chr1    123456  .       A       G       50      PASS    DP=100      GT:AD           0/1:60,40      1/1:20,80      0/0:90,10
chr1    123789  .       T       C       99      PASS    DP=80       GT:AD           0/0:70,10      0/1:40,40      0/1:50,30
chr2    234567  .       G       A       70      PASS    DP=95       GT:AD           1/1:15,80      0/1:50,45      0/1:60,35
```

### 1.3 具体实例
```
# gvcf.gz 单样本
bcftools view -H ./M99.gvcf.gz | head
# CM001582.2	1	.	G	<NON_REF>	.	.	END=4	GT:DP:GQ:MIN_DP:PL	0/0:10:0:10:0,0,0
# CM001582.2	5	.	A	<NON_REF>	.	.	END=5	GT:DP:GQ:MIN_DP:PL	0/0:14:24:14:0,24,360
# CM001582.2	6	.	A	<NON_REF>	.	.	END=6	GT:DP:GQ:MIN_DP:PL	0/0:14:30:14:0,30,450

# gvcf.gz 多样本
bcftools view -H ./qmrb_50sheep.gvcf.gz | head
# CM001582.2	1	.	G	<NON_REF>	.	.	END=3	GT:DP:GQ:MIN_DP:PL	./.:12:0:12:0,0,0	./.:11:0:11:0,0,0	./.:11:0:11:0,0,0 ./.:9:0:9:0,0,0	./.:7:0:7:0,0,0	./.:24:0:24:0,0,0	./.:20:0:20:0,0,0	./.:16:0:16:0,0,0	./.:8:0:8:0,0,0	./.:22:0:22:0,0,0	./.:8:0:8:0,0,0	./.:17:0:17:0,0,0	./.:5:0:5:0,0,0	./.:26:0:26:0,0,0	./.:7:0:7:0,0,0	./.:11:0:11:0,0,0	./.:8:0:8:0,0,0	./.:14:0:14:0,0,0	./.:9:0:9:0,0,0	./.:12:0:12:0,0,0	./.:13:0:13:0,0,0	./.:9:0:9:0,0,0	./.:10:0:10:0,0,0	./.:19:0:19:0,0,0	./.:14:0:14:0,0,0	./.:11:0:11:0,0,0	./.:9:0:9:0,0,0	./.:8:0:8:0,0,0	./.:17:0:17:0,0,0	./.:13:0:13:0,0,0	./.:15:0:15:0,0,0	./.:13:0:13:0,0,0	./.:14:0:13:0,0,0	./.:19:0:19:0,0,0	./.:18:0:18:0,0,0	./.:16:0:16:0,0,0 ./.:12:0:12:0,0,0	./.:8:0:8:0,0,0	./.:6:0:6:0,0,0	./.:9:0:9:0,0,0	./.:12:0:12:0,0,0	./.:14:0:14:0,0,0	./.:14:0:14:0,0,0	./.:10:0:10:0,0,0	./.:17:0:17:0,0,0	./.:15:0:15:0,0,0	./.:8:0:8:0,0,0	./.:5:0:5:0,0,0	./.:9:0:9:0,0,0	./.:10:0:10:0,0,0

# vcf.gz 多样本
bcftools view -H ./qmrb_50sheep.vcf.gz | head
# CM001582.2	84	.	A	AC	10042	.	AC=49;AF=0.49;AN=100;BaseQRankSum=-0.242;DP=2331;ExcessHet=128.468;FS=3.503;MLEAC=49;MLEAF=0.49;MQ=59.77;MQRankSum=0;QD=10.79;ReadPosRankSum=-0.148;SOR=0.484	GT:AD:DP:GQ:PL	0/1:13,10:23:99:353,0,511	0/1:7,5:12:99:165,0,269	0/1:7,10:17:99:341,0,264	0/1:9,8:17:99:277,0,344	0/1:8,4:12:99:126,0,324	0/1:25,9:34:99:266,0,1018	0/1:11,4:15:99:121,0,450	0/1:15,11:26:99:366,0,597	0/1:8,4:12:99:126,0,324	0/1:17,12:29:99:396,0,678	0/1:19,4:23:99:111,0,781	0/1:16,9:25:99:290,0,645	0/1:7,2:9:54:54,0,273	0/1:26,7:33:99:183,0,1071	0/1:7,5:12:99:166,0,274	0/1:16,6:22:99:195,0,644	0/0:15,1:16:3:0,3,630	0/1:17,6:23:99:172,0,676	0/1:6,5:11:99:168,0,237	0/1:13,4:17:99:120,0,529	0/1:10,5:15:99:166,0,405	0/1:13,6:19:99:190,0,528	0/1:16,3:19:64:64,0,658	0/1:15,8:23:99:272,0,606	0/1:10,6:16:99:199,0,402	0/1:13,4:17:99:111,0,529	0/1:4,6:10:99:212,0,150	0/1:5,7:12:99:247,0,189	0/1:8,5:13:99:168,0,321	0/1:15,4:19:99:118,0,603	0/1:9,3:12:85:85,0,369	0/1:7,11:18:99:389,0,261	0/1:16,5:21:99:143,0,657	0/1:20,7:27:99:176,0,819	0/1:26,11:37:99:351,0,1054	0/1:13,6:19:99:184,0,528	0/1:18,6:24:99:170,0,723	0/1:7,8:15:99:283,0,270	0/1:4,6:10:99:217,0,145	0/1:9,3:12:85:85,0,369	0/1:20,6:26:99:173,0,822	0/1:14,6:20:99:182,0,570	0/1:21,11:32:99:357,0,844	0/1:9,6:15:99:206,0,355	0/1:15,6:21:99:188,0,607	0/1:11,8:19:99:266,0,433	0/1:18,6:24:99:164,0,728	0/1:6,7:13:99:244,0,231	0/1:16,5:21:99:148,0,647	0/1:5,5:10:99:181,0,195
```


## 2. joint call
将经过比对、变异检测生成的单样本 gvcf.gz 文件，合并为多样本、大群体类型的 vcf.gz 文件

### 2.1 传统方法
  
  
### 2.2 赛乐服务器方式
```
# 1. 全部个体 gvcf.gz 文件路径
# 保存至输入文件：jointcall.list
find /workspace/public/zongzhan/qmrb_sheep/snp_calling -type f -name "*.gvcf.gz" > jointcall.list

# 样本数量
wc -l test.list 

# jointcall.list 内容如下：
# /workspace/public/zongzhan/qmrb_sheep/snp_calling/G31.gvcf.gz
# /workspace/public/zongzhan/qmrb_sheep/snp_calling/G33.gvcf.gz
# ..................
# /workspace/public/zongzhan/qmrb_sheep/snp_calling/G28.gvcf.gz


# 2. 若通过 Windows 手动输入 jointcall.list 文件，可能需要额外转为 unix 格式
dos2unix jointcall.list


# 3. 合并（joint call）
# -o 输出 gvcf.gz 文件，-g 输出 vcf.gz 文件
# 内存允许尽量同时保存  gvcf.gz 和 vcf.gz 文件，vcf.gz 可直接用于下游分析，如果需增添群体样本，则需要 gvcf.gz 文件
/data/saile/cmd/slmgvcf_gpu  $(cat jointcall.list) -o /workspace/public/zongzhan/qmrb_sheep/jointcall/qmrb_50sheep.gvcf.gz -g /workspace/public/zongzhan/qmrb_sheep/jointcall/qmrb_50sheep.vcf.gz

# 4. 成功获得 qmrb_50sheep.vcf.gz 和 qmrb_50sheep.gvcf.gz 文件
```





















