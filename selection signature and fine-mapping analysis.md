

**选择信号分析本质上是利用不同的统计方法评估不同群体之间的分化程度**

### 待解决疑问
- 最多可以指定多少群体进行分析，即同时分析的群体的最大数量


## 1. 群体 vcf 文件和质控？？
准备经质控后的 population_qc.vcf 文件。**具体哪些步骤？？质控到哪些程度？？**


## 2. 选择信号分析
### 2.1 population 划分为 case group（目标性状组）和 control group（对照性状组）
提取 case control control 对应样本子集数据
```
# --recode 重新输出文件 --out 输出文件路径和前缀
vcftools --vcf /PATH/TO/population_qc.vcf --keep  /PATH/TO/id_case.txt --recode --out /PATH/TO/trait_case        # 输出 trait_case.recode.vcf 文件
vcftools --vcf /PATH/TO/population_qc.vcf --keep  /PATH/TO/id_control.txt --recode --out /PATH/TO/trait_control  # 输出 trait_control.recode.vcf 文件
```
id_case.txt 文件内容如下：
```
4Jm200	4Jm200
4Jm201	4Jm201
4Jm203	4Jm203
4Jm401	4Jm401
4Jm405	4Jm405
4JM500	4JM500
4JM503	4JM503
4Jm801	4Jm801
4JM1701	4JM1701
4JM1703	4JM1703
```


### 2.2 群体分化指数 Fst




### 2.3 多态性指数（π）和 π ratio
对case control group 分别计算 pi value
```
# --window-pi 窗口大小，--window-pi-step 窗口步长，单位 bp
vcftools --vcf /PATH/TO/trait_case.recode.vcf --window-pi 100000 --window-pi-step 15000 --out /PATH/TO/trait_case        # 输出 trait_case.windowed.pi
vcftools --vcf /PATH/TO/trait_control.recode.vcf --window-pi 100000 --window-pi-step 15000 --out /PATH/TO/trait_control  # 输出 trait_control.windowed.pi
```
.windowed.pi 文件内容如下
```
CHROM	BIN_START	BIN_END	N_VARIANTS	PI
1	1	100000	1264	0.00464854
1	15001	115000	1041	0.00378065
1	30001	130000	848	0.0027974
1	45001	145000	865	0.00282258
1	60001	160000	745	0.0023026
1	75001	175000	659	0.00185558
1	90001	190000	644	0.00178819
```
计算 π ratio（需明确 case control group 信号一致区域，之后在两个群体中分别提取，最后计算 pi ration）
```
# 该步骤为升级代码，需要验证可靠性

# 筛选 case 和 control group 信号重叠区域
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a) { print $1"\t"$2 }' /PATH/TO/trait_case.windowed.pi /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_shared_regions
# 提取 case group 信号子集
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_case.windowed.pi > /PATH/TO/pi_case_overlap
# 提取 control group 信号子集
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_control_overlap
# 提取 case group 第五列，即 PI value
awk '{ print $5 }' /PATH/TO/pi_case_overlap > /PATH/TO/pi_case_val
# 合并两个群体 pi value 文件：顺序 control group + case pi value，共计6列
paste -d '\t' /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val >  /PATH/TO/pi_merge
# 删除标题行
sed -i '1d' /PATH/TO/pi_merge
# 计算 pi ration (control/case)
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($5/$6)}' /PATH/TO/pi_merge > /PATH/TO/trait_pi_merge_ratio

# 删除过程文件
rm /PATH/TO/pi_shared_regions /PATH/TO/pi_case_overlap /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val
```
上方代码修改自源码，源码是正确的，就是逻辑不太清晰，如下
```
# control
awk -F '\t' 'NR==FNR{a[$1"\t"$2]}NR>FNR{if($1"\t"$2 in a)print}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.windowed.pi /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.windowed.pi > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/1
# case
awk -F '\t' 'NR==FNR{a[$1"\t"$2]}NR>FNR{if($1"\t"$2 in a)print}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/1 /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.windowed.pi > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/2
# 提取 case 第五列
awk '{print $5}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/2 > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/3
# 合并 control+case
paste -d "\t" /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/1 /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/3 > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/a
# 删除第一行
sed -i '1d' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/a
# control/case # 值越大表明在 case 中越受选择
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($5/$6)}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/a > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case_pi-ratio
```
### 2.4 计算 Zhp
对 case group 和 control group 的 VCF 文件进行等位基因计数
```
vcftools --vcf /PATH/TO/trait_case.recode.vcf --counts --out /PATH/TO/trait_count_case
vcftools --vcf /PATH/TO/trait_control.recode.vcf --counts --out /PATH/TO/trait_count_control
```
输出 .frq.count 文件内容如下
```
CHROM	POS	N_ALLELES	N_CHR	{ALLELE:COUNT}
1	87	2	366	A:193	C:173
1	88	2	366	A:191	G:175
1	94	2	366	T:190	C:176
1	95	2	366	G:324	A:42
1	137	2	366	T:219	C:147
1	138	2	366	G:303	A:63
1	191	2	366	C:298	T:68
```
**case group 数据预处理，包括子集提取、字符修改、排序等**
```
# 提取等位基因计数文件 .frq.count 的第1，2，5，6列
awk '{print $1"\t"$2"\t"$5"\t"$6}' /PATH/TO/trait_count_case.frq.count > /PATH/TO/trait_count_case.count

# 删除 ATCG 字符
sed -i 's/A://g' /PATH/TO/trait_count_case.count
sed -i 's/T://g' /PATH/TO/trait_count_case.count
sed -i 's/C://g' /PATH/TO/trait_count_case.count
sed -i 's/G://g' /PATH/TO/trait_count_case.count

# 针对第3，4列数值大小逐行进行排序，但是否会丢失等位基因信息？？？
awk -F '\t' '{if($3>=$4){print $1"\t"$2"\t"$3"\t"$4} else {print $1"\t"$2"\t"$4"\t"$3}}'/PATH/TO/trait_count_case.count > /PATH/TO/trait_count_case.count.input
```
trait_count_case.frq.count（6列）转为 trait_count_case.count（4列），同时删除 ATCG 等字符串，并针对第3、4列逐行进行排序
```
CHROM	POS	{ALLELE:COUNT}	
1	87	193	173
1	88	191	175
1	94	190	176
1	95	324	42
1	137	219	147
1	138	303	63
1	191	298	68
1	198	291	75
1	238	322	44
```
**计算 Zhp**
```
perl calz.pl /PATH/TO/trait_count_case.count.input 100000 15000
# 输出文件？？
```
**control group 数据预处理，与上方类似**
```
awk '{print $1"\t"$2"\t"$5"\t"$6}' /PATH/TO/control_counts.frq.count > /PATH/TO/trait_count_control.count
sed -i 's/A://g' /PATH/TO/trait_count_control.count
sed -i 's/T://g' /PATH/TO/trait_count_control.count
sed -i 's/C://g' /PATH/TO/trait_count_control.count
sed -i 's/G://g' /PATH/TO/trait_count_control.count
awk -F '\t' '{if($3>=$4){print $1"\t"$2"\t"$3"\t"$4} else {print $1"\t"$2"\t"$4"\t"$3}}' /PATH/TO/trait_count_control.count > /PATH/TO/trait_count_control.count.input
perl calz.pl /PATH/TO/trait_count_control.count.input 100000 15000
```




## 3. fine mapping
**对某一目标基因上下游区域发生的选择信号进行精细定位**
以 PDGFD 基因为例，指定PDGFD基因上下游目标区域
```
vcftools --vcf /PATH/TO/population_qc.vcf --chr 15 --from-bp 3785000 --to-bp 3985000 --recode --out /PATH/TO/trait_pdgfd_200kb     
# 输出 trait_pdgfd_200kb.recode.vcf 文件
```
### 3.1 群体分化指数 Fst
```
# 不指定窗口和步长，对每一碱基逐帧计算
vcftools --vcf /PATH/TO/trait_pdgfd_200kb.recode.vcf --weir-fst-pop id_case.txt --weir-fst-pop id_control.txt --out /PATH/TO/trait_pdgfd_fst  
```
输出 trait_pdgfd_fst.weir.fst 文件内容如下
```
CHROM	POS	WEIR_AND_COCKERHAM_FST
15	3785025	0.180823
15	3785026	0.0163474
15	3785166	0.198976
15	3785218	0.282564
15	3785252	0.27378
15	3785373	0.0718726
15	3785436	0.0819801
15	3785460	0.0886902
15	3785488	0.0914016
```
### 3.2 多态性指数（π）和 π ratio
先提取子集，之后分别计算 π value，合并文件后计算π ratio
```
# 提取子集 case control group
vcftools --vcf /PATH/TO/trait_pdgfd_200kb.recode.vcf --keep id_case.txt --recode --out /PATH/TO/trait_pdgfd_200kb_case
vcftools --vcf /PATH/TO/trait_pdgfd_200kb.recode.vcf --keep id_control.txt --recode --out /PATH/TO/trait_control_200kb_control

# 对 case control group 分别计算 π value
vcftools --vcf /PATH/TO/trait_pdgfd_200kb_case.recode.vcf --window-pi 20000 --window-pi-step 20000 --out  /PATH/TO/trait_pdgfd_200kb_case            # 输出 trait_pdgfd_200kb_case.windowed.pi
vcftools --vcf  /PATH/TO/trait_pdgfd_200kb_control.recode.vcf --window-pi 20000 --window-pi-step 20000 --out  /PATH/TO/trait_pdgfd_200kb_control     # 输出 trait_pdgfd_200kb_case.windowed.pi

# 计算 π ratio（需明确 case control group 信号一致区域，之后在两个群体中分别提取，最后计算 pi ratio

# 筛选 case 和 control group 信号重叠区域
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a) { print $1"\t"$2 }' /PATH/TO/trait_case.windowed.pi /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_shared_regions
# 提取 case group 信号子集
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_case.windowed.pi > /PATH/TO/pi_case_overlap
# 提取 control group 信号子集
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_control_overlap
# 提取 case group 第五列，即 PI value
awk '{ print $5 }' /PATH/TO/pi_case_overlap > /PATH/TO/pi_case_val
# 合并两个群体 pi value 文件：顺序 control group + case pi value，共计6列
paste -d '\t' /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val >  /PATH/TO/pi_merge
# 删除标题行
sed -i '1d' /PATH/TO/pi_merge
# 计算 pi ration (control/case)
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($5/$6)}' /PATH/TO/pi_merge > /PATH/TO/trait_pdgfd_pi_merge_ratio

# 删除过程文件
rm /PATH/TO/pi_shared_regions /PATH/TO/pi_case_overlap /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val
```
### 3.3 TajimaD 
对提取的 case control group 先分别计算 TajimaD value，预处理后计算两个群体之间的差值 control group - case group
```
# 20k 窗口，若不指定步长，则窗口与步长同等大小
vcftools --vcf /PATH/TO/trait_pdgfd_200kb_case.recode.vcf --TajimaD 20000 --out /PATH/TO/trait_pdgfd_200kb_case           # 输出 trait_pdgfd_200kb_case.Tajima.D 文件
vcftools --vcf /PATH/TO/trait_pdgfd_200kb_control.recode.vcf --TajimaD 20000 --out /PATH/TO/trait_pdgfd_200kb_control     # 输出 trait_pdgfd_200kb_control.Tajima.D 文件
```
.Tajima.D 文件内容如下
```
# case group
CHROM	BIN_START	N_SNPS	TajimaD
15	3780000	188	-0.583854
15	3800000	205	-0.55682
15	3820000	137	-1.01208
15	3840000	136	-0.729032
15	3860000	254	-0.235477
15	3880000	202	-0.955357
15	3900000	235	-0.0111044
15	3920000	209	1.52266
15	3940000	274	3.4678
15	3960000	221	2.52908
15	3980000	33	2.27151

# control group
CHROM	BIN_START	N_SNPS	TajimaD
15	3780000	188	2.74521
15	3800000	205	4.17985
15	3820000	137	2.3483
15	3840000	136	3.09756
15	3860000	255	4.53153
15	3880000	203	3.08861
15	3900000	235	3.68404
15	3920000	209	3.52589
15	3940000	274	4.12145
15	3960000	221	3.74258
15	3980000	33	3.8598
```
提取两个群体信号重叠区域，预处理并计算差值
```
# control overlap
awk -F '\t' 'NR==FNR{a[$1"\t"$2]}NR>FNR{if($1"\t"$2 in a)print}' /PATH/TO/trait_pdgfd_200kb_case.Tajima.D /PATH/TO/trait_pdgfd_200kb_control.Tajima.D > PATH/TO/TD_control
# case overlap
awk -F '\t' 'NR==FNR{a[$1"\t"$2]}NR>FNR{if($1"\t"$2 in a)print}' PATH/TO/TD_control PATH/TO/trait_pdgfd_200kb_case.Tajima.D > PATH/TO/TD_case
# case TajimaD value
awk '{print $4}' PATH/TO/TD_case > PATH/TO/TD_case_value
# 合并 control + case TD value
paste -d "\t"  PATH/TO/TD_control PATH/TO/TD_case_value > PATH/TO/TD_merge

# 删除首行
sed -i '1d' PATH/TO/TD_merge
# 计算两组差值：control - case
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($4-$5)}' PATH/TO/TD_merge > PATH/TO/trait_pdgfd_200kb_TD.txt

# 删除过程文件
rm /PATH/TO/trait_pdgfd_200kb_case.Tajima.D /PATH/TO/trait_pdgfd_200kb_control.Tajima.D /PATH/TO/TD_control /PATH/TO/TD_case PATH/TO/TD_merge
```
trait_pdgfd_200kb_TD.txt 文件内容如下
```
15	3780000	188	2.74521	-0.583854	3.32906
15	3800000	205	4.17985	-0.55682	4.73667
15	3820000	137	2.3483	-1.01208	3.36038
15	3840000	136	3.09756	-0.729032	3.82659
15	3860000	255	4.53153	-0.235477	4.76701
15	3880000	203	3.08861	-0.955357	4.04397
15	3900000	235	3.68404	-0.0111044	3.69514
15	3920000	209	3.52589	1.52266	2.00323
15	3940000	274	4.12145	3.4678	0.65365
15	3960000	221	3.74258	2.52908	1.2135
15	3980000	33	3.8598	2.27151	1.58829
```










