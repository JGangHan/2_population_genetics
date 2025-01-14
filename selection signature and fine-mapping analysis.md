
--recode 重新输出文件
### 1. 群体文件和质控

### 
population.vcf 需要经过质控
population 划分为 case group（目标性状组）和 control group（对照性状组）
在 population.vcf 文件中提取 case control 对应样本子集数据
# --recode 重新输出文件 --out 输出文件路径和前缀
vcftools --vcf /PATH/TO/population_qc.vcf --keep  /PATH/TO/id_case.txt --recode --out /PATH/TO/trait_case        # 输出 trait_case.recode.vcf 文件
vcftools --vcf /PATH/TO/population_qc.vcf --keep  /PATH/TO/id_control.txt --recode --out /PATH/TO/trait_control  # 输出 trait_control.recode.vcf 文件

计算多态性指数（π）和 π ratio
# --window-pi 窗口大写，--window-pi-step 窗口步长，单位 bp
vcftools --vcf /PATH/TO/trait_case.recode.vcf --window-pi 100000 --window-pi-step 15000 --out /PATH/TO/trait_case        # 输出 trait_case.windowed.pi
vcftools --vcf /PATH/TO/trait_control.recode.vcf --window-pi 100000 --window-pi-step 15000 --out /PATH/TO/trait_control  # 输出 trait_control.windowed.pi
.windowed.pi 文件内容
CHROM	BIN_START	BIN_END	N_VARIANTS	PI
1	1	100000	1264	0.00464854
1	15001	115000	1041	0.00378065
1	30001	130000	848	0.0027974
1	45001	145000	865	0.00282258
1	60001	160000	745	0.0023026
1	75001	175000	659	0.00185558
1	90001	190000	644	0.00178819

# 该步骤为升级代码，需要验证可靠性
# 筛选 case 和 control 重叠区域
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a) { print $1"\t"$2 }' /PATH/TO/trait_case.windowed.pi /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_shared_regions
# 提取 case group 重叠区域
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_case.windowed.pi > /PATH/TO/pi_case_overlap
# 提取 control group 重叠区域
awk -F '\t' 'NR==FNR { a[$1"\t"$2]; next }($1"\t"$2 in a)' /PATH/TO/pi_shared_regions /PATH/TO/trait_control.windowed.pi > /PATH/TO/pi_control_overlap
# 提取 case group 第五列，即 PI value
awk '{ print $5 }' /PATH/TO/pi_case_overlap > /PATH/TO/pi_case_val
# 合并顺序 control group + case pi value，共计6列
paste -d '\t' /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val >  /PATH/TO/pi_merge
# 删除标题行
sed -i '1d' /PATH/TO/pi_merge
# 计算 pi ration (control/case)
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($5/$6)}' /PATH/TO/pi_merge > /PATH/TO/trait_pi_merge_ratio
# 删除过程文件
rm /PATH/TO/pi_shared_regions /PATH/TO/pi_case_overlap /PATH/TO/pi_control_overlap /PATH/TO/pi_case_val

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"($5/$6)}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/a > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case_pi-ratio


```
# 源码是正确的，就是逻辑不太清晰
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








vcftools --vcf /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.recode.vcf --counts --out /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case_counts
vcftools --vcf /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.recode.vcf --counts --out /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control_counts
awk '{print $1"\t"$2"\t"$5"\t"$6}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case_counts.frq.count > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count
awk '{print $1"\t"$2"\t"$5"\t"$6}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control_counts.frq.count  > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count
sed -i 's/A://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count
sed -i 's/T://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count
 sed -i 's/C://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count
sed -i 's/G://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count
sed -i 's/A://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count
sed -i 's/T://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count
sed -i 's/C://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count
sed -i 's/G://g' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count
awk -F '\t' '{if($3>=$4){print $1"\t"$2"\t"$3"\t"$4} else {print $1"\t"$2"\t"$4"\t"$3}}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count.input
awk -F '\t' '{if($3>=$4){print $1"\t"$2"\t"$3"\t"$4} else {print $1"\t"$2"\t"$4"\t"$3}}' /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count > /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count.input
perl /data/liangbm/newsheep/liangbm/keti2/zyhtzhp/calz.pl /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/case.count.input 100000 15000
perl /data/liangbm/newsheep/liangbm/keti2/zyhtzhp/calz.pl /data/liangbm/newsheep/liangbm/ancient/Analysis/390new/344/fst/control.count.input 100000 15000












