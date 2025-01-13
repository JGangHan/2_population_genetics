从一个大群中提取子群，之后进行 SNP 一般性质控和连锁不平衡质控，最后进行群体结构分析：PCA analysis、admixture analysis、phylogenetic analysis

### 1. plink 软件提取子集
```
# 注：需要重新润色代码，去掉“#”符号和不同行
plink --file population_group       # 输入 population_group.ped 和 population_group.map 文件
      --keep subpopulation1.txt     # 保留 subpopulation1.txt 文件中出现的样品名
      --out subpopulation1          # 输出文件前缀为 subpopulation1
      --recode                      # 输出文件格式为 .ped 和 .map 
      --sheep                       # 指示数据类型为羊（sheep）的基因组数据，可能会启用一些物种特定的处理选项
# 输出文件为 subpopulation1.ped 和 subpopulation1.map 文件
```
subpopulation1.txt 文件中第一列为 Family ID，第二列为 Individual ID
```
# 无表头，两列可以相同也可以不同
4Jm200  4Jm200
4Jm201  4Jm201
4Jm203  4Jm203
4Jm401  4Jm401
4Jm405  4Jm405
4JM500  4JM500
```

### 2. plink 软件对提取子集质控
```
plink --file subpopulation1   # 输入 subpopulation1.ped 和 subpopulation1.map 文件
      --maf 0.05              # 剔除次要等位基因频率（ Minor Allele Frequency）低于 5%（即 0.05）的标记，这些位点变异频率非常低，不具有统计学意义
      --mind 0.1              # 基因型数据缺失比例，0.1 代表每个个体最多允许有 10% 的基因型数据缺失
      --geno 0.1              # 标记缺失程度，0.1 表示如果某个标记的缺失数据超过 10%，那么该标记将会被剔除
      --out subpopulation1_qc # 过滤后输出文件前缀 subpopulation1_qc
      --recode                # 输出文件格式为 .ped 和 .map
      --sheep                 # 物种
# 输出文件为 subpopulation1_qc.ped 和 subpopulation1_qc.map 文件
```


### 3. plink 软件基于连锁不平衡进一步质控
计算 LD
```
plink --file subpopulation1_qc       # 输入 subpopulation1_qc.ped 和 subpopulation1_qc.map 文件
      --indep 50 5 2                 # 删除在群体中存在高LD的变异标记，减少计算复杂度。
                                     # 50：在每个 LD 窗口内考察 50 个变异标记；
                                     # 5：该窗口内的变异标记之间最大距离为 5 Kb；
                                     # 2：当两个标记之间的相关性大于 0.2（通常是 r²）时，删除其中一个标记
      --out subpopulation1_qc_ld     # 过滤后输出文件前缀 subpopulation1_qc_ld
      --recode                      
      --noweb                        # 不进行 Web 更新检查
# 输出文件为 subpopulation1_qc_ld.ped 和 subpopulation1_qc_ld.map 文件；subpopulation1_qc_ld.prune.in 和 subpopulation1_qc_ld.prune.out 文件
```
提取目标 SNP 位点
```
plink --file subpopulation1_qc                 # 输入文件前缀
      --extract subpopulation1_qc_ld.prune.in  # 需要保留的 SNP ID
      --out subpopulation1_qc_prune            # subpopulation1_qc_prune.ped 和 subpopulation1_qc_prune.map 文件
      --recode


# subpopulation1_qc_ld.prune.in 文件内容如下
1:94
1:95
...
26:44046956
26:44046973
26:44047006
26:44047044
```
**上方两行代码感觉非常相似，不确定是不是可以删除某一步骤**



### 4. PCA 分析
```
plink --file subpopulation1_qc_prune     # 输入文件 subpopulation1_qc_prune.ped 和 subpopulation1_qc_prune.map
      --pca 3                            # 计算前3个主成分
      --out subpopulation1_qc_prune      # 输出文件 subpopulation1_qc_prune.eigenval 和 subpopulation1_qc_prune.eigenvec

# subpopulation1_qc_prune.eigenval 中包含每个主成分的特征值（方差贡献），对应3个主成分
11.0884
7.71976
6.40445

# subpopulation1_qc_prune.eigenvec 中包含每个个体在每个主成分上的得分（相当于坐标）
# 从左到右依次为 Family ID，Individual ID，PC1，PC2，PC3
4Jm200 4Jm200 0.0576596 0.201919 -0.0276461
4Jm201 4Jm201 0.0553672 0.211469 -0.028295
4Jm203 4Jm203 0.0498517 0.168678 -0.0202663
4Jm401 4Jm401 0.056855 0.237961 -0.03153
4Jm405 4Jm405 0.0532857 0.178357 -0.0223027
4JM500 4JM500 0.0570752 0.224826 -0.0298816
4JM503 4JM503 0.0544851 0.201213 -0.0279637
```

### 5. MEGA 系统发生树
plink 计算遗传距离
```
plink1 --file subpopulation1_qc_prune   # 输入文件 subpopulation1_qc_prune.ped 和 subpopulation1_qc_prune.map
       --cluster                        # 进行基于遗传距离的样本聚类
       --distance-matrix                # 计算每对样本之间的遗传距离，并输出距离矩阵
       --out subpopulation1             # 输出文件前缀
       --sheep  --noweb
# 输出 subpopulation1.cluster 和 subpopulation1.mdist 文件

# subpopulation1.cluster：包含群体结构聚类的结果
# subpopulation1.mdist：包含样本两两之间的遗传距离矩阵。
```
perl 脚本将遗传距离矩阵数据文件转换为 .mega 格式
```
perl plink.distance.matrix.to.mega.pl            # 
     subpopulation1_prunename.txt                # 所包含的样本名
     subpopulation1.mdist                        # 矩阵文件 
     341                                         # 样本数量
     subpopulation1                              # 输出文件前缀
# 输出文件 subpopulation1.meg 文件



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
4JM1708	4JM1708
```





