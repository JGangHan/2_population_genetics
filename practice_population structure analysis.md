从一个大群中提取子群，并进行群体结构分析：PCA analysis、admixture analysis、phylogenetic analysis

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


```





