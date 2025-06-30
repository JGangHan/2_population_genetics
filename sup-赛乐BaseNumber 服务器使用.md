该服务器内部集成一套自己的生信分析流程，涵盖：测序数参考基因组比对，SNP calling，gvcf merging
### 1. 原始测序数据文件 list
- sample_path.sh，bash 脚本文件，根据测序文件后缀名称输出全部测序文件的绝对路径
```
#!/bin/bash

# 设置目录路径
directory="/workspace/public/zongzhan/qmrb_sheep/JZP202504GM006-01-F001-0/JZP202504GM006-01-F001/00_Clean"

# 声明一个关联数组来存储文件对
declare -A file_pairs

# 根据测序数据后缀命名特征修改下方代码，常见的有
# 1. .clean.R1.fq.gz  .clean.R2.fq.gz
# 2. _1.fq.gz          _2.fq.gz
# 3. .R1.fq.gz         .R2.fq.gz 
# 在对应位置修改代码

# 遍历目录中所有以.clean.R1.fq.gz和.clean.R2.fq.gz结尾的文件
for file in "$directory"/*.clean.R1.fq.gz; do
    # 提取文件的基础编号（去掉.clean.R1.fq.gz后面的部分）
    base_name=$(basename "$file")
    id=${base_name%%.clean.R1.fq.gz}

    # 将文件按编号存入关联数组
    if [[ -n "${file_pairs[$id]}" ]]; then
        # 如果数组中已存在该编号的文件，输出文件对
        echo "${file_pairs[$id]} $file"
    else
        # 存储第一个文件的路径
        file_pairs[$id]=$file
    fi
done

# 遍历目录中所有以.clean.R2.fq.gz结尾的文件
for file in "$directory"/*.clean.R2.fq.gz; do
    # 提取文件的基础编号（去掉.clean.R2.fq.gz后面的部分）
    base_name=$(basename "$file")
    id=${base_name%%.clean.R2.fq.gz}

    # 将文件按编号存入关联数组
    if [[ -n "${file_pairs[$id]}" ]]; then
        # 如果数组中已存在该编号的文件，输出文件对
        echo "${file_pairs[$id]} $file"
    fi
done
```
- 修改文件权限为可执行 bash 脚本
- 如果在Windows中编辑脚本文件，**Unix中可能无法识别“回车符”**，最好通过 dos2unix 命令转换一下
```
chmod +x sample_path.sh  # 修改权限
dos2unix sample_path.sh  # 转为 linux 格式
./sample_path.sh  # 运行
```
- sample_path.txt：输出为目标测序文件的绝对路径，将其保存至sample_path.txt
![8f50f68f8956b3a50584180647201e32_image](https://github.com/user-attachments/assets/5aedb488-be85-4dcd-94eb-5210794eb0cb)


### 2. 建立参考基因组索引
- 新参考基因组必须提前建立索引，时间较长
```
/data/saile/cmd/slaidx /workspace/public/goat_ref/T2T_ref/CAU_T2T/GCA_040806595.1_T2T-goat1.0_genomic.fna
```


### 3. 参考基因组比对并生成 gvcf 文件
- gvcf.sh
```
#bam=""
curDir=$(pwd)  # 保存当前目录的路径到变量curDir

# 切换到脚本所在的目录 $0：这是一个特殊的变量，它在Bash脚本中表示当前执行的脚本的名称,dirname 是一个常用的Unix命令，它用于从完整的文件路径中提取出目录部分。
# 例如，如果你给 dirname 一个完整的文件路径 /path/to/my/script.sh，它将返回 /path/to/my，即文件所在的目录。

cd $(dirname $0) 

# 以下是一些配置变量
path=/workspace/public/zongzhan/qmrb_sheep/snp_calling  # 指定输出目录
file=/workspace/public/zongzhan/qmrb_sheep/sample_path.txt  # 指定输入文件的名称
fasta=/data1/project/goat/reference/Oar_4.0/GCA_000298735.2_Oar_v4.0_genomic.fna   # 指定fasta文件的路径

bamSplit=/data1/project/goat/bam_split  # 定义用于分割BAM文件的路径
exe="/data/saile/cmd"  # 定义执行程序的路径

if [[ ! -d $path ]];then  # 如果输出目录不存在
	mkdir -p $path  # 则创建该目录
fi

cmdLog=/workspace/public/zongzhan/qmrb_sheep/snp_calling.log  # 定义一个临时文件来存储命令日志
echo "------------------" >> $cmdLog  # 在日志文件中添加分隔线

cnt=3  # 设置一个计数变量
i=0  # 初始化一个迭代器
cat $file | while read line  # 逐行读取文件file中的内容
do
	a=($line)  # 将读取的行分割成数组
	r1=${a[0]}  # 第一个元素，通常是读取的第一个fastq文件
	r2=${a[1]}  # 第二个元素，通常是读取的第二个fastq文件
	bed=${a[2]}  # 第三个元素，通常是bed文件的路径

	fastqName=`basename $r1`  # 从第一个fastq文件的路径中提取文件名
	prefix=${fastqName%%.*}  # 从文件名中提取前缀

	if [[ ! -z $bed ]];then  # 如果bed文件路径不为空
		metrics_bed="--metrics-bait-bed-file $bed"  # 设置bed文件的参数
		addbed=" -B $bed"  # 添加-bed参数
		gdepth_bed="--gdepth-bed-file $bed"  # 设置gdepth-bed文件的参数
	fi

	# 接下来是一系列变量和参数的定义，用于后续的命令行操作
	low=""
	slaArgus=""
	illumina_adapter_r1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	illumina_adapter_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	bgi_adapter_r1="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
	bgi_adapter_r2="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"
	fastp_r1=${path}/${prefix}_fastp_r1.fq.gz
	fastp_r2=${path}/${prefix}_fastp_r2.fq.gz
	fastp_output_argus=""
	no_alt=""
	metrics=${path}/${prefix}_sla.metrics

	# 构建命令行命令
	cmd="${exe}/slat -R \"@RG\tID:TEST\tSM:${prefix}\tPL:ILLUMINA\" $no_alt -B $bamSplit $slaArgus $fasta $r1 $r2 $low $fastp_output_argus -i $metrics"
	echo $cmd >> $cmdLog
	bts=`date +%s`
	eval $cmd
	if [[ $? -ne 0 ]];then
		echo "sla exec failed"
		rm -rf $bamSplit
		continue
	fi
	ets=`date +%s`
	slats=`expr $ets - $bts`

	# 接下来的部分构建了更多的命令行命令，用于处理生物信息学数据
	bamName=${path}/${prefix}.slc.bam
	gvcf=${path}/${prefix}${no_alt}.vcf.gz
	vcf_metrics=${path}/${prefix}_slc.metrics
	gdepth_args=""
	#cmd="${exe}/slc --no-mark-supplementary --keep-split$slbcArgus -P $bamSplit -b $bamName -R $fasta -o $gvcf $addbed $low --with-metrics --metrics-file $vcf_metrics --vcf-index $metrics_bed $static_quantized_quals_args $gdepth_args"
	echo $cmd >> $cmdLog
	bts=`date +%s`
	#eval $cmd
	ets=`date +%s`
	slcts=`expr $ets - $bts`

	# 生成GVCF文件的命令
	gvcf=${path}/${prefix}${no_alt}.gvcf.gz
	cmd="${exe}/slc -e GVCF --keep-split --no-mark-supplementary$slbcArgus -P $bamSplit -R $fasta -o $gvcf $addbed --with-metrics --metrics-file $vcf_metrics $metrics_bed --vcf-index $static_quantized_quals_args $gdepth_args"
	echo $cmd >> $cmdLog
	bts=`date +%s`
	eval $cmd

	ets=`date +%s`
	slcgvcfts=`expr $ets - $bts`
	echo "germline,$prefix,$slats,$slcts,$slcgvcfts" >> ret.csv

done

cd $curDir  # 返回最初的目录

```









