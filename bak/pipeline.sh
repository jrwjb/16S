#!/usr/bin/env bash

set -e;  # 出现错误立即退出

help(){
	cat <<- EOF
	Desc: 16S
	Usage: sh pipeline.sh -p 项目路径 -db 注释数据库
	参数选项如下:
	-h --help 打印帮助文档并退出
	-p --project 项目路径
	-db --database (参数为g或s;g为Greengene数据库,s为Silva数据库,若无该参数，默认使用Greengene)
EOF
	exit 0
}
## 不输入参数打印帮助文档
[ $# -lt 1 ] && help

while [ -n "$1" ]; do
	case $1 in
		-p|--project) project=$2
		if [ ! "$project" ];then printf "\033[31m Error:参数错误,未输入路径\033[0m\n" ;exit 1; fi
		shift 2;;  # 将参数后移2个，进入下一个参数的判别
		-db|--database) db=$2
		shift 2;;
		-h|--help) help;; # function help is called
		--) shift;break;; # end of options
		*) printf "\033[31m Error: 没有该参数 $1.\033[0m\n"; exit 1;
	esac
done

# 判断 db 参数
[ ! "$db" ] && db="g"

if [ "$db" = "g" ];then
	echo "该项目使用Greengene数据库注释"
elif [ "$db" = "s" ]; then
	echo "该项目使用Silva数据库注释"
else
	printf "\033[31m Error:数据库参数错误\033[0m\n"
	exit 1;
fi

project=${project%/}   # 去掉斜杠
# ori_path=$(pwd)
[ $project = "." ] && project=$(pwd)
echo "输入的路径为：$project"
# echo $db

# ================================== #
# 程序运行提示
# sp='/-\|'
# sc=0
# spin(){
#     printf "\b${sp:sc++:1}"
#     ((sc==${#sp})) && sc=0
#     sleep 0.1
# }
# ================ 主程序 ================== #
source activate qiime1  ## 激活工作环境

pre_process(){
	sed -i.bak 's/_[0-9]*//g;s/-/_/g;s/\./_/g' ${project}/01_CleanData/all.fa
}

if [ ! -e ${project}/pre_process.success ];then
	pre_process
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: pre_process is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/pre_process.success
	fi
fi

if [ ! -e ${project}/main.success ];then
	bash /home/jbwang/code/main.sh -p $project -db $db &
	bash /home/jbwang/code/picrust/picrust.sh $project &
	wait
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: main programm is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/main.success
	fi
fi

## 画图前备份处理数据
biom convert -i ${project}/03_Diversity/otu_table_even.biom -o ${project}/03_Diversity/otu_table_even.txt --table-type="OTU table" --to-tsv
biom convert -i ${project}/03_Diversity/otu_table_even.biom -o ${project}/03_Diversity/otu_table_even_tax.txt --table-type="OTU table" --header-key taxonomy --to-tsv
sed -i.bak 's/#//' ${project}/02_OTU/otu_table.txt
sed -i.bak '/# Const/d;s/#//' ${project}/02_OTU/otu_table_tax.txt     
sed -i.bak '/# Const/d;s/#//g' ${project}/03_Diversity/otu_table_even.txt
sed -i.bak '/# Const/d;s/#//' ${project}/03_Diversity/otu_table_even_tax.txt

python3 /home/jbwang/code/tax.py $project
python3 /home/jbwang/code/otu.py $project

source deactivate qiime1  ## 关闭工作环境
## ================ 绘图 ===============#
if [ ! -e ${project}/plot.success ];then
	bash /home/jbwang/code/plot.sh $project
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: plot is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/plot.success
	fi
fi

## =============== Result ==============#
bash /home/jbwang/code/result.sh ${project}
