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

[[ "$db" != "g" && "$db" != "s" ]] && printf "\033[31m Error:数据库参数错误\033[0m\n"
# if [ "$db" = 'g' ];then
# 	echo "该项目使用Greengene数据库注释"
# elif [ "$db" = 's' ]; then
# 	echo "该项目使用Silva数据库注释"
# else
# 	printf "\033[31m Error:数据库参数错误\033[0m\n"
# 	exit 1;
# fi

project=${project%/}   # 去掉斜杠
[ $project = "." ] && project=$(pwd)
# echo "输入的路径为：$project"
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
# ====================================================== #
pre_process(){
	sed -i.bak 's/_[0-9]*//g;s/-/_/g;s/\./_/g' ${project}/01_CleanData/all.fa
}

## ===================== OTU ===========================#
otu(){
	printf "\033[32m OTU聚类\033[0m\n"
	# Usearch 方法
	mkdir -p ${project}/02_OTU
	Usearch11="/home/jbwang/soft/Usearch/Usearch11"
	# sed -i.bak 's/_[0-9]*//g;s/-/_/g;s/\./_/g' ${project}/01_CleanData/all.fa
	$Usearch11 -fastx_uniques ${project}/01_CleanData/all.fa -sizeout -relabel Uniq -fastaout ${project}/01_CleanData/all_norep.fa
	$Usearch11 -cluster_otus ${project}/01_CleanData/all_norep.fa -relabel OTU_ -uparseout ${project}/01_CleanData/uparse.txt -otus ${project}/02_OTU/otus.fa
	# 生成 Otu_table
	$Usearch11 -usearch_global ${project}/01_CleanData/all.fa -db ${project}/02_OTU/otus.fa -otutabout ${project}/02_OTU/otu_table.txt -strand plus -id 0.97
	# biom转换
	biom convert -i ${project}/02_OTU/otu_table.txt -o ${project}/02_OTU/otu_table.biom --table-type="OTU table" --to-json
	# 统计
	biom summarize-table -i ${project}/02_OTU/otu_table.biom -o ${project}/02_OTU/biom_table_summary.txt
	printf "\033[32m OTU聚类完成\033[0m\n"
}

## ===================== 分类 ===========================#
classify(){
	printf "\033[32m 物种注释\033[0m\n"
	silva_fa="/home/jbwang/refedata/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna"
	silva_tax="/home/jbwang/refedata/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt"
	# Greengene base(默认)
	if [ "$db" = "g" ];then
	    assign_taxonomy.py -i ${project}/02_OTU/otus.fa -m rdp -c 0.8 -o ${project}/02_OTU
	else
	    assign_taxonomy.py -i {project}/02_OTU/otus.fa -r $silva_fa -t $silva_tax -m rdp -c 0.8 -o {project}/02_OTU --rdp_max_memory 40000
	fi

	biom add-metadata -i ${project}/02_OTU/otu_table.biom --observation-metadata-fp ${project}/02_OTU/otus_tax_assignments.txt -o ${project}/02_OTU/otu_table_tax.biom --sc-separated taxonomy --observation-header OtuId,taxonomy
	biom convert -i ${project}/02_OTU/otu_table_tax.biom -o ${project}/02_OTU/otu_table_tax.txt --table-type="OTU table" --header-key taxonomy --to-tsv
	printf "\033[32m 物种注释完成\033[0m\n"
}

## ===================== 建树 ===========================#
build_tree(){
	printf "\033[32m 建树\033[0m\n"
	mkdir -p ${project}/02_OTU/backup
	# clustalo软件比对
	clustalo="/home/jbwang/soft/clustalo/clustalo-1.2.4-Ubuntu-x86_64"
	$clustalo -i ${project}/02_OTU/otus.fa -o ${project}/02_OTU/backup/otus.afa --seqtype=DNA --full --force

	# muscle软件比对
	# muscle="/home/jbwang/soft/muscle/muscle3.8.31_i86linux64"
	# $muscle -in ${project}/02_OTU/otus.fa -out ${project}/02_OTU/tmp/otus.afa

	# 筛选结果中保守序列和保守区
	filter_alignment.py -i ${project}/02_OTU/backup/otus.afa -o ${project}/02_OTU/backup/
	# 基于fasttree建树
	make_phylogeny.py -i ${project}/02_OTU/backup/otus_pfiltered.fasta -o ${project}/02_OTU/otus.tre
	printf "\033[32m 建树完成\033[0m\n"
}

## ===================== 多样性 ==========================#
diversity(){
	printf "\033[32m 多样性分析\033[0m\n"
	[ -d ${project}/03_Diversity ] && rm -rf ${project}/03_Diversity
	# 添加#到sample_info.txt中
	sed -i '1s/SampleID/#SampleID/' ${project}/sample_info.txt
	## -e 参数 最小数据量
	alpha_params=/home/jbwang/code/config/alpha_params.txt
	core_diversity_analyses.py -o ${project}/03_Diversity -i ${project}/02_OTU/otu_table_tax.biom -m ${project}/sample_info.txt -t ${project}/02_OTU/otus.tre -e ${min} -p $alpha_params
	# # 用完修改回来
	sed -i '1s/#//' ${project}/sample_info.txt

	# 计算常用的Alpha多样性指数
	min=${min# }   ## 删除空格(莫名出现的空格o(╥﹏╥)o)
	mv ${project}/03_Diversity/arare_max${min} ${project}/03_Diversity/Alpha
	mv ${project}/03_Diversity/bdiv_even${min} ${project}/03_Diversity/Beta
	gzip -d ${project}/03_Diversity/table_even${min}.biom.gz
	mv ${project}/03_Diversity/table_even${min}.biom ${project}/03_Diversity/otu_table_even.biom
	alpha_diversity.py -i ${project}/02_OTU/otu_table_tax.biom -o ${project}/03_Diversity/Alpha/alpha_diversity_index.txt -t ${project}/02_OTU/otus.tre -m shannon,simpson,ace,goods_coverage,chao1,observed_species,PD_whole_tree

	# 计算贝塔多样性距离指数
	mkdir -p ${project}/03_Diversity/Beta/Beta_div
	sed -i '1s/SampleID/#SampleID/' ${project}/sample_info.txt
	weighted_dm=${project}/03_Diversity/Beta/weighted_unifrac_dm.txt
	unweighted_dm=${project}/03_Diversity/Beta/unweighted_unifrac_dm.txt
	make_distance_boxplots.py -d ${weighted_dm} -m ${project}/sample_info.txt -f Group -o ${project}/03_Diversity/Beta/Beta_div/weighted_unifrac --save_raw_data --suppress_all_with --suppress_all_between --suppress_individual_between &
	make_distance_boxplots.py -d ${unweighted_dm} -m ${project}/sample_info.txt -f Group -o ${project}/03_Diversity/Beta/Beta_div/unweighted_unifrac --save_raw_data --suppress_all_with --suppress_all_between --suppress_individual_between &
	wait   ##并行计算

	sed -i '1s/#//' ${project}/sample_info.txt
	printf "\033[32m 多样性分析完成\033[0m\n"
}

## ================= lefse =============================#
lefse(){
	printf "\033[32m LefSe分析\033[0m\n"
	sed -i 's/SampleID/#SampleID/' ${project}/sample_info.txt
	summarize_taxa.py -i ${project}/03_Diversity/otu_table_even.biom -o ${project}/03_Diversity/Beta/lefse/tmp -m ${project}/sample_info.txt --delimiter "|"
	sed -i 's/#//' ${project}/sample_info.txt
	sed -i 's/#*//' ${project}/03_Diversity/Beta/lefse/tmp/*.txt
	python3 /home/jbwang/code/pre_lefse.py ${project}   ## lefse分析前预处理
	for vs in $(ls ${project}/03_Diversity/Beta/lefse |grep vs)
	do
		# vs=${vs##*/}   # 返回斜杠最后的字符
		L6=${project}/03_Diversity/Beta/lefse/${vs}/L6.txt
		out=${project}/03_Diversity/Beta/lefse/${vs}
		awk -F'\t' '{{$2=null;$3=null;$5=null;print $0}}' ${L6} > ${out}/${vs}.txt   ## 删除不必要的列
		sed -i 's/\s\+/\t/g' ${out}/${vs}.txt   # 替换连续空格为tab
		lefse-format_input.py ${out}/${vs}.txt ${out}/${vs}.in -f c -c 2 -s -1 -u 1 -o 1000000
		run_lefse.py ${out}/${vs}.in ${out}/${vs}.res -l 2 > /dev/null
		lefse-plot_res.py ${out}/${vs}.res ${out}/${vs}.png --dpi 600 > /dev/null
		lefse-plot_cladogram.py ${out}/${vs}.res ${out}/${vs}.cla.png --format png --dpi 600 > /dev/null
		lefse-plot_features.py -f diff --archive zip ${out}/${vs}.in ${out}/${vs}.res ${out}/biomarkers.zip > /dev/null
	done
	printf "\033[32m LefSe分析完成\033[0m\n"
}


## ================================================== #
source activate qiime1  ## 激活工作环境

## pre process
if [ ! -e ${project}/pre_process.success ];then
	pre_process
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: pre_process is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/pre_process.success
	fi
fi

## otu
if [ ! -e ${project}/otu.success ];then
	otu
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: otu is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/otu.success
	fi
fi

## check
min=$(grep "Min" ${project}/02_OTU/biom_table_summary.txt|cut -f2 -d ":")
min=${min//,/} ## 替换所有的逗号
min=${min%.000} # 去掉末尾的0
printf "\033[32m 最小测序深度为:$min\033[0m\n"
# minTF=$(echo $min|awk '{if ($0 < 30000) print "True";else print "False"}')
# if [ "$minTF" = "True" ];then
if [ $min -lt 30000 ];then
	printf "\033[33m Warnning:测序最小深度为:$min, 少于三万，请重新测序.\033[0m\n";
	exit 1;
fi

## classify
if [ ! -e ${project}/classify.success ];then
	classify
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: classify is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/classify.success
	fi
fi

## build_tree
if [ ! -e ${project}/build_tree.success ];then
	build_tree
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: build_tree is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/build_tree.success
	fi
fi

## diversity
if [ ! -e ${project}/diversity.success ];then
	diversity
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: diversity is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/diversity.success
	fi
fi

## lefse
if [ ! -e ${project}/lefse.success ];then
	lefse
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: lefse is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/lefse.success
	fi
fi
