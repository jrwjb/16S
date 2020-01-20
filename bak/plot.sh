#!/usr/bin/env bash

set -e;  # 出现错误立即退出

path=$1
path=${path%/}     #  '%' 从后向前删除, '#' 从前向后删除

heatmap(){
	for file in $(ls ${path}/02_OTU/taxa_heatmap/cluster_group/*.csv)
	do
		Rscript /home/jbwang/code/otu/heatmap_group.R ${path}/sample_info.txt ${file}
	done

	for file in $(ls ${path}/02_OTU/taxa_heatmap/cluster/*.csv)
	do
		Rscript /home/jbwang/code/otu/heatmap_group.R ${path}/sample_info.txt ${file}
	done
}

venn(){
	Rscript /home/jbwang/code/otu/venn.R $path
}

barplot(){
	for file in $(ls ${path}/02_OTU/Top10*/*/*.csv)
	do
		Rscript /home/jbwang/code/otu/top10_barplot.R ${file}
	done
}

krona_tree(){
	python3 /home/jbwang/code/generateFile.py $path
	file_list=$(ls ${path}/02_OTU/krona/file_list/*.txt | sort)
	ktImportText ${file_list} -o ${path}/02_OTU/krona/krona.html &
	python3 /home/jbwang/code/otu/level_tree_sample.py -f ${path}/02_OTU/otu_table_tax.txt -t 10 -o ${path}/02_OTU/taxa_tree/sample_tree &
	python3 /home/jbwang/code/otu/level_tree.py -f ${path}/02_OTU/otu_table_tax.txt -t 10 -g ${path}/02_OTU/taxa_tree/group.txt -o ${path}/02_OTU/taxa_tree &
	wait
}

alpha(){
	Rscript /home/jbwang/code/alpha/AlphaDiversity.R ${path} > /dev/null &
	Rscript /home/jbwang/code/alpha/rareplot2.R ${path} &
	wait
}

beta(){
	pcoa=${path}/03_Diversity/Beta/PCoA
	tree_barplot=${path}/03_Diversity/Beta/Tree_barplot
	mkdir -p $pcoa ${tree_barplot}/weighted_unifrac ${tree_barplot}/unweighted_unifrac
	cp ${path}/03_Diversity/Beta/weighted_unifrac_dm.txt ${tree_barplot}/weighted_unifrac
	cp ${path}/03_Diversity/Beta/unweighted_unifrac_dm.txt ${tree_barplot}/unweighted_unifrac
	Rscript /home/jbwang/code/beta/BetaDiversity.R ${path} >/dev/null

	for file in $(ls ${path}/03_Diversity/Beta/Tree_barplot/*/*_dm.txt)
	do
		Rscript /home/jbwang/code/beta/upgma.R ${path} ${file}
	done
}

main(){
	heatmap &
	venn &
	barplot &
	krona_tree &
	alpha &
	beta &
	wait
}

# 关闭qiime1环境
# conda deactivate
source deactivate qiime1
main
