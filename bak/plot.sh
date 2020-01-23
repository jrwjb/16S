#!/usr/bin/env bash

set -e;  # 出现错误立即退出

project=$1
project=${project%/}     #  '%' 从后向前删除, '#' 从前向后删除
[ $project = "." ] && project=$(pwd)

heatmap(){
	for file in $(ls ${project}/02_OTU/taxa_heatmap/cluster_group/*.csv)
	do
		Rscript /home/jbwang/code/otu/heatmap_group.R ${project}/sample_info.txt ${file}
	done

	for file in $(ls ${project}/02_OTU/taxa_heatmap/cluster/*.csv)
	do
		Rscript /home/jbwang/code/otu/heatmap_group.R ${project}/sample_info.txt ${file}
	done
}

venn(){
	Rscript /home/jbwang/code/otu/venn.R $project
}

barplot(){
	for file in $(ls ${project}/02_OTU/Top10*/*/*.csv)
	do
		Rscript /home/jbwang/code/otu/top10_barplot.R ${file}
	done
}

krona_tree(){
	python3 /home/jbwang/code/generateFile.py $project
	file_list=$(ls ${project}/02_OTU/krona/file_list/*.txt | sort)
	ktImportText ${file_list} -o ${project}/02_OTU/krona/krona.html &
	python3 /home/jbwang/code/otu/level_tree_sample.py -f ${project}/02_OTU/otu_table_tax.txt -t 10 -o ${project}/02_OTU/taxa_tree/sample_tree &
	python3 /home/jbwang/code/otu/level_tree.py -f ${project}/02_OTU/otu_table_tax.txt -t 10 -g ${project}/02_OTU/taxa_tree/group.txt -o ${project}/02_OTU/taxa_tree &
	wait
}

alpha(){
	Rscript /home/jbwang/code/alpha/AlphaDiversity.R ${project} > /dev/null &
	Rscript /home/jbwang/code/alpha/rareplot2.R ${project} &
	wait
}

beta(){
	pcoa=${project}/03_Diversity/Beta/PCoA
	tree_barplot=${project}/03_Diversity/Beta/Tree_barplot
	mkdir -p $pcoa ${tree_barplot}/weighted_unifrac ${tree_barplot}/unweighted_unifrac
	cp ${project}/03_Diversity/Beta/weighted_unifrac_dm.txt ${tree_barplot}/weighted_unifrac
	cp ${project}/03_Diversity/Beta/unweighted_unifrac_dm.txt ${tree_barplot}/unweighted_unifrac
	Rscript /home/jbwang/code/beta/BetaDiversity.R ${project} >/dev/null

	for file in $(ls ${project}/03_Diversity/Beta/Tree_barplot/*/*_dm.txt)
	do
		Rscript /home/jbwang/code/beta/upgma.R ${project} ${file}
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
