#!/usr/bin/env bash

set -e;  # 出现错误立即退出

project=$1
project=${project%/}     #  '%' 从后向前删除, '#' 从前向后删除

source activate qiime1
ori_path=$(pwd)
mkdir -p ${project}/04_Picrust
cp ${project}/*.txt ${project}/04_Picrust

# 默认按照gg138版本pick otu
predict(){
	cd ${project}/04_Picrust
	pick_closed_reference_otus.py -i ../01_CleanData/all.fa.bak -o out -f && \
	biom convert -i out/otu_table.biom -o out/otu_table_tax.txt --table-type="OTU table" --header-key taxonomy --to-tsv
	# 标准化otu_table
	normalize_by_copy_number.py -i out/otu_table.biom -o out/otu_table_normalized.biom
	# 预测(默认ko)
	mkdir -p KEGG
	predict_metagenomes.py -i out/otu_table_normalized.biom -o ./KEGG/ko_prediction.biom
	# 统计不同水平
	# kegg的3个水平, -c 指输出类型,有KEGG_Pathways, COG_Category, RFAM三种,-l是级别
	categorize_by_function.py -i KEGG/ko_prediction.biom -o ./KEGG/kegg_predicted_L3.biom -c KEGG_Pathways -l 3
	categorize_by_function.py -i KEGG/ko_prediction.biom -o ./KEGG/kegg_predicted_L2.biom -c KEGG_Pathways -l 2
	categorize_by_function.py -i KEGG/ko_prediction.biom -o ./KEGG/kegg_predicted_L1.biom -c KEGG_Pathways -l 1 

	## 预测COG
	mkdir -p COG
	predict_metagenomes.py --type_of_prediction cog -i out/otu_table_normalized.biom -o ./COG/cog_prediction.biom

	categorize_by_function.py -i COG/cog_prediction.biom -o ./COG/cog_predicted_L2.biom -c COG_Category -l 2
	categorize_by_function.py -i COG/cog_prediction.biom -o ./COG/cog_predicted_L1.biom -c COG_Category -l 1
	cd ori_path
}

if [ ! -e ${project}/predict.success ];then
	predict
	if [ "$?" -ne 0 ]; then
		printf "\033[31m Error: function predict is failed and exit, please try it again !!! \033[0m\n"
		exit 1
	else
		touch ${project}/predict.success
	fi
fi

cd ${project}/04_Picrust
#biom=${i##*/}    #返回 / 后的字符
for i in $(ls ./*/*.biom)
do
	sh /home/jbwang/code/picrust/biom2txt.sh $i
done

sed -i '/# Const/d;s/#//;1s/\./_/g' ./*/*.txt
#summarize_taxa.py -i ./kegg/kegg_predicted_L1.biom -o level1 --md_identifier "KEGG_Pathways" --level 1

mkdir -p KEGG/PCA KEGG/Barplot KEGG/Heatmap KEGG/STAMP KEGG/LEfSe COG/PCA COG/Barplot COG/Heatmap COG/STAMP COG/LEfSe

# pca
Rscript /home/jbwang/code/picrust/pca_sw.R KEGG/PCA/ KEGG/ko_prediction.txt && \
Rscript /home/jbwang/code/picrust/pca_sw.R COG/PCA/ COG/cog_prediction.txt

# barplot
Rscript /home/jbwang/code/picrust/barplot.R KEGG/Barplot/ KEGG/kegg_predicted_L2.txt 10 && \
Rscript /home/jbwang/code/picrust/barplot.R COG/Barplot/ COG/cog_predicted_L2.txt 10

# heatmap
Rscript /home/jbwang/code/picrust/heatmap.R KEGG/Heatmap/ KEGG/kegg_predicted_L2.txt && \
Rscript /home/jbwang/code/picrust/heatmap.R COG/Heatmap/ COG/cog_predicted_L2.txt

## lefse-kegg、cog
##kegg_predicted_L3_lefse.txt 由 kegg_predicted_L3.txt转化过来
python3 /home/jbwang/code/picrust/lefse-data.py KEGG/kegg_predicted_L2.txt sample_info.txt KEGG/LEfSe/kegg_predicted_L2_lefse.txt && \
python3 /home/jbwang/code/picrust/lefse-data.py KEGG/kegg_predicted_L3.txt sample_info.txt KEGG/LEfSe/kegg_predicted_L3_lefse.txt && \
python3 /home/jbwang/code/picrust/lefse-data.py COG/cog_predicted_L2.txt sample_info.txt COG/LEfSe/cog_predicted_L2_lefse.txt

python3 /home/jbwang/code/picrust/lefse-function_prediction.py ../
