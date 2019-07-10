#!/bin/bash


project=$1

if [ -d "${project}/Result" ]; then
	rm -rf ${project}/Result
fi

##### QC #####
mkdir -p ${project}/Result/01_QC
cp ${project}/QcStatic.csv ${project}/Result/01_QC

## ##### OTU ####
mkdir -p ${project}/Result/02_OTU_Taxa
cp -r ${project}/02_OTU/Venn ${project}/Result/02_OTU_Taxa
rm -rf ${project}/Result/02_OTU_Taxa/Venn/*.log
cp ${project}/02_OTU/otu_table_tax.txt ${project}/Result/02_OTU_Taxa
sed -i 's/\t/,/g' ${project}/Result/02_OTU_Taxa/otu_table_tax.txt
mv ${project}/Result/02_OTU_Taxa/otu_table_tax.txt ${project}/Result/02_OTU_Taxa/otu_table_tax.csv

# ##### community ###
mkdir -p ${project}/Result/03_Community/community
cp -r ${project}/02_OTU/Top10* ${project}/Result/03_Community/community
mv ${project}/Result/03_Community/community/Top10 ${project}/Result/03_Community/community/Top10_sample
mkdir -p ${project}/Result/03_Community/krona
cp -r ${project}/02_OTU/krona/*.html ${project}/Result/03_Community/krona
cp -r ${project}/02_OTU/taxa_heatmap ${project}/Result/03_Community/
mv ${project}/Result/03_Community/taxa_heatmap/cluster ${project}/Result/03_Community/taxa_heatmap/cluster_sample
mkdir -p ${project}/Result/03_Community/taxa_tree/group_tree
cp ${project}/02_OTU/taxa_tree/level_tree.* ${project}/02_OTU/taxa_tree/tax_percent* ${project}/Result/03_Community/taxa_tree/group_tree
cp -r ${project}/02_OTU/taxa_tree/sample_tree ${project}/Result/03_Community/taxa_tree
sed -i 's/\t/,/g' ${project}/Result/03_Community/taxa_tree/sample_tree/*.txt
rename 's/.txt/.csv/' ${project}/Result/03_Community/taxa_tree/sample_tree/*.txt
cp -r ${project}/03_Diversity/Beta/Tree_barplot ${project}/Result/03_Community
mv ${project}/Result/03_Community/Tree_barplot ${project}/Result/03_Community/UPGMA
rm -rf ${project}/Result/03_Community/UPGMA/*/*.txt ${project}/Result/03_Community/UPGMA/*/*.tre

# #### alpha ###
mkdir -p ${project}/Result/04_Alpha_diversity/Rank_Abundance
cp ${project}/03_Diversity/Alpha/rank_sampleID.* ${project}/Result/04_Alpha_diversity/Rank_Abundance
cp -r ${project}/03_Diversity/Alpha/Specaccum ${project}/Result/04_Alpha_diversity
mkdir ${project}/Result/04_Alpha_diversity/Shannon
mkdir ${project}/Result/04_Alpha_diversity/Rarefaction_Curve
cp ${project}/03_Diversity/Alpha/alpha_div_collated/observed_species_sample.* ${project}/Result/04_Alpha_diversity/Rarefaction_Curve
cp ${project}/03_Diversity/Alpha/alpha_div_collated/shannon_sample.* ${project}/Result/04_Alpha_diversity/Shannon
cp -r ${project}/03_Diversity/Alpha/alpha_div_collated ${project}/Result/04_Alpha_diversity
mv ${project}/Result/04_Alpha_diversity/alpha_div_collated ${project}/Result/04_Alpha_diversity/Alpha_div_diff
rm -rf ${project}/Result/04_Alpha_diversity/Alpha_div_diff/*.png ${project}/Result/04_Alpha_diversity/Alpha_div_diff/*.pdf ${project}/Result/04_Alpha_diversity/Alpha_div_diff/*.txt
cp ${project}/03_Diversity/Alpha/alpha_diversity_index.txt ${project}/Result/04_Alpha_diversity
sed -i 's/\t/,/g' ${project}/Result/04_Alpha_diversity/alpha_diversity_index.txt
mv ${project}/Result/04_Alpha_diversity/alpha_diversity_index.txt ${project}/Result/04_Alpha_diversity/alpha_diversity_index.csv

# #### beta ###
mkdir -p ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/anosim ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/adonis ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/PCA ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/PCoA ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/NMDS ${project}/Result/05_Beta_diversity
cp -r ${project}/03_Diversity/Beta/Beta_div ${project}/Result/05_Beta_diversity/Beta_div_diff

#  ### diff_analysis ###
mkdir -p ${project}/Result/06_Differential_analysis/LefSe
cp -r ${project}/03_Diversity/Beta/lefse/*_vs_* ${project}/Result/06_Differential_analysis/LefSe
rm -rf ${project}/Result/06_Differential_analysis/LefSe/*/*.txt ${project}/Result/06_Differential_analysis/LefSe/*.in
