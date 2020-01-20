#!/bin/bash

project=$1
project=${project%/}
cd $project

if [ -d "Result" ]; then
	rm -rf Result
fi

##### QC #####
mkdir -p Result/01_QC
# cp QcStatic.csv /Result/01_QC

## ##### OTU ####
mkdir -p Result/02_OTU_Taxa/Venn
cp -r 02_OTU/Venn/*_vs_* Result/02_OTU_Taxa/Venn
rm -rf Result/02_OTU_Taxa/Venn/*/*.log
cp 02_OTU/otu_table_tax.txt Result/02_OTU_Taxa
sed -i 's/\t/,/g' Result/02_OTU_Taxa/otu_table_tax.txt
mv Result/02_OTU_Taxa/otu_table_tax.txt Result/02_OTU_Taxa/otu_table_tax.csv

# ##### community ###
mkdir -p Result/03_Community/community
cp -r 02_OTU/Top10* Result/03_Community/community
mv Result/03_Community/community/Top10 Result/03_Community/community/Top10_sample
mkdir -p Result/03_Community/krona
cp -r 02_OTU/krona/*.html Result/03_Community/krona
cp -r 02_OTU/taxa_heatmap Result/03_Community/
mv Result/03_Community/taxa_heatmap/cluster Result/03_Community/taxa_heatmap/cluster_sample
mkdir -p Result/03_Community/taxa_tree/group_tree
cp 02_OTU/taxa_tree/level_tree.* 02_OTU/taxa_tree/tax_percent* Result/03_Community/taxa_tree/group_tree
cp -r 02_OTU/taxa_tree/sample_tree Result/03_Community/taxa_tree
sed -i 's/\t/,/g' Result/03_Community/taxa_tree/sample_tree/*.txt
rename 's/.txt/.csv/' Result/03_Community/taxa_tree/sample_tree/*.txt
cp -r 03_Diversity/Beta/Tree_barplot Result/03_Community
mv Result/03_Community/Tree_barplot Result/03_Community/UPGMA
rm -rf Result/03_Community/UPGMA/*/*.txt Result/03_Community/UPGMA/*/*.tre

# #### alpha ###
mkdir -p Result/04_Alpha_diversity/Rank_Abundance
cp 03_Diversity/Alpha/rank_sampleID.* Result/04_Alpha_diversity/Rank_Abundance
cp -r 03_Diversity/Alpha/Specaccum Result/04_Alpha_diversity
mkdir Result/04_Alpha_diversity/Shannon
mkdir Result/04_Alpha_diversity/Rarefaction_Curve
cp 03_Diversity/Alpha/observed_species_sample.* Result/04_Alpha_diversity/Rarefaction_Curve
cp 03_Diversity/Alpha/shannon_sample.* Result/04_Alpha_diversity/Shannon
cp -r 03_Diversity/Alpha/alpha_div_collated Result/04_Alpha_diversity
mv Result/04_Alpha_diversity/alpha_div_collated Result/04_Alpha_diversity/Alpha_div_diff
rm -rf Result/04_Alpha_diversity/Alpha_div_diff/*.png Result/04_Alpha_diversity/Alpha_div_diff/*.pdf Result/04_Alpha_diversity/Alpha_div_diff/*.txt
cp 03_Diversity/Alpha/alpha_diversity_index.txt Result/04_Alpha_diversity
sed -i 's/\t/,/g' Result/04_Alpha_diversity/alpha_diversity_index.txt
mv Result/04_Alpha_diversity/alpha_diversity_index.txt Result/04_Alpha_diversity/alpha_diversity_index.csv

# #### beta ###
mkdir -p Result/05_Beta_diversity
cp -r 03_Diversity/Beta/anosim Result/05_Beta_diversity
cp -r 03_Diversity/Beta/adonis Result/05_Beta_diversity
cp -r 03_Diversity/Beta/PCA Result/05_Beta_diversity
cp -r 03_Diversity/Beta/PCoA Result/05_Beta_diversity
cp -r 03_Diversity/Beta/NMDS Result/05_Beta_diversity
cp -r 03_Diversity/Beta/Beta_div Result/05_Beta_diversity/Beta_div_diff
rm -rf Result/05_Beta_diversity/Beta_div_diff/*/Group* Result/05_Beta_diversity/*/all

#  ### diff_analysis ###
mkdir -p Result/06_Differential_analysis/LefSe Result/06_Differential_analysis/STAMP
cp -r 03_Diversity/Beta/lefse/*_vs_* Result/06_Differential_analysis/LefSe
rm -rf Result/06_Differential_analysis/LefSe/*/*.txt Result/06_Differential_analysis/LefSe/*/*.in

cp -r 03_Diversity/Beta/PCA/*_vs_* Result/06_Differential_analysis/STAMP
rm -rf Result/06_Differential_analysis/STAMP/*/*
cp 02_OTU/Abundance/Relative/Genus.csv Result/06_Differential_analysis/STAMP

#   function prediction
mkdir -p Result/07_FunctionPrediction
cp -r 04_Picrust/KEGG Result/07_FunctionPrediction
cp -r 04_Picrust/COG Result/07_FunctionPrediction
cp -r 03_Diversity/Beta/PCA/*_vs_* Result/07_FunctionPrediction/KEGG/STAMP
cp -r 03_Diversity/Beta/PCA/*_vs_* Result/07_FunctionPrediction/COG/STAMP
rm -rf Result/07_FunctionPrediction/*/*.biom Result/07_FunctionPrediction/*/*/all Result/07_FunctionPrediction/*/*/*.txt Result/07_FunctionPrediction/*/*/*/L3* Result/07_FunctionPrediction/*/*/*/L2.txt Result/07_FunctionPrediction/*/*/*/L2.in Result/07_FunctionPrediction/KEGG/STAMP/*/* Result/07_FunctionPrediction/COG/STAMP/*/*
sed -i 's/,/;/g;s/\t/,/g' Result/07_FunctionPrediction/*/*.txt
rename 's/.txt/.csv/' Result/07_FunctionPrediction/*/*.txt
