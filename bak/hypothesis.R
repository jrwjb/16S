library(vegan)
library(RColorBrewer)

args <- commandArgs(T)
path <- as.character(args[1])

sample_info <- read.table(paste0(path, '/sample_info.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- read.table(paste0(path, '/03_Diversity/otu_table_even.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- otu_table[,rownames(sample_info)]
dm <- vegdist(t(otu_table), method = 'bray')

## anosim可视化
dir.create(paste0(path, '/03_Diversity/Beta/anosim'))
anosim_result<-anosim(dm, sample_info$Group, permutations = 999)
group_num <- length(unique(sample_info$Group))
cols<-brewer.pal(n=group_num + 1, name="Set1")
p_value <- anosim_result$signif
group <- unique(sample_info$Group)
group <- paste0(group[1], '_vs_', group[2])
anosim_df <- data.frame(group, p_value)
write.table(anosim_df, paste0(path, '/03_Diversity/Beta/anosim/anosim_result.txt'), sep = '\t',row.names = F)
png(filename = paste0(path, "/03_Diversity/Beta/anosim/anosim_result.png"), width = 26, height = 22, units = "cm", res = 300)
plot(anosim_result, col=cols)
dev.off()
#pdf
pdf(filename = paste0(path, "/03_Diversity/Beta/anosim/anosim_result.pdf"))
plot(anosim_result, col=cols)
dev.off()

## adonis
dir.create(paste0(path, '/03_Diversity/Beta/adonis'))
adonis_result <- adonis(dm~sample_info$Group, permutations = 999)
adonis_result <- adonis_result$aov.tab
write.table(adonis_result, paste0(path, '/03_Diversity/Beta/adonis/adonis_result.txt'), sep = '\t')


rm(list = ls())
