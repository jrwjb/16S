library(phangorn)
library(ggtree)
library(vegan)
library(phytools)
library(RColorBrewer)

args <- commandArgs(T)
path <- as.character(args[1])
file <- as.character(args[2])
out <- sapply(strsplit(file, '\\.'), '[', 1)

## 距离矩阵（由qiime软件生成)
dm <- read.table(file, row.names = 1, header = T, sep = '\t')

###方法一
#### 使用距离矩阵生成树文件并作图
# h<-hclust(as.dist(dm), method='average')
# png(filename = "UPGMA.png",width = 26,height = 26,units = "cm",res = 300)
# plot(as.dendrogram(h), main = 'Sample Cluster',horiz=T)
# dev.off()

###方法二
###使用phangorn中的upgma函数生成树文件并作图
tree <- upgma(as.matrix(dm))
write.tree(tree, file = paste0(out, '.tre'))
png(filename = paste0(out, "_upgma.png"), width = 26, height = 26, units = "cm", res = 300)
plot(tree, main = '')
dev.off()

###方法三
##### 直接调用qiime脚本生成的树文件(upgma_cluster.py -i weighted_unifrac_sorted_otu_table.txt -o qiime_weighted_upgma.tre)
# qiime_tre <- read.tree('qiime_weighted_upgma.tre')
# png(filename = "UPGMA-3.png",width = 26,height = 26,units = "cm",res = 300)
# plot(qiime_tre)
# dev.off()

####tree_barplot
data <- read.csv(paste0(path, '/02_OTU/Top10/Phylum/Phylum_top10.csv'), header = T, row.names = 1)
# data <- t(data)
mycol <- c(brewer.pal(9, "Set1"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))
png(filename = paste0(out, "_tree_barplot.png"), width = 25, height = 22, units = "cm", res = 300)

plotTree.barplot(tree, data,
                 # args.plotTree=list(ftype="off"),
                 args.plotTree=list(fsize=1, ftype="reg", lwd=1.5, xlim=c(-0.2,0.6)),
                 args.barplot=list(col=mycol, border=mycol, xlab="",
                                   xlim=c(0.2, 1.5)), args.axis=list(at=seq(0, 1, by=0.2)))
legend(x="right", legend=names(data), pch=22, pt.cex=1.5, pt.bg=mycol, cex = 0.7, box.col="transparent")

mtext("Relative abundance in Phylum level", 1, at=0.5, line=2.5)
dev.off()

rm(list = ls())
