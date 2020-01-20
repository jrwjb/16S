library(pheatmap)
library(gtools)

args <- commandArgs(T)
path <- as.character(args[1])
file <- as.character(args[2])
# file <- 'C:\\Users\\jbwang\\Desktop\\test\\02_OTU\\taxa_heatmap\\cluster\\Class.csv'
# path <- sapply(strsplit(file, '\\.'), '[', 1)

data <- read.table(file,quote = "", row.names = 1, header = T, sep = '\t', check.names = F)
data <- data[mixedsort(colnames(data), decreasing = T)] # 字母数字混合排序
data <- t(t(data)/colSums(data))  # 相对值


# if (row_height < 15){
#   pdf(file = paste0(path, '.pdf'), width = 15, height = 10, family="GB1")                     
#   png(filename = paste0(path, '.png'), width = 30, height = 18, units = "cm",res = 600)   #row_height/3  
# }else if (row_height < 30){
#   pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/3, family="GB1")             
#   png(filename = paste0(path, '.png'), width = 30, height = row_height/1.5, units = "cm",res = 600)
# }else if (row_height < 45){
#   pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/4.5, family="GB1")             
#   png(filename = paste0(path, '.png'), width = 30, height = row_height/2.25, units = "cm",res = 600)
# }else if (row_height < 60){
#   pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/5, family="GB1")           
#   png(filename = paste0(path, '.png'), width = 30, height = row_height/2.5, units = "cm",res = 600)
# }else if (row_height > 800){
#   pdf(file = paste0(path, '.pdf'), width = 40, height = row_height/6, family="GB1")                                 #row_height/6  width = 10
#   png(filename = paste0(path, '.png'), width = 50, height = row_height/3, units = "cm",res = 200)   #row_height/3  width = 20
# }else{
#   pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/6, family="GB1")                                 #row_height/6  width = 10
#   png(filename = paste0(path, '.png'), width = 30, height = row_height/3, units = "cm",res = 300)   #row_height/3  width = 20
# }

pdf(file = paste0(path, sapply(strsplit(file, '[./]'), '[', 2), '.pdf'), width = 10, height = 8, family="GB1")
png(paste0(path, sapply(strsplit(file, '[./]'), '[', 2), '.png'), width = 30, height = 25, units = "cm",res = 300)

p <- pheatmap(as.matrix(data), scale="row", cluster_rows=F, cluster_cols=F)

p
dev.off()
p
dev.off()

rm(list = ls())








