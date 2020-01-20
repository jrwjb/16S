library(pheatmap)

args <- commandArgs(T)
sample_info <- as.character(args[1])
# sample_info <- 'C:\\Users\\jbwang\\Desktop\\sample_info.txt'

file <- as.character(args[2])
# file <- 'C:\\Users\\jbwang\\Desktop\\test\\02_OTU\\taxa_heatmap\\cluster\\Class.csv'
path <- sapply(strsplit(file, '\\.'), '[', 1)

sample_info_data <- read.table(sample_info, row.names = 1, header = T, sep = '\t')
group <- as.data.frame(sample_info_data['Group'])
# row.names(group) <- row.names(sample_info_data)

data <- read.csv(file, sep = ",", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
# data <- read.table(file, row.names = 1, header = T, sep = '\t')
tax <- as.data.frame(rownames(data))
colnames(tax) <- 'tax'
row_height <- as.numeric(nrow(data))
row_height

if (row_height < 15){
  pdf(file = paste0(path, '.pdf'), width = 15, height = 10, family="GB1")                     
  png(filename = paste0(path, '.png'), width = 30, height = 18, units = "cm",res = 600)   #row_height/3  
}else if (row_height < 30){
  pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/3, family="GB1")             
  png(filename = paste0(path, '.png'), width = 30, height = row_height/1.5, units = "cm",res = 600)
}else if (row_height < 45){
  pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/4.5, family="GB1")             
  png(filename = paste0(path, '.png'), width = 30, height = row_height/2.25, units = "cm",res = 600)
}else if (row_height < 60){
  pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/5, family="GB1")           
  png(filename = paste0(path, '.png'), width = 30, height = row_height/2.5, units = "cm",res = 600)
}else if (row_height > 800){
  pdf(file = paste0(path, '.pdf'), width = 40, height = row_height/6, family="GB1")                                 #row_height/6  width = 10
  png(filename = paste0(path, '.png'), width = 50, height = row_height/3, units = "cm",res = 200)   #row_height/3  width = 20
}else{
  pdf(file = paste0(path, '.pdf'), width = 15, height = row_height/6, family="GB1")                                 #row_height/6  width = 10
  png(filename = paste0(path, '.png'), width = 30, height = row_height/3, units = "cm",res = 300)   #row_height/3  width = 20
}

# pdf(file = paste0(path, '.pdf'), width = 10, height = 8, family="GB1")
# png(filename = paste0(path, '.png'), width = 25, height = row_height/2, units = "cm",res = 300)

p <- pheatmap(as.matrix(data), scale="row", cluster_rows=T, cluster_cols=FALSE, annotation_col=group)

p
dev.off()
p
dev.off()


rm(list = ls())








