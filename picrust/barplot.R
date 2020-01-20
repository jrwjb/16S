library(reshape2)
library(gtools)
library(ggplot2)
library(RColorBrewer)

args <- commandArgs(T)
path <- as.character(args[1])
file <- as.character(args[2])
top <- as.numeric(args[3])
# setwd(path)
# path <- sapply(strsplit(file, '\\.'), '[', 1)
# data <- read.csv(file = file, header = T, row.names = 1)

data <- read.table(file, quote = "", row.names = 1, header = T, sep = '\t', check.names = F)
data <- data[mixedsort(colnames(data), decreasing = T)] # 字母数字混合排序

## kegg
data['Total'] <- rowSums(data)
data <- data[order(data['Total'], decreasing = T),]

data <- data[-ncol(data)]
data_top <- data[1:top,]
data_other <- data[top + 1:nrow(data),]
data_other[is.na(data_other)] <- 0
data_top_other <- rbind(data_top, colSums(data_other))
rownames(data_top_other)[top + 1] <- 'Others'
# data <- data[,-ncol(data)]
data_top_other <- t(t(data_top_other)/colSums(data_top_other))  # 相对值
data_top_other <- as.matrix(t(data_top_other))
# data <- data[order(data[,1],decreasing = F),]
newdata <- melt(data_top_other)
newdata$Var1 <- as.character(newdata$Var1)

# data <- t(t(data)/colSums(data))
# data <- as.matrix(t(data))
# newdata <- melt(data)

##设置颜色
mycol<-c(brewer.pal(9, "Set1"),brewer.pal(12, "Paired"),brewer.pal(12, "Set3"))

ggplot(newdata, aes(x=Var1, y=value, fill=Var2)) + 
    geom_bar(stat = "identity", position = 'fill') +
    theme_bw() +
    theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
    # # Remove x y axis title
    # theme(axis.title.x = element_blank()) +
    # theme(axis.title.y = element_blank()) +
    labs(x = "Sample Name", y = "Relative abundance", title = "") +
    theme(legend.title=element_blank(), legend.text = element_text(size = 6.5)) +
    scale_fill_manual(values = mycol) +
    theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 45))

ggsave(paste0(path, sapply(strsplit(file, '[./]'), '[', 2), '.png'), width = 28, height = 16, units ="cm", bg="white", dpi=600)
ggsave(paste0(path, sapply(strsplit(file, '[./]'), '[', 2), '.pdf'), width = 12, height = 7)

rm(list = ls())
