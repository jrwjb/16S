library(reshape2)
library(ggplot2)
library(RColorBrewer)

args <- commandArgs(T)
file <- as.character(args[1])
path <- sapply(strsplit(file, '\\.'), '[', 1)
data <- read.csv(file = file, header = T, row.names = 1)
data <- as.matrix(t(data))
data <- data[order(data[,1],decreasing =F),]
data <- melt(data)
data$Var2 <- as.character(data$Var2)
##设置颜色
mycol<-c(brewer.pal(9, "Set1"),brewer.pal(12, "Paired"),brewer.pal(12, "Set3"))

p <- ggplot(data, aes(x=Var2, y=value, fill=Var1)) + 
    geom_bar(stat = "identity", position = 'fill') +
    theme_bw() +
    theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
    # # Remove x y axis title
    # theme(axis.title.x = element_blank()) +
    # theme(axis.title.y = element_blank()) +
    labs(x = "Sample Name", y = "Relative abundance", title = "") +
    theme(legend.title=element_blank()) +
    scale_fill_manual(values = mycol) +
    theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 45))

ggsave(paste0(path, '.png'), width = 28, height = 16, units ="cm", bg="white", dpi=600)
ggsave(paste0(path, '.pdf'), width = 12, height = 7)

rm(list = ls())





