library(ggplot2)

# 输入参数1：路径
args <- commandArgs(T)
path <- as.character(args[1])
# path <- 'C:\\Users\\jbwang\\Desktop\\1'
# setwd(path)

############ Alpha多样性 #############
##### 稀释曲线 #########
# 处理后的数据observed_species_T.txt，原数据observed_species.txt由qiime1生成(其他曲线如shannon曲线均可由此生成)
# data1 <- read.table(paste0(path, '/Diversity/Alpha/alpha_div_collated/observed_species_T.txt'), header = T, sep = '\t')
file <- as.character(args[2])
out <- sapply(strsplit(file, '\\.'), '[', 1)
out <- substr(out, start = 1, stop = nchar(out) - 2)
index <- sapply(strsplit(out, '/'), '[', length(unlist(strsplit(out, '/'))))
data <- read.table(file, header = T, sep = '\t')

alpha_func <- function(index){
    # 按sampleID画图
    p1 <- ggplot(data,aes(x = Seqs, y = data[, c(index)], color = sampleID)) +
        geom_smooth(se = F, method = 'lm', formula = y ~ log(x)) +
        theme_bw() +
        labs(x = "Seqs Num", y = paste0(index, ' Num'), title = "") +
        theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
        theme(legend.text = element_text(size = 6))
    ggsave(paste0(out, '_sample.png'), width = 25, height = 18, units ="cm", dpi = 600)
    ggsave(paste0(out, '_sample.pdf'), width = 12, height = 7)

    # 按组画图
    p2 <- ggplot(data, aes(x = Seqs, y = data[, c(index)], color = Group, group = Group)) +
        geom_smooth(se = F, method = 'lm', formula = y ~ log(x)) +
        theme_bw() +
        labs(x = "Seqs Num", y = paste0(index, ' Num'), title = "") +
        theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
        theme(legend.text = element_text(size = 6))
    ggsave(paste(out, '_group.png', sep=""), width = 25, height = 18, units ="cm", dpi = 600)
    ggsave(paste(out, '_group.pdf', sep=""))
}

alpha_func(index)


rm(list = ls())
