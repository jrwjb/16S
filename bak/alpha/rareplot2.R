library(vegan)    #用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)    #用于计算 PD_whole_tree，若不计算它就无需加载。事实上，picante 包加载时默认同时加载 vegan
library(ggplot2)
library(doBy)
# library(ggalt)

args <- commandArgs(T)
path <- as.character(args[1])

otu <- read.table(paste0(path, '/02_OTU/otu_table.txt'),header = T, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)

##定义函数
#计算多种 Alpha 多样性指数，结果返回至向量
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
    if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
    else if (method == 'chao1') result <- estimateR(x)[3, ]    #Chao1 指数
    else if (method == 'ace') result <- estimateR(x)[5, ]    #ACE 指数
    else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
    else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
    else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
    else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
    else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
        pd <- pd(x, tree, include.root = FALSE)
        result <- pd[ ,1]
        names(result) <- rownames(pd)
    }
    return(result)
}


#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
    x_nrow <- nrow(x)
    if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
    alpha_rare <- list()
    
    for (i in 1:x_nrow) {
        step_num <- seq(0, rare[i], step)
        if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
        
        alpha_rare_i <- NULL
        for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
        names(alpha_rare_i) <- step_num
        alpha_rare <- c(alpha_rare, list(alpha_rare_i))
    }
    
    names(alpha_rare) <- rownames(x)
    return(alpha_rare)
}



# #统计 OTU 丰度表中各样本的 Shannon 指数，对数底数使用 e
# alpha_index <- function(otu, method = 'shannon', base = exp(1))
# #以 1000 条序列为抽样步长，依次对 OTU 表稀释抽样，直到最大序列深度；并统计各抽样梯度下的 OTU 丰度表中各样本的 Shannon 指数，对数底数使用 e
# alpha_curves <- function(otu, step = 1000, method = 'shannon', base = exp(1))
    
  
    
    
#以下以物种丰富度指数为例绘制 Alpha 多样性曲线（当为丰富度指数时，另一个名称即为常说的稀释曲线，或物种累计曲线）
#以 2000 步长（step=2000）为例统计
curves <- alpha_curves(otu, step = 2000, method = 'richness')
##Shannon 公式的对数底数默认为 e，若有需要可更改（例如 2）
# curves <- alpha_curves(otu, step = 2000, method = 'shannon', base = 2)

#获得 ggplot2 作图文件
plot_data <- data.frame()
for (i in names(curves)) {
    curves_i <- (curves[[i]])
    curves_i <- data.frame(rare = names(curves_i), alpha = curves_i, sample = i, stringsAsFactors = FALSE)
    plot_data <- rbind(plot_data, curves_i)
}

rownames(plot_data) <- NULL
plot_data$rare <- as.numeric(plot_data$rare)
plot_data$alpha <- as.numeric(plot_data$alpha)

max_num <- max(rowSums(otu))

ggplot(plot_data, aes(rare, alpha, color = sample)) +
    # geom_line() +
    # geom_xspline() +
    geom_smooth(se = F,method = 'loess', formula = y ~ x) +
    # geom_smooth(se = F, method = 'lm', formula = y ~ log(x)) +
    labs(x = 'Number of sequences', y = 'Observed Species', color = NULL) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
    scale_x_continuous(breaks = seq(0, max_num, 5000), labels = as.character(seq(0, max_num, 5000)))

ggsave(paste0(path, '/03_Diversity/Alpha/observed_species_sample.png'), width = 25, height = 18, units ="cm", dpi = 600)
ggsave(paste0(path, '/03_Diversity/Alpha/observed_species_sample.pdf'), width = 12, height = 7)

# ##多计算几次以获取均值 ± 标准差，然后再展示出也是一个不错的选择
# #重复抽样 5 次
# plot_richness <- data.frame()
# 
# for (n in 1:5) {
#     richness_curves <- alpha_curves(otu, step = 2000, method = 'richness')
#     
#     for (i in names(richness_curves)) {
#         richness_curves_i <- (richness_curves[[i]])
#         richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
#         plot_richness <- rbind(plot_richness, richness_curves_i)
#     }
# }
# 
# #计算均值 ± 标准差（doBy 包中的 summaryBy() 函数）
# plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
# plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
# plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA
# 
# #ggplot2 作图
# ggplot(plot_richness_stat, aes(rare, alpha.mean, color = sample)) +
#     geom_line() +
#     # geom_point() +
#     # geom_errorbar(aes(ymin = alpha.mean - alpha.sd, ymax = alpha.mean + alpha.sd), width = 500) +
#     labs(x = 'Number of sequences', y = 'Richness', color = NULL) +
#     theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
#     geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
#     scale_x_continuous(breaks = seq(0, max_num, 5000), as.character(seq(0, max_num, 5000)))



shannon_curves1 <- alpha_curves(otu, step = 200, rare = 150, method = 'shannon')
shannon_curves2 <- alpha_curves(otu, step = 2000, method = 'shannon')
shannon_curves <- c(shannon_curves1, shannon_curves2)

plot_shannon <- data.frame()
for (i in 1:length(shannon_curves)) {
    shannon_curves_i <- shannon_curves[[i]]
    shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = names(shannon_curves)[i], stringsAsFactors = FALSE)
    plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)
plot_shannon <- plot_shannon[order(plot_shannon$sample, plot_shannon$rare), ]

#ggplot2 作图（使用到 ggalt 包的 geom_xspline() 绘制平滑拟合线）

ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
    geom_line() +
    # geom_xspline() +
    # geom_smooth() +
    labs(x = 'Number of sequences', y = 'Shannon', color = NULL) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
    scale_x_continuous(breaks = seq(0, max_num, 5000), labels = as.character(seq(0, max_num, 5000)))

ggsave(paste0(path, '/03_Diversity/Alpha/shannon_sample.png'), width = 25, height = 18, units ="cm", dpi = 600)
ggsave(paste0(path, '/03_Diversity/Alpha/shannon_sample.pdf'), width = 12, height = 7)

