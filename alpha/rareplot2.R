library(vegan)    #���ڼ��� Shannon ��ָ����Simpson ָ����Chao1 ָ����ACE ָ���ȣ�ͬʱ���ڳ���
library(picante)    #���ڼ��� PD_whole_tree��������������������ء���ʵ�ϣ�picante ������ʱĬ��ͬʱ���� vegan
library(ggplot2)
library(doBy)
# library(ggalt)

args <- commandArgs(T)
path <- as.character(args[1])

otu <- read.table(paste0(path, '/02_OTU/otu_table.txt'),header = T, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)

##���庯��
#������� Alpha ������ָ�����������������
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
    if (method == 'richness') result <- rowSums(x > 0)    #�ḻ��ָ��
    else if (method == 'chao1') result <- estimateR(x)[3, ]    #Chao1 ָ��
    else if (method == 'ace') result <- estimateR(x)[5, ]    #ACE ָ��
    else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon ָ��
    else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson ָ��
    else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou ���ȶ�
    else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
    else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
        pd <- pd(x, tree, include.root = FALSE)
        result <- pd[ ,1]
        names(result) <- rownames(pd)
    }
    return(result)
}


#���ݳ���������step����ͳ��ÿ��ϡ���ݶ��µ� Alpha ������ָ��������������б�
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



# #ͳ�� OTU ��ȱ��и������� Shannon ָ������������ʹ�� e
# alpha_index <- function(otu, method = 'shannon', base = exp(1))
# #�� 1000 ������Ϊ�������������ζ� OTU ��ϡ�ͳ�����ֱ�����������ȣ���ͳ�Ƹ������ݶ��µ� OTU ��ȱ��и������� Shannon ָ������������ʹ�� e
# alpha_curves <- function(otu, step = 1000, method = 'shannon', base = exp(1))
    
  
    
    
#���������ַḻ��ָ��Ϊ������ Alpha ���������ߣ���Ϊ�ḻ��ָ��ʱ����һ�����Ƽ�Ϊ��˵��ϡ�����ߣ��������ۼ����ߣ�
#�� 2000 ������step=2000��Ϊ��ͳ��
curves <- alpha_curves(otu, step = 2000, method = 'richness')
##Shannon ��ʽ�Ķ�������Ĭ��Ϊ e��������Ҫ�ɸ��ģ����� 2��
# curves <- alpha_curves(otu, step = 2000, method = 'shannon', base = 2)

#��� ggplot2 ��ͼ�ļ�
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

# ##����㼸���Ի�ȡ��ֵ �� ��׼�Ȼ����չʾ��Ҳ��һ��������ѡ��
# #�ظ����� 5 ��
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
# #�����ֵ �� ��׼�doBy ���е� summaryBy() ������
# plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
# plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
# plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA
# 
# #ggplot2 ��ͼ
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

#ggplot2 ��ͼ��ʹ�õ� ggalt ���� geom_xspline() ����ƽ������ߣ�

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
