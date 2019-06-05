library(vegan)
library(ggplot2)
library(tidyverse)
library(agricolae)
# library(ggthemes)
# library(ggsci)

# 输入参数1：路径
args <- commandArgs(T)
path <- as.character(args[1])
# setwd(path)

# 读取数据
sample_info <- read.table(paste0(path, '/sample_info.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- read.table(paste0(path, '/03_Diversity/otu_table_even.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- otu_table[,rownames(sample_info)]
alpha_diversity_index <- read.table(paste0(path, '/03_Diversity/Alpha/alpha_diversity_index.txt'), header = T, row.names = 1, sep = '\t')

############ Alpha多样性 #############
#### 物种累积曲线 #########
dir.create(paste0(path, '/03_Diversity/Alpha/Specaccum'))
otu_table_T <- t(otu_table)
sp <- specaccum(otu_table_T, method = 'random')
# 保存pdf
pdf(paste0(path, '/03_Diversity/Alpha/Specaccum/specaccum.pdf'))
plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'lightblue', xlab = 'Number of samples', ylab = 'OTUs')
boxplot(sp, col = 'yellow', add = TRUE, pch = '+')
dev.off()
# 保存png
png(paste(path, '/03_Diversity/Alpha/Specaccum/specaccum.png', sep=""), width = 25, height = 18,units ="cm",bg="white",res=300)
plot(sp, ci.type = 'poly', col = 'blue', lwd = 2, ci.lty = 0, ci.col = 'lightblue', xlab = 'Number of samples', ylab = 'OTUs')
boxplot(sp, col = 'yellow', add = TRUE, pch = '+')
dev.off()

# ##### rank abundancd曲线 #########
# ####  构造函数
rank_abundance_func <- function(x1, x2) {
    rank_func <- function(x){
        res <- length(x) - rank(x, ties.method = 'first') + 1
        res[x == 0] <- NA
        return(res)
    }
    
    # 生成rank矩阵
    otu_tab_rank <- map_dfr(x1, rank_func)
    otu_tab_rank <- as.data.frame(otu_tab_rank)
    rownames(otu_tab_rank) <- rownames(x1)
    
    # 重朔数据
    otu_tab <- x1 %>% rownames_to_column('FeatureID') %>% gather(x2,'Count',-FeatureID)
    otu_tab_rank2 <- otu_tab_rank %>% rownames_to_column('FeatureID') %>% gather(x2, 'Rank', -FeatureID)

    # 合并otu_tab数据和rank数据
    data <- left_join(otu_tab, otu_tab_rank2) %>% group_by(x2) %>% mutate(Abundance = Count / sum(Count)) %>% na.omit()

}

# # 按sampleID画图
data <- rank_abundance_func(otu_table, 'SampleID')
colnames(data)[2] <- 'SampleID'
p <- ggplot(data, aes(x = Rank, y = Abundance, color = SampleID)) +
            geom_line() +
            scale_y_log10() +
            labs(title = '') +
            ylab('Abundance') +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(lineheight = 2.5, face = 'bold', hjust = 0.5))

ggsave(paste(path, '/03_Diversity/Alpha/rank_sampleID.png', sep=""), width = 25, height = 18, units ="cm", dpi = 600)
ggsave(paste(path, '/03_Diversity/Alpha/rank_sampleID.pdf', sep=""))

# # 按组画图
group <- unique(sample_info$Group)
total_df <- c(1:nrow(otu_table))
for (i in 1:length(group)){
    tmp <- otu_table[, grep(group[i], colnames(otu_table))]
    tmp$total <- rowSums(tmp)
    total_df <- cbind(total_df, tmp$total)
}
total_df <- total_df[,-1]
colnames(total_df) <- group
rownames(total_df) <- rownames(otu_table)

data2 <- rank_abundance_func(as.data.frame(total_df), 'Group')
colnames(data2)[2] <- 'Group'
p <- ggplot(data2, aes(x = Rank, y = Abundance, color = Group)) +
            geom_line() +
            scale_y_log10() +
            labs(title = '') +
            ylab('Abundance') +
            theme_bw() +
            theme(panel.grid = element_blank(), plot.title = element_text(lineheight = 2.5, face = 'bold', hjust = 0.5))
ggsave(paste(path, '/03_Diversity/Alpha/rank_Group.png', sep=""), width = 25, height = 18, units ="cm", dpi = 600)
ggsave(paste(path, '/03_Diversity/Alpha/rank_Group.pdf', sep=""))     

##### 多样性指数箱线图  #########
index_list <- c('shannon', 'simpson', 'ace', 'goods_coverage', 'chao1', 'observed_species', 'PD_whole_tree')

box_plot_func <- function(index) {
    data <- cbind(alpha_diversity_index[, c(index)], sample_info[match(rownames(alpha_diversity_index), rownames(sample_info)),])
    colnames(data)[1] <- index
    p <- ggplot(data, aes(x = Group, y = data[, c(index)], color = Group)) +
        theme_bw() +
        theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
        theme(legend.text = element_text(size = 8)) +
        # geom_boxplot(alpha = 1, outlier.size = 0, size = 0.7, width = 0.5, fill = 'transparent') +
        geom_boxplot() +
        labs(x = 'Groups', y = index)
    # p
    ggsave(paste(path, '/03_Diversity/Alpha/alpha_div_collated/', index, '/', index, '.png', sep=""), width = 25, height = 18, units ="cm", dpi = 600)
    ggsave(paste(path, '/03_Diversity/Alpha/alpha_div_collated/', index, '/', index, '.pdf', sep=""))
}

##### 多样性指数统计检验 ########
cbind_data <- cbind(alpha_diversity_index, sample_info[match(rownames(alpha_diversity_index), rownames(sample_info)),])
group_num <- length(unique(sample_info$Group))

for (index in index_list) {
    dir.create(paste0(path, '/03_Diversity/Alpha/alpha_div_collated/', index))
    box_plot_func(index)

    data <- cbind_data[,c(index, 'Group')]
    colnames(data)[1] <- 'x'
    #### 2组ttest和wilcox, 3组及以上tukey和wilcox
    if (group_num == 2){
        groupvs <- paste0(group[1], '_vs_', group[2])       
        ## t-test
        ttest <- t.test(x ~ Group, data)  # 假设俩个变量方差不齐
        p_value <- ttest$p.value
        ttest_result <- data.frame(groupvs, p_value)
        write.table(ttest_result, file=paste0(path, '03_Diversity/Alpha/alpha_div_collated/', index, '/ttest.txt'), sep = '\t', row.names = F)
        ## wilcox
        ## wilcox.test()函数可以用来做Wilcoxon秩和检验，也可以用于做Mann-Whitney U检验。当参数为单个样本，或者是两个样本相减，或者是两个参数，paired=TRUE时，是Wilcoxon秩和检验。当paired = FALSE（独立样本）时，就是Mann-Whitney U检验
        wilcox <- wilcox.test(x ~ Group, data)
        p_value <- wilcox$p.value
        wilcox_result <- data.frame(groupvs, p_value)
        write.table(wilcox_result, file=paste0(path, '03_Diversity/Alpha/alpha_div_collated/', index, '/wilcox.txt'), sep = '\t', row.names = F)
    }else{
        # groupvs <- paste(unique(sample_info$Group), collapse = '_vs_')
        ## tukey
        ## 最小显著差数检验法(LSD法)
        ## Tukey氏固定差距检验法(Tukey HSD)
        ## 在单因素方差分析ANOVA中，如果该因素影响比较显著，那么需要进一步利用多重比较方法比较该因素不同水平的影响，确定不同水平下该因素的影响是否显著。常见的多重比较方法主要有两种，LSD法和Tukey HSD法, 无论是LSD还是Tukey进行多重比较，都需要先进行ANOVA分析
        fit <- aov(x ~ Group, data)
        tukey <- TukeyHSD(fit)
        tukey_result <- tukey$Group
        write.table(tukey_result, file=paste0(path, '03_Diversity/Alpha/alpha_div_collated/', index, '/tukey.txt'), sep = '\t')
        ##Kurskal-Wallis检验是Wilcoxon方法（其实是Mann-Whitney检验）用于多个样本。当对两个样本进行比较的时候，Kurskal-Wallis检验与Mann-Whitney检验是等价的。
        # wilcox <- kruskal.test(x ~ Group, data)
        wilcox <- kruskal(data$x, data$Group, group = F)
        wilcox_result <- wilcox$comparison
        write.table(wilcox_result, file=paste0(path, '03_Diversity/Alpha/alpha_div_collated/', index, '/wilcox.txt'), sep = '\t')

        # ## 每俩组wilcox
        # group2vs <- t(combn(unique(sample_info$Group), 2))
        # result_wilcox <- data.frame(groupvs = character(), p_value = character(), stringsAsFactors = F)
        # for (i in 1:nrow(group2vs)){
        #     groupvs <- paste0(group2vs[i, 1], '_vs_', group2vs[i, 2])
        #     tmp_data <- data[which((data$Group == group2vs[i, 1]) | (data$Group == group2vs[i, 2])),]
        #     wilcox <- wilcox.test(x ~ Group, tmp_data)
        #     p_value <- wilcox$p.value
        #     wilcox_result <- data.frame(groupvs, p_value)
        #     result_wilcox <- rbind(result_wilcox, wilcox_result, stringsAsFactors = F)
        # }
        # write.table(result_wilcox, file=paste0(path, '03_Diversity/Alpha/alpha_div_collated/', index, '/wilcox.txt'), sep = '\t', row.names = F)
    }
}


rm(list = ls())
# # 合并alpha指数和实验设计
# cbind_data <- cbind(alpha_diversity_index, sample_info[match(rownames(alpha_diversity_index), rownames(sample_info)),])
# group <- unique(sample_info$Group)
# group <- paste0(group[1], '_vs_', group[2])
# ### t test 检验 ###
# t_result <- t.test(simpson ~ Group, data = cbind_data)   # 假设俩个变量方差不齐(该结果准确)
# simpson_df <- data.frame(group, t_result$p.value)
# write.table(simpson_df, file="C:/Users/jbwang/Desktop/2/simpson.txt", sep = '\t', row.names = F)
# ### wilcox 非参检验 ####
# wilcox_result <- wilcox.test(simpson ~ Group, data = cbind_data)
# wilcox_df <- data.frame(group, wilcox_result$p.value)
# write.table(wilcox_df, file="C:/Users/jbwang/Desktop/2/wilcox.txt", sep = '\t', row.names = F)




# ######以observed_otus为例
# # 统计组间是否差异显著
# observed_otus_stats <- aov(observed_otus ~ Type, data = index)
# # 使用TukeyHSD对组间进行检验，效正pvalue
# TukeyHSD_observed_otus <- TukeyHSD(observed_otus_stats, ordered = FALSE, conf.level = 0.95)
# TukeyHSD_observed_otus_table <- as.data.frame(TukeyHSD_observed_otus$Type)
# # 保存结果到文件，按Pvaule值由小到大排序
# write.table(TukeyHSD_observed_otus_table[order(TukeyHSD_observed_otus_table$p, decreasing=FALSE), ], file="alpha_observed_otus_stats.txt",append = FALSE, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)


# ### T-test检验
# t1 <- t.test(observed_species ~ Group, data = cbind_data)  # 假设俩个变量方差不齐
# t1
# # t2 <- t.test(observed_otus ~ Type, data = index, var.equal = T)  # 假设方差齐性

# # ### Wilcoxon非参数秩检验
# w1 <- wilcox.test(observed_species ~ Group, data = cbind_data)
# w1

# ### 正态性检验
# opar<-par(no.readonly = TRUE)
# par(mfrow=c(1,2))

# dat1<-subset(index, Type=='N')
# dat2 <- subset(index, Type == 'BC')
# hist(dat1$observed_otus)
# hist(dat2$observed_otus)
# qqnorm(dat1$observed_otus)
# qqnorm(dat2$observed_otus)
# par(opar)
# shapiro.test(dat1$observed_otus)   # P>0.05 :正态分布，P<0.05 :非正态分布
# shapiro.test(dat2$observed_otus)

# boxplot(observed_otus ~ Type, index)

# ## 方差齐性检验
# library(car)
# leveneTest(index$observed_otus, index$Type, mean)  # P>0.05 :方差齐，P<0.05 :方差不齐


