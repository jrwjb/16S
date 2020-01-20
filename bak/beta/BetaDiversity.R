suppressMessages(library(phyloseq))
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(vegan))
library(ropls)
library(stringr)
suppressMessages(library(car))
library(reshape2)
suppressMessages(library(agricolae))
library(RColorBrewer)

args <- commandArgs(T)
path <- as.character(args[1])

sample_info <- read.table(paste0(path, '/sample_info.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- read.table(paste0(path, '/03_Diversity/otu_table_even.txt'), row.names = 1, header = T, sep = '\t', check.names = F)
otu_table <- otu_table[, rownames(sample_info)]
tax <- read.table(paste0(path, '/02_OTU/tax.txt'), row.names = 1, header = F, sep = '\t')
colnames(tax) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
groupVs <- read.table(paste0(path, '/groupvs.txt'),header = FALSE,sep = '\t',fill = TRUE,quote = "",check.names = TRUE)

######################## PCA ############################
dir.create(paste0(path, '/03_Diversity/Beta/PCA'))
colpalette <- c("#1d953f","#102b6a","#c77eb5", "#fcf16e", "#2585a6", "purple", "#e0861a", "#d71345", "#6b473c", "#78a355", "#fdb933", "#5e7c85", "#411445", "#c37e00", "#bed742","#009ad6","#9d9087")
# data <- otu_table
# datasum <- apply(data[1:ncol(data)],2,sum,na.rm=T)
pca_func <- function(data, sample_info, groupvs){
    dir.create(paste0(path, '/03_Diversity/Beta/PCA/', groupvs))
    xMN <- t(data)
    # xMN <- xMN/datasum
    xMN2 <- xMN[,apply(xMN,2,function(x){is.na(x) %>% sum()})<(nrow(xMN)/2)]

    if (nrow(xMN2) < 7 ){
        compute_pca <- opls(x = xMN2,predI = NA,orthoI = 0,testL = FALSE,scaleC= "standard",plotL =FALSE, crossvalI = nrow(xMN2))
    }else{
        compute_pca <- opls(x = xMN2,predI = NA,orthoI = 0,testL = FALSE,scaleC= "standard",plotL =FALSE)
    }
    if (getSummaryDF(compute_pca)$pre == 1){
      if (nrow(xMN2) < 7 ){
          compute_pca <- opls(x = xMN2,predI = 2,orthoI = 0,testL = FALSE,scaleC= "standard",plotL =FALSE, crossvalI = nrow(xMN2))
      }else{
          compute_pca <- opls(x = xMN2,predI = 2,orthoI = 0,testL = FALSE,scaleC= "standard",plotL =FALSE)
      }
      # compute_pca  <- opls(xMN2,  predI = 2,orthoI = 0,testL = FALSE,scaleC= "standard",plotL = FALSE)
    }
    
    plot(compute_pca, typeVc ="x-score", parAsColFcVn = factor(sample_info$Group), parCexN = 0.8, parCompVi = c(1, 2), 
         parDevNewL = TRUE,parEllipsesL = FALSE, parLabVc =NA, parTitleL = TRUE,
         file.pdfC = paste0(path, '/03_Diversity/Beta/PCA/', groupvs, '/pca.pdf'),.sinkC = NULL) 
    modC_pca <-compute_pca@typeC  
    sumDF_pca <- getSummaryDF(compute_pca)
    desMC_pca <- compute_pca@descriptionMC
    scoreMN_pca <- getScoreMN(compute_pca)
    modelDF_pca <- compute_pca@modelDF
    loadingMN_pca <- getLoadingMN(compute_pca) 

    out_pca <- data.frame(Type = modC_pca,A = sumDF_pca$pre, N = desMC_pca[1],"R2X(cum)" = sumDF_pca$`R2X(cum)`, "R2Y(cum)"="-", "Q2(cum)"="-", Title = "dele-QC" ,stringsAsFactors = F) 
    OUT1 <- rbind(out_pca,stringsAsFactors = F) 
    colnames(OUT1) <- c("Type","A","N","R2X(cum)","R2Y(cum)",   "Q2(cum)","Title")
    write.table(OUT1, paste0(path, '/03_Diversity/Beta/PCA/', groupvs, '/pca.csv'), sep = ',',row.names = F, quote = F)

    samDF <- data.frame(name=rownames(xMN),group=sample_info$Group)
    ploMN0 <- scoreMN_pca
    hotFisN0 <- (nrow(ploMN0 ) - 1) * 2 * (nrow(ploMN0 )^2 - 1) / (nrow(ploMN0) * nrow(ploMN0) * (nrow(ploMN0) - 2)) * qf(0.95, 2, nrow(ploMN0) - 2)
    colorss0 <- as.numeric(samDF[,2])
    uniq.cols0 <- unique(colpalette[colorss0])
    if (length(unique(as.numeric(samDF[,2]))) >6){
      uniq.pchs0 <- unique(as.numeric(samDF[,2]))
    }else{
      uniq.pchs0 <- unique(as.numeric(samDF[,2]) + 14)
    }

    par(mar=c(12, 6, 4, 3)+0.1,family="RMN",font.lab=2,font.axis=1,cex.lab=1.5,cex.axis=1,las=2,xpd=TRUE) 
    png(filename = paste0(path, '/03_Diversity/Beta/PCA/', groupvs, '/pca.png'),width = 26,height = 22,units = "cm",res = 300)
    if (length(unique(as.numeric(samDF[,2]))) >6){
      plot(ploMN0,main=paste0("Scores (PCA)"),
           xlab = paste0("t", 1, " (",round(modelDF_pca[1, "R2X"] * 100),"%)"),
           ylab = paste0("t", 2, " (",round(modelDF_pca[2, "R2X"] * 100),"%)"),
           xlim = c(-1, 1) * max(sqrt(var(ploMN0[, 1]) * hotFisN0), max(abs(ploMN0[, 1]))),
           ylim =c(-1, 1) *max(sqrt(var(ploMN0[, 2]) * hotFisN0), max(abs(ploMN0[, 2]))),
           col = colpalette[colorss0],
           pch=as.numeric(samDF[,2]),cex=1.2)  
    }else{
      plot(ploMN0,main=paste0("Scores (PCA)"),
           xlab = paste0("t", 1, " (",round(modelDF_pca[1, "R2X"] * 100),"%)"),
           ylab = paste0("t", 2, " (",round(modelDF_pca[2, "R2X"] * 100),"%)"),
           xlim = c(-1, 1) * max(sqrt(var(ploMN0[, 1]) * hotFisN0), max(abs(ploMN0[, 1]))),
           ylim =c(-1, 1) *max(sqrt(var(ploMN0[, 2]) * hotFisN0), max(abs(ploMN0[, 2]))),
           col = colpalette[colorss0],
           pch=as.numeric(samDF[,2]) + 14,cex=1.2)   
    }

    abline(v = 0)
    abline(h = 0)
    radVn0 <- seq(0, 2 * pi, length.out = 100)
    lines(sqrt(var(ploMN0[, 1]) * hotFisN0) * cos(radVn0),sqrt(var(ploMN0[, 2]) * hotFisN0) * sin(radVn0))  ## Tenenhaus98, p87
    legend.nm0 <- unique(as.character(samDF[,2])) 
    if ( length(uniq.cols0) > 1 ) {
      names(uniq.cols0) <- legend.nm0;
    }
    legend("topright", legend = legend.nm0,  col=uniq.cols0,pch=uniq.pchs0)
    dev.off()
}


######################### NMDS ###########################
dir.create(paste0(path, '/03_Diversity/Beta/NMDS'))

nmds_func <- function(data, tax, sample_info, groupvs){
    dir.create(paste0(path, '/03_Diversity/Beta/NMDS/',groupvs))
    count_phy <- otu_table(data, taxa_are_rows = T)
    tax_phy <- tax_table(as.matrix(tax))
    sample_phy <- sample_data(sample_info)
    otu_phy <- phyloseq(count_phy, tax_phy, sample_phy)
    nmds_of_bray_curtis <- ordinate(physeq = otu_phy, distance = 'bray', method = 'NMDS')
    p <- plot_ordination(otu_phy, nmds_of_bray_curtis, type = 'samples', color = 'Group')
    p <- p + geom_point(size = 3) + ggtitle('NMDS') + stat_ellipse() + theme(text = element_text(size = 15)) + theme_bw() +
    	       theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
            geom_hline(aes(yintercept=0),linetype="dashed") + geom_vline(aes(xintercept=0),linetype="dashed")

    ggsave(paste0(path, '/03_Diversity/Beta/NMDS/',groupvs,'/NMDS.png'), width = 25, height = 18, units ="cm", dpi = 600)
    ggsave(paste0(path, '/03_Diversity/Beta/NMDS/',groupvs,'/NMDS.pdf'))
}


######################### PcoA ###########################
# dir.create(paste0(path, '/03_Diversity/Beta/PCoA'))
# 读入距离矩阵
weighted_unifrac_file <- paste0(path, '/03_Diversity/Beta/weighted_unifrac_dm.txt')
unweighted_unifrac_file <- paste0(path, '/03_Diversity/Beta/unweighted_unifrac_dm.txt')
weighted_unifrac <- read.table(weighted_unifrac_file, sep = '\t', row.names = 1, header = T, check.names = F)
unweighted_unifrac <- read.table(unweighted_unifrac_file, sep = '\t', row.names = 1, header = T, check.names = F)

# 对距离矩阵进行主坐标轴分析
pcoa_func <- function(data, file, groupvs){
    dir.create(paste0(path, '/03_Diversity/Beta/PCoA/',groupvs))
    # out <- sapply(strsplit(file, '\\.'), '[', 1)
    tmp <- unlist(strsplit(file, '/'))
    out <- paste(tmp[-length(tmp)], collapse = '/')
    index <- sapply(strsplit(tmp[length(tmp)], '\\.'), '[', 1)
    # index <- sapply(strsplit(index, '/'), '[', length(unlist(strsplit(index, '/'))))
    pcoa <- cmdscale(data, k = 3, eig = T) # k is dimension, 3 is recommended, eig is eigenvalue
    points <- as.data.frame(pcoa$points)
    colnames(points) <- c('x', 'y', 'z')
    eig <- pcoa$eig
    points <- cbind(points, sample_info[match(rownames(points), rownames(sample_info)),])

    p <- ggplot(points, aes(x = x, y = y, color = Group)) +
        geom_point(alpha = .7, size = 2) +
        stat_ellipse(type = 't', linetype = 2) +      ## 添加聚类圆
        theme_bw() +
        theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
        geom_hline(aes(yintercept=0),linetype="dashed") + geom_vline(aes(xintercept=0),linetype="dashed") +
        labs(x = paste('PCoA 1 (', format(100 * eig[1] / sum(eig), digits = 4), '%)', sep = ''),
             y = paste('PCoA 2 (', format(100 * eig[2] / sum(eig), digits = 4), '%)', sep = ''),
             title = '')

    ggsave(paste0(out, '/PCoA/', groupvs, '/', index, 'PCoA.png'), width = 25, height = 18, units ="cm", dpi = 600)
    ggsave(paste0(out, '/PCoA/', groupvs, '/', index, 'PCoA.pdf'))
}


dir.create(paste0(path, '/03_Diversity/Beta/adonis'))
dir.create(paste0(path, '/03_Diversity/Beta/anosim'))

# pca pcoa nmds所有组
num <- length(unique(sample_info$Group))
if (num > 2){
  pca_func(otu_table, sample_info, 'all')
  nmds_func(otu_table, tax, sample_info, 'all')
  pcoa_func(weighted_unifrac, weighted_unifrac_file, 'all')
  pcoa_func(unweighted_unifrac, unweighted_unifrac_file, 'all')
}
## 比较组
result_anosim <- data.frame(groupvs = character(), R = character(), p_value = character(), stringsAsFactors = F)
result_adonis <- data.frame(groupvs = character(), R2 = character(), p_value = character(), stringsAsFactors = F)
for (i in 1:nrow(groupVs)) {
    groupvs <- groupVs[i,1]
    group_list <- unlist(strsplit(as.character(groupVs[i,1]),"_vs_"))
    tmp_sample_info <- sample_info[sample_info$Group %in% group_list,]
    tmp_sample_info <- droplevels(tmp_sample_info)
    tmp_otu_table <- otu_table[,rownames(tmp_sample_info)]

    pca_func(tmp_otu_table, tmp_sample_info, groupvs)
    nmds_func(tmp_otu_table, tax, tmp_sample_info, groupvs)
    tmp_weighted_unifrac <- weighted_unifrac[rownames(tmp_sample_info), rownames(tmp_sample_info)]
    pcoa_func(tmp_weighted_unifrac, weighted_unifrac_file, groupvs)
    tmp_unweighted_unifrac <- unweighted_unifrac[rownames(tmp_sample_info), rownames(tmp_sample_info)]
    pcoa_func(tmp_unweighted_unifrac, unweighted_unifrac_file, groupvs)

    ## anosim
    dm <- vegdist(t(tmp_otu_table), method = 'bray')
    anosim <- anosim(dm, tmp_sample_info$Group, permutations = 999)
    R <- anosim$statistic
    p_value <- anosim$signif
    anosim_result <- data.frame(groupvs, R, p_value)
    result_anosim <- rbind(result_anosim, anosim_result, stringsAsFactors = F)
    ### anosim plot
    cols <- brewer.pal(n=3, name="Set1")
    png(filename = paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.png'), width = 26, height = 22, units = "cm", res = 300)
    plot(anosim, col=cols)
    dev.off()
    #pdf
    pdf(paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.pdf'))
    plot(anosim, col=cols)
    dev.off()
    ### adonis
    adonis <- adonis(dm~tmp_sample_info$Group, permutations = 999)
    adonis <- adonis$aov.tab
    R2 <- adonis$R2[1]
    p_value <- adonis$`Pr(>F)`[1]
    adonis_result <- data.frame(groupvs, R2, p_value)
    result_adonis <- rbind(result_adonis, adonis_result, stringsAsFactors = F)

}
write.table(result_anosim, paste0(path, '/03_Diversity/Beta/anosim/anosim_result.csv'), sep = ',', row.names = F, quote = F)
write.table(result_adonis, paste0(path, '/03_Diversity/Beta/adonis/adonis_result.csv'), sep = ',', row.names = F, quote = F)


########### distance boxplot ###########
distance_unweight <- read.table(paste0(path, '/03_Diversity/Beta/Beta_div/unweighted_unifrac/Group_Distances.txt'), sep = '\t', fill = T)
distance_weight <- read.table(paste0(path, '/03_Diversity/Beta/Beta_div/weighted_unifrac/Group_Distances.txt'), sep = '\t', fill = T)

boxplot_func <- function(data, index) {
    group <- c()
    for (i in 1:length(data[, 1])){
        group[i] <- sapply(strsplit(as.character(data[i, 1]), '_vs._'), '[', 1)
    }
    data[, 1] <- group
    colnames(data)[1] <- 'Group'
    # data2 <- melt(data, id.vars = c('Group'))

    result <- data.frame(groupvs = character(), p_value = character(), method = character(), stringsAsFactors = F)
    for (i in 1:nrow(groupVs)) {
      groupvs <- groupVs[i,1]
      group_list <- unlist(strsplit(as.character(groupVs[i,1]),"_vs_"))
      group_num <- length(group_list)

      if (group_num == 2){
        tmp_data <- data[data$Group %in% group_list,]
        data2 <- melt(tmp_data, id.vars = c('Group'))
        normality_test <- shapiro.test(data2$value)   # P>0.05 :正态分布，P<0.05 :非正态分布
        if (normality_test$p.value > 0.05){
          ttest <- t.test(value ~ Group, data2)  # 假设俩个变量方差不齐
          p_value <- ttest$p.value
          ttest_result <- data.frame(groupvs, p_value)
          ttest_result$method <- 't-test'
          result <- rbind(result, ttest_result, stringsAsFactors = F)
          # write.table(ttest_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/ttest.csv'), sep = ',', row.names = F)
        }else{
          wilcox <- wilcox.test(value ~ Group, data2)
          p_value <- wilcox$p.value
          wilcox_result <- data.frame(groupvs, p_value)
          wilcox_result$method <- 'Wilcoxon'
          result <- rbind(result, wilcox_result, stringsAsFactors = F)
          # write.table(wilcox_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/wilcox.csv'), sep = ',', row.names = F)
        }       
      }else{
        tmp_data <- data[data$Group %in% group_list,]
        data2 <- melt(tmp_data, id.vars = c('Group'))
        normality_test <- shapiro.test(data2$value)
        if (normality_test$p.value > 0.05){
          fit <- aov(value ~ Group, data2)
          tukey <- TukeyHSD(fit)
          tukey_result <- tukey$Group
          tukey_result <- as.data.frame(tukey_result)
          tukey_result$method <- 'tukey'
          # result <- rbind(result, tukey_result, stringsAsFactors = F)
          write.table(tukey_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/', groupvs, '.csv'), sep = ',', quote = F)
        }else{
          ##Kurskal-Wallis检验是Wilcoxon方法（其实是Mann-Whitney检验）用于多个样本。当对两个样本进行比较的时候，Kurskal-Wallis检验与Mann-Whitney检验是等价的。
          # wilcox <- kruskal.test(x ~ Group, data)
          wilcox <- kruskal(data2$value, data2$Group, group = F)
          wilcox_result <- wilcox$comparison
          wilcox_result$method <- 'kruskal'
          # result <- rbind(result, tukey_result, stringsAsFactors = F)
          write.table(wilcox_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index,  '/', groupvs, '.csv'), sep = ',', quote = F)
        }
    }
  }
  write.table(result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/result.csv'), sep = ',', row.names = F, quote = F)

  data_box <- melt(data, id.vars = c('Group'))
  p <- ggplot(as.data.frame(na.omit(data_box)), aes(Group, value, color = Group)) +
        geom_boxplot() +
        theme_bw() +
        theme(panel.grid=element_blank(),axis.line=element_line(size=0.5,colour="black")) +
        theme(legend.text = element_text(size = 8)) +
        # geom_boxplot(alpha = 1, outlier.size = 0, size = 0.7, width = 0.5, fill = 'transparent') +
        labs(y = index)
  ggsave(paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/', index, '.png'), width = 25, height = 18, units ="cm", dpi = 600)
  ggsave(paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/', index, '.pdf'))
}

boxplot_func(distance_weight, 'weighted_unifrac')
boxplot_func(distance_unweight, 'unweighted_unifrac')


rm(list = ls())

    # # 统计检验
    # if (nrow(data) == 2){
    #     groupvs <- paste(data$Group, collapse = '_vs_')

    #     ## 判断是否是正太分布
    #     normality_test <- shapiro.test(data2$value)   # P>0.05 :正态分布，P<0.05 :非正态分布
    #     if (normality_test$p.value > 0.05){
    #       ttest <- t.test(value ~ Group, data2)  # 假设俩个变量方差不齐
    #       p_value <- ttest$p.value
    #       ttest_result <- data.frame(groupvs, p_value)
    #       write.table(ttest_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/ttest.csv'), sep = ',', row.names = F)
    #     }
    #     else{
    #       wilcox <- wilcox.test(value ~ Group, data2)
    #       p_value <- wilcox$p.value
    #       wilcox_result <- data.frame(groupvs, p_value)
    #       write.table(wilcox_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/wilcox.csv'), sep = ',', row.names = F)
    #     }       
    # }
    # else{
    #     normality_test <- shapiro.test(data2$value)
    #     if (normality_test$p.value > 0.05){
    #       fit <- aov(value ~ Group, data2)
    #       tukey <- TukeyHSD(fit)
    #       tukey_result <- tukey$Group
    #       write.table(tukey_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/tukey.csv'), sep = ',')
    #     }
    #     else{
    #       ##Kurskal-Wallis检验是Wilcoxon方法（其实是Mann-Whitney检验）用于多个样本。当对两个样本进行比较的时候，Kurskal-Wallis检验与Mann-Whitney检验是等价的。
    #       # wilcox <- kruskal.test(x ~ Group, data)
    #       wilcox <- kruskal(data2$value, data2$Group, group = F)
    #       wilcox_result <- wilcox$comparison
    #       write.table(wilcox_result, file=paste0(path, '/03_Diversity/Beta/Beta_div/', index, '/wilcox.csv'), sep = ',')
    #     }
    # }

# ######################### Anosim 、Adonis ###########################
# dir.create(paste0(path, '/03_Diversity/Beta/adonis'))
# dir.create(paste0(path, '/03_Diversity/Beta/anosim'))
# ####每两组比较
# group2vs <- t(combn(unique(sample_info$Group), 2))
# result_anosim <- data.frame(groupvs = character(), R = character(), p_value = character(), stringsAsFactors = F)
# result_adonis <- data.frame(groupvs = character(), R2 = character(), p_value = character(), stringsAsFactors = F)
# for (i in 1:nrow(group2vs)){
#     group1 <- grep(paste0('^', group2vs[i, 1], '$'), sample_info$Group, value = T)
#     group2 <- grep(paste0('^', group2vs[i, 2], '$'), sample_info$Group, value = T)
#     tmp_groupvs <- c(group1, group2)
#     groupvs <- paste(unique(tmp_groupvs), collapse = '_vs_')
#     tmp_sample_info <- sample_info[which((sample_info$Group == group2vs[i, 1]) | (sample_info$Group == group2vs[i, 2])),]
#     tmp_otu_table <- otu_table[, rownames(tmp_sample_info)]
#     dm <- vegdist(t(tmp_otu_table), method = 'bray')
#     ### anosim
#     anosim <- anosim(dm, tmp_groupvs, permutations = 999)   
#     R <- anosim$statistic
#     p_value <- anosim$signif
#     anosim_result <- data.frame(groupvs, R, p_value)
#     result_anosim <- rbind(result_anosim, anosim_result, stringsAsFactors = F)
#     ### anosim plot
#     cols <- brewer.pal(n=3, name="Set1")
#     png(filename = paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.png'), width = 26, height = 22, units = "cm", res = 300)
#     plot(anosim, col=cols)
#     dev.off()
#     #pdf
#     pdf(paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.pdf'))
#     plot(anosim, col=cols)
#     dev.off()
#     ### adonis
#     adonis <- adonis(dm~tmp_groupvs, permutations = 999)
#     adonis <- adonis$aov.tab
#     R2 <- adonis$R2[1]
#     p_value <- adonis$`Pr(>F)`[1]
#     adonis_result <- data.frame(groupvs, R2, p_value)
#     result_adonis <- rbind(result_adonis, adonis_result, stringsAsFactors = F)
# } 
# write.table(result_anosim, paste0(path, '/03_Diversity/Beta/anosim/anosim_result.csv'), sep = ',', row.names = F)
# write.table(result_adonis, paste0(path, '/03_Diversity/Beta/adonis/adonis_result.csv'), sep = ',', row.names = F)

# #####多组一起比较
# ## anosim
# dm <- vegdist(t(otu_table), method = 'bray')
# groupvs <- paste(unique(sample_info$Group), collapse = '_vs_')
# anosim <- anosim(dm, sample_info$Group, permutations = 999)
# group_num <- length(unique(sample_info$Group))
# cols <- brewer.pal(n=group_num + 1, name="Set1")
# R <- anosim$statistic
# p_value <- anosim$signif
# anosim_result <- data.frame(groupvs, R, p_value)
# write.table(anosim_result, paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.csv'), sep = ',', row.names = F)
# png(filename = paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.png'), width = 26, height = 22, units = "cm", res = 300)
# plot(anosim, col=cols)
# dev.off()
# #pdf
# pdf(paste0(path, '/03_Diversity/Beta/anosim/', groupvs, '.pdf'))
# plot(anosim, col=cols)
# dev.off()

# ## adonis
# adonis <- adonis(dm~sample_info$Group, permutations = 999)
# adonis <- adonis$aov.tab
# R2 <- adonis$R2[1]
# p_value <- adonis$`Pr(>F)`[1]
# adonis_result <- data.frame(groupvs, R2, p_value)
# write.table(adonis_result, paste0(path, '/03_Diversity/Beta/adonis/', groupvs, '.csv'), sep = ',', row.names = F)
