library(VennDiagram)
library(plotrix)

args <- commandArgs(T)
path <- as.character(args[1])
# setwd(path)

sample_info <- read.table(paste0(path, '/sample_info.txt'), row.names = 1, header = T, sep = '\t', quote = "")
otu_table <- read.table(paste0(path, '/02_OTU/otu_table_tax.txt'), row.names = 1, header = T, sep = '\t', quote = "", check.names = F)
groupVs <- read.table(paste0(path, '/groupvs.txt'),header = FALSE,sep = '\t',fill = TRUE,quote = "",check.names = TRUE)

dir.create(paste0(path, '/02_OTU/Venn'))
out <- paste0(path, '/02_OTU/Venn')
thresh <- 0

ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')

flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
    par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
    plot(c(0,10),c(0,10),type='n')
    n   <- length(sample)
    deg <- 360 / n
    res <- lapply(1:n, function(t){
        draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                     y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                     col = ellipse_col[t],
                     border = ellipse_col[t],
                     a = a, b = b, angle = deg * (t - 1))
        text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
             otu_num[t])
        
        if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                 sample[t],
                 srt = deg * (t - 1) - start,
                 adj = 1,
                 cex = 1
            )
        } else {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                 sample[t],
                 srt = deg * (t - 1) + start,
                 adj = 0,
                 cex = 1
            )
        }            
    })
    draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
    text(x = 5, y = 5, paste('Core:', core_otu))
}


venn_func <- function(otu_table, sample_info, out){

    tax <- otu_table[,dim(otu_table)[2]]
    data <- otu_table[, -dim(otu_table)[2]]
    
    group <- sample_info['Group']
    
    sum <- colSums(t(data))
    
    func1 <- function(x) {
        return(tapply(x, INDEX = group, sum) > thresh)
    }
    
    tb1 <- apply(data, 1, func1)
    unig <- rownames(tb1)
    d1 <- colSums(tb1)
    ll <- length(unig)
    com <- sum(d1 == ll)
    com_spe <- data.frame(OTUID = colnames(tb1)[d1 == ll], taxonomy = tax[d1 == ll], Sum_Abundance = sum[d1 == ll])
    
    write.table(com_spe, paste(out,'/common.csv', sep = ''), row.names = F, sep = ',')
    
    flower_data <- c()
    group_data <- c()
    for (i in 1:ll){
        # is_uniq <- (d1 == 1 & tb1[i, 1] == 1)
        # flower_data[i] <- sum(is_uniq)
        flower_data[i] <- sum(tb1[i,])
        group_data[i] <- rownames(tb1)[i]
        uniq_spe <- data.frame(OTUID = colnames(tb1)[d1 == 1 & tb1[i,] == 1], taxonomy = tax[d1 == 1 & tb1[i,] == 1], Sum_Abundance = sum[d1 == 1 & tb1[i,] == 1])
        write.table(uniq_spe, paste0(out,'/', 'group_', unig[i],'.csv'), row.names = F, sep = ',')
    }
    
    if (length(unig) <= 5){
        ls1 <- list()
        for (i in 1:ll){
            ls1[[i]] <- colnames(tb1)[tb1[i,]]
        }
        names(ls1) <- rownames(tb1)
        venn.diagram(ls1, filename=paste0(out,'/venn.png'), alpha = 0.5, lwd = 1.2, cat.cex = 1.4, fill = rainbow(length(ls1)), margin = 0.15)
    }else{
        png(filename=paste0(out,'/venn.png'), width = 18,height = 18,units ="cm",bg="white",res=600)
        flower_plot(group_data, flower_data, nrow(com_spe), start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
        dev.off()
    }
}


num <- length(unique(sample_info$Group))
if (num == 2){
    out <- paste0(out, '/', groupVs[1,1])
    dir.create(out)
    venn_func(otu_table, sample_info, out)
}else{
    venn_func(otu_table, sample_info, out)
    ## 按组画图
    for (i in 1:nrow(groupVs)) {
        groupvs <- groupVs[i,1]
        dir.create(paste0(path, '/02_OTU/Venn/', groupvs))
        group_list <- unlist(strsplit(as.character(groupVs[i,1]),"_vs_"))
        tmp_sample_info <- sample_info[sample_info$Group %in% group_list,]
        tmp_sample_info <- droplevels(tmp_sample_info)
        tmp_otu_table <- otu_table[,c(rownames(tmp_sample_info), 'taxonomy')]
        tmp_out <- paste0(path, '/02_OTU/Venn/', groupvs)
        venn_func(tmp_otu_table, tmp_sample_info, tmp_out)
    }
    
}


rm(list = ls())
