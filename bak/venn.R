library(VennDiagram)

args <- commandArgs(T)
path <- as.character(args[1])
# setwd(path)

sample_info <- read.table(paste0(path, '/sample_info.txt'), row.names = 1, header = T, sep = '\t')
otu_table <- read.table(paste0(path, '/03_Diversity/otu_table_even_tax.txt'), row.names = 1, header = T, sep = '\t')
dir.create(paste0(path, '/02_OTU/Venn'))
out <- paste0(path, '/02_OTU/Venn')
thresh <- 0

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

write.table(com_spe, paste(out,'/Group_common.csv', sep = ''), row.names = F, sep = ',')

flower_data <- c()
for (i in 1:ll){
    is_uniq <- (d1 == 1 & tb1[i, 1] == 1)
    flower_data[i] <- sum(is_uniq)
    uniq_spe <- data.frame(OTUID = colnames(tb1)[is_uniq], taxonomy = tax[is_uniq], Sum_Abundance = sum[is_uniq])
    write.table(uniq_spe, paste0(out,'/Group_', unig[i],'.csv'), row.names = F, sep = ',')
}

if (length(unig) <= 5){
    ls1 <- list()
    for (i in 1:ll){
        ls1[[i]] <- colnames(tb1)[tb1[i,]]
    }
    names(ls1) <- rownames(tb1)
    venn.diagram(ls1, filename = paste0(out,'/Group_venn.png'), alpha = 0.5, lwd = 1.2, cat.cex = 1.4, fill = rainbow(length(ls1)), margin = 0.15)
}

rm(list = ls())
