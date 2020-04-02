# library(ggbiplot)
options(warn =-1)
suppressMessages(library(ropls))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(car))

args <- commandArgs(T)
path <- as.character(args[1])
file <- as.character(args[2])
# setwd(path)
colpalette <- c("#1d953f","#102b6a","#c77eb5", "#fcf16e", "#2585a6", "purple", "#e0861a", "#d71345", "#6b473c", "#78a355", "#fdb933", "#5e7c85", "#411445", "#c37e00", "#bed742","#009ad6","#9d9087")
sample_info <- read.table('sample_info.txt', header = T, row.names = 1, sep = '\t')
groupVs <- read.table('groupvs.txt',header = FALSE,sep = '\t',fill = TRUE,quote = "",check.names = TRUE)
data <- read.table(file, row.names = 1, header = T, sep = '\t', check.names = F)
# data <- data[,-c(ncol(data))]
#datasum <- apply(data[1:ncol(data)],2,sum,na.rm=T)
pca_func <- function(data, sample_info, groupvs){
  dir.create(paste0(path, groupvs))
  xMN <- t(data)
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
  plot(compute_pca, typeVc ="x-score", parAsColFcVn = factor(sample_info$Group),parCexN = 0.8, parCompVi = c(1, 2), 
       parDevNewL = TRUE,parEllipsesL = FALSE, parLabVc =NA, parTitleL = TRUE,
       file.pdfC = paste0(path, groupvs, "/pca.pdf"),.sinkC = NULL) 
  modC_pca <-compute_pca@typeC  
  sumDF_pca <- getSummaryDF(compute_pca)
  desMC_pca <- compute_pca@descriptionMC
  scoreMN_pca <- getScoreMN(compute_pca)
  modelDF_pca <- compute_pca@modelDF
  loadingMN_pca <- getLoadingMN(compute_pca) 

  out_pca <- data.frame(Type = modC_pca,A = sumDF_pca$pre, N = desMC_pca[1],"R2X(cum)" = sumDF_pca$`R2X(cum)`, "R2Y(cum)"="-", "Q2(cum)"="-", Title = "dele-QC" ,stringsAsFactors = F) 
  OUT1 <- rbind(out_pca,stringsAsFactors = F) 
  colnames(OUT1) <- c("Type","A","N","R2X(cum)","R2Y(cum)",	"Q2(cum)","Title")
  write.table(OUT1,paste0(path, groupvs, "/pca.csv"), sep=',')
  #############################################
  # samDF <- data.frame(name=rownames(xMN),group=str_extract(rownames(xMN), '[A-Z]*'))
  samDF <- data.frame(name=rownames(xMN),group=sample_info$Group)
  ###########################################
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
  png(filename = paste0(path, groupvs, "/pca.png"),width = 26,height = 22,units = "cm",res = 300)
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
  legend.nm0 <- unique(as.character(samDF[,2])) #ȡ??sort
  if ( length(uniq.cols0) > 1 ) {
    names(uniq.cols0) <- legend.nm0;
  }
  legend("topright", legend = legend.nm0,  col=uniq.cols0,pch=uniq.pchs0)
  #dataEllipse(ploMN0[1:8,1],ploMN0[1:8,2],add = TRUE,col = colpalette[1],levels = 0.8)
  #dataEllipse(ploMN0[9:27,1],ploMN0[9:27,2],add = TRUE,col = colpalette[2],levels = 0.8)
  #dataEllipse(ploMN0[28:32,1],ploMN0[28:32,2],add = TRUE,col = colpalette[3],levels = 0.8)
  dev.off()
}

# pca所有组
num <- length(unique(sample_info$Group))
if (num > 2){
  pca_func(data, sample_info, 'all')
}

for (i in 1:nrow(groupVs)) {
    groupvs <- groupVs[i,1]
    group_list <- unlist(strsplit(as.character(groupVs[i,1]),"_vs_"))
    tmp_sample_info <- sample_info[sample_info$Group %in% group_list,]
    tmp_sample_info <- droplevels(tmp_sample_info)
    tmp_data <- data[,rownames(tmp_sample_info)]

    pca_func(tmp_data, tmp_sample_info, groupvs)
}

rm(list = ls())