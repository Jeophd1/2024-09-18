rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages("ggplot2")

library(limma)
library(ggplot2)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig6\\fig6G\\回复审稿人意见")  
clusterFile="ACEcluster.txt"   

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
ACEcluster=as.vector(rt[,ncol(rt)])

data.pca=prcomp(data, scale.=TRUE)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], ACEcluster= ACEcluster)
PCA.mean=aggregate(PCA[,1:2], list(ACEcluster=PCA$ACEcluster), mean)

variance_explained = data.pca$sdev^2 / sum(data.pca$sdev^2)
PC1_var_explained = round(variance_explained[1] * 100, 2)
PC2_var_explained = round(variance_explained[2] * 100, 2)

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ACECluCol=bioCol[1:length(levels(factor(ACEcluster)))]

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(factor(PCA$ACEcluster))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$ACEcluster==g,],
                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),
                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,
                                                                   center=c(mean(PC1),mean(PC2))))), ACEcluster=g))
}

pdf(file="PCA.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = ACEcluster)) +
  scale_colour_manual(name="ACEcluster", values =ACECluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=ACEcluster), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$ACEcluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(paste0("PC1 (", PC1_var_explained, "% variance)")) +
  ylab(paste0("PC2 (", PC2_var_explained, "% variance)"))
dev.off()

print(paste("PC1 variance explained:", PC1_var_explained, "%"))
print(paste("PC2 variance explained:", PC2_var_explained, "%"))
