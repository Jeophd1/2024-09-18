rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("ggpubr")

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig6")  
clusterFile="ACEcluster.txt"   

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$ACEcluster),]
data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
ACECluCol=bioCol[1:length(levels(factor(Type$ACEcluster)))]
names(ACECluCol)=levels(factor(Type$ACEcluster))
ann_colors[["ACEcluster"]]=ACECluCol

pdf("heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

sum(is.na(data))
sum(is.infinite(data))

svg("heatmap.svg", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

pdf("C:/Users/JeoPhD/Desktop/PCOS2/heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()

data=melt(rt, id.vars=c("ACEcluster"))
colnames(data)=c("ACEcluster", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "ACEcluster", 
            ylab="Gene expression",
            xlab="",
            legend.title="ACEcluster",
            palette = ACECluCol,
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=ACEcluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=12, height=5)
print(p1)
dev.off()
