rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

install.packages("pheatmap")
library(limma)
library(pheatmap)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig8\\fig8E")     
expFile = "101interGeneExp.txt"     
clusterFile = "geneCluster.txt"    

exp = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
exp = t(exp)

cluster = read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

sameSample = intersect(row.names(exp), row.names(cluster))
exp = exp[sameSample, , drop = FALSE]
cluster = cluster[sameSample, , drop = FALSE]
rt = cbind(exp, cluster)
rt = rt[order(rt$geneCluster), ]

data = t(rt[, 1:(ncol(rt) - 1), drop = FALSE])
Type = rt[, ncol(rt), drop = FALSE]
bioCol = c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
ann_colors = list()
geneCluCol = bioCol[1:length(levels(factor(Type$geneCluster)))]
names(geneCluCol) = levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]] = geneCluCol
pdf(file = "heatmap.pdf", width = 8, height = 12)
pheatmap(data,
         annotation = Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 6)
dev.off()
