rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

install.packages("ggpubr")
library(limma)
library(reshape2)
library(ggpubr)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig8\\fig8G")  
clusterFile = "geneCluster.txt"        
ssgseaFile = "ssGSEA.result.txt"       

cluster = read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
ssgsea = read.table(ssgseaFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
ssgsea = t(ssgsea)
sameSample = intersect(row.names(cluster), row.names(ssgsea))
cluster = cluster[sameSample, "geneCluster", drop = FALSE]
ssgsea = ssgsea[sameSample, , drop = FALSE]
scoreCluster = cbind(ssgsea, cluster)

data = melt(scoreCluster, id.vars = c("geneCluster"))
colnames(data) = c("geneCluster", "Immune", "Fraction")

bioCol = c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol = bioCol[1:length(levels(factor(data[, "geneCluster"])))]
p = ggboxplot(data, x = "Immune", y = "Fraction", color = "geneCluster", 
              xlab = "",
              ylab = "Immune infiltration",
              legend.title = "geneCluster",
              palette = bioCol)
p = p + rotate_x_text(50)
pdf(file = "boxplot.pdf", width = 8, height = 6)                        
p + stat_compare_means(aes(group = geneCluster),
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                       label = "p.signif")
dev.off()
