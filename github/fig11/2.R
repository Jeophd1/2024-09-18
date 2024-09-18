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

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig8\\fig8F")    
expFile = "diffGeneExp.txt"     
clusterFile = "geneCluster.txt"   

exp = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
exp = t(exp)

cluster = read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

sameSample = intersect(row.names(exp), row.names(cluster))
exp = exp[sameSample, , drop = FALSE]
cluster = cluster[sameSample, , drop = FALSE]
rt = cbind(exp, cluster)
rt = rt[order(rt$geneCluster), ]

data = melt(rt, id.vars = c("geneCluster"))
colnames(data) = c("geneCluster", "Gene", "Expression")

bioCol = c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol = bioCol[1:length(levels(factor(data[, "geneCluster"])))]
p = ggboxplot(data, x = "Gene", y = "Expression", color = "geneCluster", 
              xlab = "",
              ylab = "Gene expression",
              legend.title = "geneCluster",
              palette = bioCol,
              width = 1)
p = p + rotate_x_text(60)
p1 = p + stat_compare_means(aes(group = geneCluster),
                            symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                            label = "p.signif")

pdf(file = "boxplot.pdf", width = 12, height = 5)
print(p1)
dev.off()
