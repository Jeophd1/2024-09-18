rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages("ggpubr")
library(limma)
library(ggpubr)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3G")
expFile = "GSE169007.txt"
conFile = "GSE169007_s1.txt"
treatFile = "GSE169007_s2.txt"
geneFile = "interGenes.txt"

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt = avereps(data)

qx = as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC = ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))
if (LogC) {
  rt[rt < 0] = 0
  rt = log2(rt + 1)
}
data = normalizeBetweenArrays(rt)

con = read.table(conFile, header = FALSE, sep = "\t", check.names = FALSE)
treat = read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
conData = data[, as.vector(con[, 1])]
treatData = data[, as.vector(treat[, 1])]
data = cbind(conData, treatData)
conNum = ncol(conData)
treatNum = ncol(treatData)

Type = c(rep("con", conNum), rep("treat", treatNum))
outData = rbind(id = paste0(colnames(data), "_", Type), data)
write.table(outData, file = "test.normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)

geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data = data[as.vector(geneRT[, 1]), , drop = FALSE]

Type = c(rep("Con", conNum), rep("Treat", treatNum))
my_comparisons = list()
my_comparisons[[1]] = levels(factor(Type))

newGeneLists = c()
outTab = data.frame()
for (i in row.names(data)) {
  rt1 = data.frame(expression = data[i, ], Type = Type)
  boxplot = ggboxplot(rt1, x = "Type", y = "expression", color = "Type",
                      xlab = "",
                      ylab = paste(i, "expression"),
                      legend.title = "",
                      palette = c("blue", "red"),
                      add = "jitter") +
    stat_compare_means(comparisons = my_comparisons)
  pdf(file = paste0("boxplot.", i, ".pdf"), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}
