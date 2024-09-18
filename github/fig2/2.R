install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("GSEABase")

rm(list = ls())

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig2\\原始代码")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GSEABase")
library("DOSE")

pvalueFilter = 0.05
qvalueFilter = 1

colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

rt = read.table("diff.txt", header = TRUE, sep = "\t", check.names = FALSE)

genes = as.vector(rt[, 1])
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs = as.character(entrezIDs)
rt = cbind(rt, entrezID = entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]

kk = enrichDO(gene = gene, ont = "DO", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
DO = as.data.frame(kk)
DO = DO[(DO$pvalue < pvalueFilter & DO$qvalue < qvalueFilter), ]

write.table(DO, file = "DO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

showNum = 30
if (nrow(DO) < showNum) {
  showNum = nrow(DO)
}

pdf(file = "barplot.pdf", width = 6, height = 10)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file = "bubble.pdf", width = 6, height = 10)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", color = colorSel)
dev.off()
