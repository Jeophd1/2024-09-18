install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

rm(list = ls())

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig2\\原始代码")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

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

kk = enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
KEGG = as.data.frame(kk)
KEGG$geneID = as.character(sapply(KEGG$geneID, function(x) paste(rt$id[match(strsplit(x, "/")[[1]], as.character(rt$entrezID))], collapse = "/")))
KEGG = KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter), ]

write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

showNum = 30
if (nrow(KEGG) < showNum) {
  showNum = nrow(KEGG)
}

pdf(file = "barplot.pdf", width = 8, height = 9)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

pdf(file = "bubble.pdf", width = 8, height = 9)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", color = colorSel)
dev.off()
