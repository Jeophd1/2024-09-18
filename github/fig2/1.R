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
qvalueFilter = 0.05

colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

rt = read.table("diff.txt", header = TRUE, sep = "\t", check.names = FALSE)

genes = as.vector(rt[, 1])
entrezIDs = mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs = as.character(entrezIDs)
gene = entrezIDs[entrezIDs != "NA"]

kk = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO = as.data.frame(kk)
GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ]

write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

showNum = 6
if (nrow(GO) < 30) {
  showNum = nrow(GO)
}

pdf(file = "barplot.pdf", width = 6, height = 9)
bar = barplot(kk, drop = TRUE, showCategory = showNum, split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

pdf(file = "bubble.pdf", width = 6, height = 9)
bub = dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()
