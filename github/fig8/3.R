rm(list = ls())
install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05  
qvalueFilter=0.2  
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig6\\S1")            
rt=read.table("interGene.txt", header=F, sep="\t", check.names=F)   
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=qvalueFilter)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="barplot.pdf", width=9, height=11)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()

pdf(file="bubble.pdf", width = 9, height = 11)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, color=colorSel)
dev.off()
