rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("GSEABase")
BiocManager::install("GSVA")
BiocManager::install("DelayedMatrixStats")

install.packages("ggpubr")

library(reshape2) 
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(DelayedMatrixStats)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig7\\fig7B")  
expFile="normalize.txt"        
gmtFile="immune.gmt"           
clusterFile="ACEcluster.txt"  

original_data <- read.table("修改用.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
split_data <- split(original_data, original_data$V1)
result_data <- data.frame()
for (key in names(split_data)) {
  row <- c(key, unlist(split_data[[key]][, -1]))
  result_data <- rbind(result_data, row)
}
write.table(result_data, "your_output_file.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)

geneSets = getGmt(gmtFile, geneIdType = SymbolIdentifier())
ssgseaScore = gsva(data, geneSets, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

normalize = function(x){
  return((x - min(x)) / (max(x) - min(x)))
}
ssgseaScore = normalize(ssgseaScore)

ssgseaOut = rbind(id = colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut, file = "ssGSEA.result.txt", sep = "\t", quote = FALSE, col.names = FALSE)

cluster = read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

ssgseaScore = t(ssgseaScore)
sameSample = intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore = ssgseaScore[sameSample, , drop = FALSE]
cluster = cluster[sameSample, "ACEcluster", drop = FALSE]
scoreCluster = cbind(ssgseaScore, cluster)

data = melt(scoreCluster, id.vars = c("ACEcluster"))
colnames(data) = c("ACEcluster", "Immune", "Fraction")

bioCol = c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol = bioCol[1:length(levels(factor(data[,"ACEcluster"])))]

p = ggboxplot(data, x = "Immune", y = "Fraction", color = "ACEcluster",
              xlab = "", 
              ylab = "Immune infiltration",
              legend.title = "ACEcluster",
              palette = bioCol)
p = p + rotate_x_text(50)

pdf(file = "boxplot.pdf", width = 8, height = 5)
p + stat_compare_means(aes(group = ACEcluster),
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
                       label = "p.signif")
dev.off()

summary(data$Fraction)
hist(data$Fraction)
data <- subset(data, Fraction < quantile(Fraction, 0.99))
ggsave(file = "boxplot.pdf", plot = p, width = 12, height = 8, dpi = 300)

library(Cairo)
CairoPDF(file = "boxplot.pdf", width = 20, height = 8)
print(p)
dev.off()

p <- ggboxplot(data, x = "Immune", y = "Fraction", color = "ACEcluster",
               xlab = "", ylab = "Immune infiltration",
               legend.title = "ACEcluster", palette = bioCol) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = "boxplot_adjusted.pdf", plot = p, width = 20, height = 20, dpi = 300)

p <- ggboxplot(data, x = "Immune", y = "Fraction", color = "ACEcluster",
               xlab = "", ylab = "Immune infiltration",
               legend.title = "ACEcluster", palette = bioCol) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        text = element_text(size = 12))
ggsave(file = "boxplot_optimized.pdf", plot = p, width = 20, height = 20, dpi = 300)
