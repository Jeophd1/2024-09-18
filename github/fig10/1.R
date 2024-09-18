rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

library(limma)
library(ConsensusClusterPlus)

workDir = "C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig8\\fig8ABCD"      
setwd(workDir)   
expFile = "101interGeneExp.txt"        

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
data = data[rowMeans(data) > 0, ]

maxK = 7      
results = ConsensusClusterPlus(data,
                               maxK = maxK,
                               reps = 50,
                               pItem = 0.8,
                               pFeature = 1,
                               title = workDir,
                               clusterAlg = "km",
                               distance = "euclidean",
                               seed = 123456,
                               plot = "png")

clusterNum = 2        
cluster = results[[clusterNum]][["consensusClass"]]
cluster = as.data.frame(cluster)
colnames(cluster) = c("geneCluster")
letter = c("A","B","C","D","E","F","G")
uniqClu = levels(factor(cluster$geneCluster))
cluster$geneCluster = letter[match(cluster$geneCluster, uniqClu)]
clusterOut = rbind(ID = colnames(cluster), cluster)
write.table(clusterOut, file = "geneCluster.txt", sep = "\t", quote = FALSE, col.names = FALSE)
