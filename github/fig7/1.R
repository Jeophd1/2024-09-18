rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)

workDir = "C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig6"
setwd(workDir)
expFile = "diffGeneExp.txt"

data = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data = as.matrix(data)

group = sapply(strsplit(colnames(data), "\\_"), "[", 2)
data = data[, group == "treat"]

maxK = 7

results = ConsensusClusterPlus(
  data,
  maxK = maxK,
  reps = 50,
  pItem = 0.8,
  pFeature = 1,
  title = workDir,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 123456,
  plot = "png"
)

clusterNum = 2
cluster = results[[clusterNum]][["consensusClass"]]
cluster = as.data.frame(cluster)
colnames(cluster) = c("ACEcluster")
letter = c("A", "B", "C", "D", "E", "F", "G")
uniqClu = levels(factor(cluster$ACEcluster))
cluster$ACEcluster = letter[match(cluster$ACEcluster, uniqClu)]
outTab = cbind(t(data), cluster)
outTab = rbind(ID = colnames(outTab), outTab)
write.table(outTab, file = "ACEcluster.txt", sep = "\t", quote = FALSE, col.names = FALSE)
