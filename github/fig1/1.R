rm(list = ls())
graphics.off()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
library(limma)
library(pheatmap)
setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig1A")
inputFile = "ACEGeneExp.txt"
logFCfilter = 0.8
adj.P.Val.Filter = 0.05
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
data = data[rowMeans(data) > 0, ]
sampleName1 = c()
files = dir()
files = grep("s1.txt$", files, value = TRUE)
for (file in files) {
  rt = read.table(file, header = FALSE, sep = "\t", check.names = FALSE)
  geneNames = as.vector(rt[, 1])
  uniqGene = unique(geneNames)
  sampleName1 = c(sampleName1, uniqGene)
}
sampleName2 = c()
files = dir()
files = grep("s2.txt$", files, value = TRUE)
for (file in files) {
  rt = read.table(file, header = FALSE, sep = "\t", check.names = FALSE)
  geneNames = as.vector(rt[, 1])
  uniqGene = unique(geneNames)
  sampleName2 = c(sampleName2, uniqGene)
}
if (length(sampleName1) == 0 || length(sampleName2) == 0) {
  stop("sampleName1 或 sampleName2 为空，请检查 s1.txt 和 s2.txt 文件。")
}
print(getwd())
print(list.files())
file1 = "s1.txt.txt"
if (file.exists(file1)) {
  sampleName1 = readLines(file1)
  print(sampleName1)
} else {
  print(paste("文件不存在：", file1))
}
file2 = "s2.txt.txt"
if (file.exists(file2)) {
  sampleName2 = readLines(file2)
  print(sampleName2)
} else {
  print(paste("文件不存在：", file2))
}
if (length(sampleName1) == 0 || length(sampleName2) == 0) {
  stop("sampleName1 或 sampleName2 为空，请检查 s1.txt.txt 和 s2.txt.txt 文件。")
}
conData = data[, sampleName1, drop = FALSE]
treatData = data[, sampleName2, drop = FALSE]
data = cbind(conData, treatData)
conNum = ncol(conData)
treatNum = ncol(treatData)
Type = c(rep("con", conNum), rep("treat", treatNum))
if (length(unique(Type)) < 2) {
  stop("Type 变量中的水平不足，请检查 conData 和 treatData。")
}
design = model.matrix(~0 + factor(Type))
colnames(design) = c("con", "treat")
fit = lmFit(data, design)
cont.matrix = makeContrasts(treat - con, levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
allDiff = topTable(fit2, adjust = 'fdr', number = 200000)
allDiffOut = rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "all.txt", sep = "\t", quote = FALSE, col.names = FALSE)
outData = rbind(id = paste0(colnames(data), "_", Type), data)
write.table(outData, file = "normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)
logFCfilter = 0.8
P.Value.Filter = 0.05
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & P.Value < P.Value.Filter)), ]
diffSigOut = rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = "diff.txt", sep = "\t", quote = FALSE, col.names = FALSE)
diffGeneExp = data[row.names(diffSig), ]
diffGeneExpOut = rbind(id = paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
write.table(diffGeneExpOut, file = "diffGeneExp.txt", sep = "\t", quote = FALSE, col.names = FALSE)
geneNum = 50
diffSig = diffSig[order(as.numeric(as.vector(diffSig$logFC))), ]
diffGeneName = as.vector(rownames(diffSig))
diffLength = length(diffGeneName)
hmGene = c()
if (diffLength > (2 * geneNum)) {
  hmGene = diffGeneName[c(1:geneNum, (diffLength - geneNum + 1):diffLength)]
} else {
  hmGene = diffGeneName
}
hmExp = data[hmGene, ]
Type = c(rep("Con", conNum), rep("Treat", treatNum))
names(Type) = colnames(data)
Type = as.data.frame(Type)
pdf(file = "heatmap.pdf", width = 10, height = 8)
pheatmap(hmExp,
         annotation = Type,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 9,
         fontsize_row = 9,
         fontsize_col = 9)
dev.off()
