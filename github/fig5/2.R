rm(list = ls())

if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(glmnet)
library(pROC)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig4\\训练")

expFile = "normalized_data.txt"
geneFile = "interGenes.txt"

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y == "con", 0, 1)

geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)

bioCol = rainbow(nrow(geneRT), s = 0.9, v = 0.9)

aucText = c()
k = 0

rt = rt[as.vector(geneRT[, 1]), ]
rt = as.data.frame(t(rt))

logit = glm(y ~ ., family = binomial(link = 'logit'), data = rt)

pred = predict(logit, type = "response")

roc1 = roc(y, pred)

ci1 = ci.auc(roc1, method = "bootstrap")

ciVec = as.numeric(ci1)

pdf(file = "ROC.model.pdf", width = 5, height = 4.75)

plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = "Model")

text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")

dev.off()
