rm(list = ls())

install.packages("glmnet")
install.packages("pROC")
library(glmnet)
library(pROC)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig4\\训练")

expFile = "normalized_data.txt"
geneFile = "interGenes.txt"

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

y <- gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y <- ifelse(y == "con", 0, 1)

geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)

bioCol = rainbow(nrow(geneRT), s = 0.9, v = 0.9)

aucText = c()
k = 0

print(unique(y))

for (x in as.vector(geneRT[, 1])) {
  k = k + 1
  roc1 = roc(y, as.numeric(rt[x, ]))
  if (k == 1) {
    pdf(file = "ROC.genes.pdf", width = 5, height = 4.75)
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "")
    aucText = c(aucText, paste0(x, ", AUC = ", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", add = TRUE)
    aucText = c(aucText, paste0(x, ", AUC = ", sprintf("%.3f", roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd = 2, bty = "n", col = bioCol[1:k])

dev.off()
