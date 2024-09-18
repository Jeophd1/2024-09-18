rm(list = ls())

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(dplyr)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig1c")
expFile = "diffGeneExp.txt"

rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
data = data[rowMeans(data) > 0, ]
exp = data

Type = gsub("(.*)\\_(.*)", "\\2", colnames(data))

exp = as.data.frame(t(exp))
exp = cbind(exp, Type = Type)
data = melt(exp, id.vars = c("Type"))
colnames(data) = c("Type", "Gene", "Expression")

p = ggboxplot(data, x = "Gene", y = "Expression", color = "Type",
              xlab = "",
              ylab = "Gene expression",
              legend.title = "Type",
              palette = c("blue", "red"),
              width = 0.5)
p = p + rotate_x_text(60)
p1 = p + stat_compare_means(aes(group = Type),
                            method = "wilcox.test",
                            symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05), symbols = c("***", "**", "*")),
                            label = "p.signif")

pdf(file = "boxplot.pdf", width = 10, height = 5)
print(p1)
dev.off()

test_results <- data %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Expression ~ Type)$p.value)

print(test_results)

out_of_bounds <- test_results %>%
  filter(p_value < 0 | p_value > 0.05)

print(out_of_bounds)
