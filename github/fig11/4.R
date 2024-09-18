rm(list = ls())

# 安装并加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("limma")

required_packages <- c("ggpubr", "reshape2", "RColorBrewer")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig8\\fig8H") 

ACECluFile <- "ACEcluster.txt"    
geneCluFile <- "geneCluster.txt"     
scoreFile <- "ACEscore.txt"          

ACEClu <- read.table(ACECluFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
geneClu <- read.table(geneCluFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
score <- read.table(scoreFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

twoCluster <- cbind(ACEClu, geneClu)

sameSample <- intersect(row.names(twoCluster), row.names(score))
data <- cbind(score[sameSample, , drop = FALSE], twoCluster[sameSample, c("ACEcluster", "geneCluster"), drop = FALSE])

# 确保 'geneCluster' 列存在
if (!"geneCluster" %in% colnames(data)) {
  stop("geneCluster 列在聚类数据中不存在。请检查文件内容。")
}

data$ACEcluster <- factor(data$ACEcluster, levels = levels(factor(data$ACEcluster)))
group <- levels(factor(data$ACEcluster))
comp <- combn(group, 2)
my_comparisons <- list()
for(i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[,i]
}

bioCol <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol <- bioCol[1:length(levels(factor(data$ACEcluster)))]

# 确保 Fraction 和 ACEscore 为数值型
data$Fraction <- as.numeric(as.character(data$Fraction))
data$ACEscore <- as.numeric(as.character(data$ACEscore))
if (any(is.na(data$Fraction))) {
  warning("存在 NA 值在 Fraction 列，请检查数据源。")
}
if (any(is.na(data$ACEscore))) {
  warning("存在 NA 值在 ACEscore 列，请检查数据源。")
}

# 绘制 ACEcluster 的箱线图
boxplot_ACE <- ggboxplot(data, x = "ACEcluster", y = "ACEscore", color = "ACEcluster", 
                         xlab = "ACEcluster",
                         ylab = "ACEscore",
                         legend.title = "ACEcluster",
                         palette = bioCol,
                         add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

ggsave("ACEcluster.pdf", plot = boxplot_ACE, width = 5, height = 4.5)

# 准备 geneCluster 的比较
data$geneCluster <- factor(data$geneCluster, levels = levels(factor(data$geneCluster)))
group_gene <- levels(factor(data$geneCluster))
comp_gene <- combn(group_gene, 2)
my_comparisons_gene <- list()
for(i in 1:ncol(comp_gene)) {
  my_comparisons_gene[[i]] <- comp_gene[,i]
}

bioCol_gene <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol_gene <- bioCol_gene[1:length(levels(factor(data$geneCluster)))]

# 绘制 geneCluster 的箱线图
boxplot_gene <- ggboxplot(data, x = "geneCluster", y = "ACEscore", color = "geneCluster", 
                          xlab = "geneCluster",
                          ylab = "ACEscore",
                          legend.title = "geneCluster",
                          palette = bioCol_gene,
                          add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons_gene) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

ggsave("geneCluster.pdf", plot = boxplot_gene, width = 5, height = 4.5)
