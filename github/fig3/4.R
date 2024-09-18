rm(list = ls())

library(e1071)
library(kernlab)
library(caret)
library(fastshap)
library(ggplot2)

set.seed(123)
setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3E")
inputFile <- "C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3E\\diffGeneExp.txt"

data <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data <- t(data)
group <- gsub("(.*)\\_(.*)", "\\2", rownames(data))

Profile <- rfe(
  x = data,
  y = as.numeric(as.factor(group)),
  sizes = c(2, 4, 6, 8, seq(10, 40, by = 3)),
  rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
  method = "svmRadial"
)

pdf(file = "SVM-RFE.pdf", width = 6, height = 5.5)
par(las = 1)
x <- Profile$results$Variables
y <- Profile$results$RMSE
plot(x, y, xlab = "Variables", ylab = "RMSE (Cross-Validation)", col = "darkgreen")
lines(x, y, col = "darkgreen")
wmin <- which.min(y)
wmin.x <- x[wmin]
wmin.y <- y[wmin]
points(wmin.x, wmin.y, col = "blue", pch = 16)
text(wmin.x, wmin.y, paste0('N=', wmin.x), pos = 2, col = 2)
dev.off()

featureGenes <- Profile$optVariables
write.table(file = "SVM-RFE.gene.txt", featureGenes, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

pred_fun <- function(object, newdata) {
  predict(object, newdata)
}

shap_values <- explain(
  Profile$fit,
  newdata = data,
  pred_wrapper = pred_fun,
  nsim = 100
)

plot_shap_for_feature <- function(shap_matrix, feature, data) {
  if (!(feature %in% colnames(data))) {
    stop("Feature not found in data.")
  }
  shap_feature_values <- shap_matrix[, feature]
  plot(
    data[, feature],
    shap_feature_values,
    xlab = feature,
    ylab = "SHAP value",
    main = paste("SHAP Dependence Plot for", feature)
  )
}

features_to_plot <- c("MTA2", "H1-2", "H2AC20", "BRD7", "ACTL6A", "TFPI2", "CACNA1G", "INPP5D", "DSCC1", "ZNF615", "MIR34A")

for (feature in features_to_plot) {
  plot_shap_for_feature(shap_values, feature, data)
}

mean_abs_shap_values <- apply(shap_values[, features_to_plot, drop = FALSE], 2, function(x) mean(abs(x)))
mean_abs_shap_df <- data.frame(Feature = names(mean_abs_shap_values), MeanAbsShapValue = mean_abs_shap_values)
mean_abs_shap_df <- mean_abs_shap_df[order(mean_abs_shap_df$MeanAbsShapValue, decreasing = TRUE), ]

p <- ggplot(mean_abs_shap_df, aes(x = reorder(Feature, MeanAbsShapValue), y = MeanAbsShapValue)) +
  geom_bar(stat = "identity", fill = "#0073C2FF") +
  coord_flip() +
  xlab("Feature") +
  ylab("Mean Absolute SHAP Value") +
  theme_minimal()

ggsave(
  filename = "barplot_shap_values.pdf",
  plot = p,
  width = 8,
  height = 6,
  path = "C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3E"
)
