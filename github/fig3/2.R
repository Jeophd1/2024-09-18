rm(list = ls())
library(e1071)
library(kernlab)
library(caret)
library(pROC)
set.seed(123)
setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3C\\回复审稿人")
inputFile = "diffGeneExp.txt"
data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data = t(data)
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))
group <- factor(group, levels = c("con", "treat"))
rfe_ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
Profile <- rfe(
  x = data,
  y = group,
  sizes = c(2, 4, 6, 8, seq(10, 40, by = 3)),
  rfeControl = rfe_ctrl,
  method = "svmRadial",
  trControl = trainControl(classProbs = TRUE),
  tuneLength = 5,
  verbose = FALSE
)
print(Profile$optVariables)
selected_features <- Profile$optVariables
data_selected <- data[, selected_features]
trctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE)
svm_model <- train(
  x = data_selected,
  y = group,
  method = "svmRadial",
  trControl = trctrl
)
pred <- predict(svm_model, data_selected, type = "prob")
roc_curve <- roc(response = group, predictor = pred[, "treat"], levels = c("con", "treat"))
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))
pdf(file = "ROC_SVM-RFE.pdf", width = 6, height = 5.5)
plot(roc_curve, main = "ROC Curve for SVM-RFE")
dev.off()
threshold = 0.5
pred_class <- ifelse(pred[, "treat"] > threshold, "treat", "con")
conf_matrix <- table(True = group, Predicted = pred_class)
TP <- conf_matrix["treat", "treat"]
TN <- conf_matrix["con", "con"]
FP <- conf_matrix["con", "treat"]
FN <- conf_matrix["treat", "con"]
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
accuracy <- (TP + TN) / sum(conf_matrix)
print("Confusion Matrix:")
print(conf_matrix)
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
print(paste("Accuracy:", accuracy))
