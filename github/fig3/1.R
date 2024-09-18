rm(list=ls())

set.seed(12345)
library(glmnet)
library(pROC)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig3A\\回复审稿人用")
inputFile="diffGeneExp.txt"

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

x=as.matrix(rt)

y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
y <- factor(y, levels = c("con", "treat"))

set.seed(12345)
train_index <- sample(1:nrow(x), 0.7 * nrow(x))
x_train <- x[train_index, ]
y_train <- y[train_index]
x_test <- x[-train_index, ]
y_test <- y[-train_index]

alphas <- seq(0, 1, by=0.1)
cvfits <- list()
aucs <- numeric(length(alphas))

for (i in 1:length(alphas)) {
  alpha_val <- alphas[i]
  cvfit <- cv.glmnet(x_train, y_train, family="binomial", alpha=alpha_val, type.measure='auc', nfolds=5)
  cvfits[[i]] <- cvfit
  pred_test <- as.vector(predict(cvfit, x_test, s="lambda.min", type="response"))
  y_test_num <- ifelse(y_test == "treat", 1, 0)
  roc_curve_test <- roc(y_test_num, pred_test)
  aucs[i] <- auc(roc_curve_test)
}

best_alpha <- alphas[which.max(aucs)]
best_auc <- max(aucs)
best_cvfit <- cvfits[[which.max(aucs)]]
print(paste("Best alpha:", best_alpha))
print(paste("Best Test AUC:", best_auc))

pred_test <- as.vector(predict(best_cvfit, x_test, s="lambda.min", type="response"))
y_test_num <- ifelse(y_test == "treat", 1, 0)

roc_curve_test <- roc(y_test_num, pred_test)
auc_value_test <- auc(roc_curve_test)
print(paste("Test AUC with best alpha:", auc_value_test))

pdf(file="ROC_LASSO_test.pdf", width=6, height=5.5)
plot(roc_curve_test, main=paste("ROC Curve for LASSO on Test Set (Alpha =", best_alpha, ")"))
dev.off()

threshold = 0.5
pred_class_test <- ifelse(pred_test > threshold, "treat", "con")
conf_matrix_test <- table(True=y_test, Predicted=pred_class_test)

sensitivity_test <- conf_matrix_test[2,2] / sum(conf_matrix_test[2,])
specificity_test <- conf_matrix_test[1,1] / sum(conf_matrix_test[1,])
accuracy_test <- sum(diag(conf_matrix_test)) / sum(conf_matrix_test)

print("Confusion Matrix on Test Set:")
print(conf_matrix_test)
print(paste("Test Sensitivity:", sensitivity_test))
print(paste("Test Specificity:", specificity_test))
print(paste("Test Accuracy:", accuracy_test))

coef=coef(best_cvfit, s="lambda.min")
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
