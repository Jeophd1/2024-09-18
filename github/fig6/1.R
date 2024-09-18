rm(list = ls())
install.packages("rms")
install.packages("rmda")
library(rms)
library(rmda)

setwd("C:\\Users\\JeoPhD\\Desktop\\KOA2\\fig5")
inputFile = "diffGeneExp.txt"

data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data = t(data)
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt = cbind(as.data.frame(data), Type = group)
paste(colnames(data), collapse = "+")
ddist = datadist(rt)
options(datadist = "ddist")

lrmModel = lrm(Type ~ BRD7 + ZNF615, data = rt, x = TRUE, y = TRUE)
nomo = nomogram(lrmModel, fun = plogis,
                fun.at = c(0.0001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
                lp = FALSE, funlabel = "Risk of Disease")
pdf("Nom.pdf", width = 8, height = 6)
plot(nomo)
dev.off()

cali = calibrate(lrmModel, method = "boot", B = 1000)
pdf("Calibration.pdf", width = 6, height = 6)
plot(cali,
     xlab = "Predicted probability",
     ylab = "Actual probability", sub = FALSE)
dev.off()

rt$Type = ifelse(rt$Type == "con", 0, 1)
dc = decision_curve(Type ~ BRD7 + ZNF615, data = rt,
                    family = binomial(link = 'logit'),
                    thresholds = seq(0, 1, by = 0.01),
                    confidence.intervals = 0.95)
pdf(file = "DCA.pdf", width = 6, height = 8)
plot_decision_curve(dc,
                    curve.names = "ACE genes",
                    xlab = "Threshold probability",
                    cost.benefit.axis = TRUE,
                    col = "red",
                    confidence.intervals = FALSE,
                    standardize = FALSE)
dev.off()

pdf(file = "clinical_impact.pdf", width = 6, height = 8)
plot_clinical_impact(dc,
                     confidence.intervals = TRUE,
                     col = c("red", "blue"))
dev.off()
