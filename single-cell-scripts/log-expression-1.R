#!/usr/bin/Rscript

library(scran)
library(edgeR)
library(scater)
library(Rtsne)
library(mvoutlier)
library(destiny)
library(gplots)
library(gdata)
library(openxlsx)
library(R.utils)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

var.fit <- readRDS("var.fit.1.RData")
var.out <- readRDS("var.out.1.RData")
var.out2 <- readRDS("var.out2.1.RData")
sce <- readRDS("hvgs-setup-1.sce.RData")

## ----hvgplothsc, Variance of normalized log-expression values for each gene in the HSC data set, plotted against the mean log-expression. The red line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points). The blue line represents the trend fitted to the variances of the endogenous genes."----
jpeg("log-expression-1.jpg")
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
points(var.fit$mean, var.fit$var, col="red", pch=16)
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)
lines(var.out2$mean[o], var.out2$tech[o], col="dodgerblue", lwd=2)
dev.off()

top.hvgs <- order(var.out2$bio, decreasing=TRUE)
saveRDS(top.hvgs, "top-hvgs.RData")
write.table(file="hsc-hvg.tsv", var.out2[top.hvgs,], sep="\t", quote=FALSE, col.names=NA)
sink("top-hvgs.csv")
head(var.out2[top.hvgs,])
sink()

examined <- top.hvgs[1:10]
all.names <- matrix(rownames(sce)[examined], nrow=length(examined), ncol=ncol(sce))
jpeg("normalized-log-expression.jpg")
boxplot(split(exprs(sce)[examined,], all.names), las=2, ylab="Normalized log-expression", col="grey80")
dev.off()


