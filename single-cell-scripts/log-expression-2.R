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

args <- commandArgs(TRUE)

var.out <- readRDS("var.out-2.RData")
var.fit <- readRDS("var.fit-2.RData")
sce <- readRDS(args[1])

## ----hvgplotbrain, Variance of normalized log-expression values for each gene in the brain data set, plotted against the mean log-expression. The red line represents the mean-   dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points)."----
jpeg("log-expression-2.jpg")
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
points(var.fit$mean, var.fit$var, col="red", pch=16)
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)
dev.off()

## ------------------------------------------------------------------------
top.hvgs <- order(var.out$bio, decreasing=TRUE)
write.table(file="brain-hvg.tsv", var.out[top.hvgs,], sep="\t", quote=FALSE, col.names=NA)
sink("top-hvgs-2.csv")
head(var.out[top.hvgs,])
sink()

## ----hvgboxplotbrain, fig.cap="**Figure 21:** Boxplots of normalized log-expression values for the top 10 HVGs in the brain data set. Points correspond to cells that are more than 1.5 interquartile ranges from the edge of each box."----
examined <- top.hvgs[1:10]
all.names <- matrix(rownames(sce)[examined], nrow=length(examined), ncol=ncol(sce))
jpeg("normalized-log-expression-2.jpg")
boxplot(split(exprs(sce)[examined,], all.names), las=2, ylab="Normalized log-expression", col="grey80")
dev.off()

## ------------------------------------------------------------------------
set.seed(100)
var.cor <- correlatePairs(sce[top.hvgs[1:200],], design=design)
write.table(file="brain-cor.tsv", var.cor, sep="\t", quote=FALSE, col.names=NA)
sink("gene-pairs-2.csv")
head(var.cor)
sink()

sig.cor <- var.cor$FDR <= 0.05
sink("correlated-genes-2.csv")
summary(sig.cor)
sink()

saveRDS(var.cor, "var.cor.RData")
saveRDS(sig.cor, "sig.cor.RData")
