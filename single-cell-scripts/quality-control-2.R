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
sce <- readRDS(args[1])
               
par(mfrow=c(1,2))
jpeg(args[2])
hist(sce$pct_counts_feature_controls_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()

jpeg(args[3])
hist(sce$pct_counts_feature_controls_Spike, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
dev.off()