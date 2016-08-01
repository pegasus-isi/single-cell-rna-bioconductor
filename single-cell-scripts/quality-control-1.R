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

input <- as.numeric(args[2])

par(mfrow=c(1,2))
jpeg(args[3])
hist(sce$total_counts/input, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells");
dev.off()

jpeg(args[4])
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells");

dev.off()