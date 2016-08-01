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

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- normalize(sce)

jpeg(args[2])
plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")
dev.off()

saveRDS(sce, args[3])