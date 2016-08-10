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
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

sce <- computeSumFactors(sce, sizes=c(20, 40, 60, 80))
sink("size-factors.txt")
summary(sizeFactors(sce))
sink()

## ----normplothsc, fig.cap="**Figure 5:** Size factors from deconvolution, plotted against library sizes for all cells in the HSC data set. Axes are shown on a log-scale."----
jpeg("libSize-vs-sizeFactors-1.jpg")
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
     ylab="Library size (millions)", xlab="Size factor")
dev.off()

## ------------------------------------------------------------------------
sce <- normalize(sce)
saveRDS(sce, "normalize-1.sce.RData")
saveRDS(fontsize, "normalize-1.fontsize.RData")

