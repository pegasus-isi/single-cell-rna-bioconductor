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
library(TxDb.Mmusculus.UCSC.mm10.ensGene

sce <- readRDS("normalize-1.sce.RData")

## ----pcaplothsc,PCA plot constructed from normalized log-expression values, where each point represents a cell in the HSC data set. First and second components are shown, along with the percentage of variance explained. Bars represent the coordinates of the cells on each axis. None of the cells are controls (e.g., empty wells) so the legend can be ignored."----
jpeg("pca-normalized.jpg")
plotPCA(sce, exprs_values="exprs") + fontsize
dev.off()

## ----tsneplothsc,_t_-SNE plot constructed from normalized log-expression values using a range of perplexity values. In each plot, each point represents a cell in the HSC data set. Bars represent the coordinates of the cells on each axis.", fig.width=12, fig.height=6----
set.seed(100)
out5 <- plotTSNE(sce, exprs_values="exprs", perplexity=5) + fontsize + ggtitle("Perplexity = 5")
out10 <- plotTSNE(sce, exprs_values="exprs", perplexity=10) + fontsize + ggtitle("Perplexity = 10")
out20 <- plotTSNE(sce, exprs_values="exprs", perplexity=20) + fontsize + ggtitle("Perplexity = 20")
jpeg("tSNE-plots.jpg")
multiplot(out5, out10, out20, cols=3)
dev.off()

