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

## ----pcaplotbrain, PCA plots constructed from the normalized expression values for all remaining cells in the brain data set. Left: cells are coloured according to the tissue of origin (cortex or hippocampus). Right: cells are coloured according to the sex of the mouse -- male (-1), female (1) or unassigned (0)."----
pca1 <- plotPCA(sce, exprs_values="exprs", colour_by="tissue") + fontsize
pca2 <- plotPCA(sce, exprs_values="exprs", colour_by="sex") + fontsize
jpeg("pca-tissue-sex.jpg")
multiplot(pca1, pca2, cols=2)
dev.off()

## ----tsneplotbrain, _t_-SNE plots constructed from the normalized expression values for all remaining cells in the brain data set. Left: cells are coloured according to the tissue of origin (cortex or hippocampus). Right: cells are coloured according to the sex of the mouse -- male (-1), female (1) or unassigned (0)."----
set.seed(100)
tsne1 <- plotTSNE(sce, exprs_values="exprs", colour_by="tissue") + fontsize
set.seed(100)
tsne2 <- plotTSNE(sce, exprs_values="exprs", colour_by="sex") + fontsize
jpeg("tSNE-tissue-sex.jpg")
multiplot(tsne1, tsne2, cols=2)
dev.off()

## ----pca2plotbrain, PCA plots constructed from the normalized expression values for all cells in the brain data set from the cortex (left) or hippocampus (right). Each cell is coloured according to the C1 chip on which its library was prepared."----
sce$chip <- sub("_.*", "", sce$cell_id)
saveRDS(sce, args[2])
pca1 <- plotPCA(sce[,sce$tissue=="sscortex"], exprs_values="exprs", 
                colour_by="chip", legend="none") + fontsize + ggtitle("Cortex")
pca2 <- plotPCA(sce[,sce$tissue!="sscortex"], exprs_values="exprs", 
                colour_by="chip", legend="none") + fontsize + ggtitle("Hippocampus")

jpeg("C1-chip-colored.jpg")
multiplot(pca1, pca2, cols=2)
dev.off()
## ------------------------------------------------------------------------
design <- model.matrix(~sce$tissue)

## ------------------------------------------------------------------------
var.fit <- trendVar(sce, trend="loess", design=design, span=0.4)
var.out <- decomposeVar(sce, var.fit)

saveRDS("var.fit-2.RData")
saveRDS("var.out-2.RData")


