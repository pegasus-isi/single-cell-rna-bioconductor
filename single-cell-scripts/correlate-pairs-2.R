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

var.cor <- readRDS("var.cor.RData")
sig.cor <- readRDS("sig.cor.RData")
sce <- readRDS(args[1])

chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
norm.exprs <- exprs(sce)[chosen,,drop=FALSE]
library(limma)
norm.exprs <- removeBatchEffect(norm.exprs, batch=sce$tissue)

## ----heatmapbrain, message=FALSE, fig.width=7, fig.height=10, fig.cap="**Figure 22:** Heatmap of mean-centred normalized log-expression values for correlated HVGs in the brain data set. Dendrograms are formed by hierarchical clustering on the Euclidean distances between genes (row) or cells (column). Column colours represent the cluster to which each cell is assigned after a dynamic tree cut."----
my.dist <- dist(t(norm.exprs))
my.tree <- hclust(my.dist, method="ward.D2")
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
heat.vals <- norm.exprs - rowMeans(norm.exprs)
clust.col <- rainbow(max(my.clusters))

pdf("brain.heat.pdf", width=7, height=20)
heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.8,
          ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree))
dev.off()

cluster <- factor(my.clusters)
design <- model.matrix(~0 + cluster + sce$tissue)

## ------------------------------------------------------------------------
y <- convertTo(sce)

## ------------------------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)

sink("dispersion-summary.csv")
summary(y$tagwise.dispersion)
sink()

saveRDS(cluster, "cluster.RData")
saveRDS(design, "design.RData")
saveRDS(y, "y.RData")
saveRDS(fit, "fit.RData")
saveRDS(my.tree, "my.tree.RData")
saveRDS(my.clusters, "my.clusters.RData")
saveRDS(clust.col, "clust.col.RData")