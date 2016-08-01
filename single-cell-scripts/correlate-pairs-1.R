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

sce <- readRDS("hvgs-setup-1.sce.RData")
top.hvgs <- readRDS("top.hvgs.RData")

## ------------------------------------------------------------------------
set.seed(100)
var.cor <- correlatePairs(sce[top.hvgs[1:200],])
write.table(file="hsc-cor.tsv", var.cor, sep="\t", quote=FALSE, col.names=NA)
sink("gene-pairs-1.csv")
head(var.cor)
sink()

## ------------------------------------------------------------------------
sig.cor <- var.cor$FDR <= 0.05
sink("correlated-genes-1.csv")
summary(sig.cor)
sink()

library(RBGL)
g <- ftM2graphNEL(cbind(var.cor$gene1, var.cor$gene2)[sig.cor,], W=NULL, V=NULL, edgemode="undirected")
cl <- highlyConnSG(g)$clusters
cl <- cl[order(lengths(cl), decreasing=TRUE)]
cl <- cl[lengths(cl) > 2]
sink("pair-wise-correlations.csv")
cl
sink()

## ------------------------------------------------------------------------
chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
norm.exprs <- exprs(sce)[chosen,,drop=FALSE]

## ------------------------------------------------------------------------
my.dist <- dist(t(norm.exprs))
my.tree <- hclust(my.dist, method="ward.D2")

## ------------------------------------------------------------------------
my.clusters <- unname(cutree(my.tree, h=50))

## ----heatmaphsc, message=FALSE, fig.width=7, fig.height=7, fig.cap="**Figure 11:** Heatmap of mean-centred normalized log-expression values for correlated HVGs in the HSC data set. Dendrograms are formed by hierarchical clustering on the Euclidean distances between genes (row) or cells (column). Column colours represent the cluster to which each cell is assigned after a tree cut."----
library(gplots)
heat.vals <- norm.exprs - rowMeans(norm.exprs)
clust.col <- rainbow(max(my.clusters))
jpeg("heatmap.jpg")
heatmap.2(heat.vals[chosen,], col=bluered, symbreak=TRUE, trace='none', cexRow=0.8,
          ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree))
dev.off()
## ----pcareduxhsc, fig.width=10, fig.height=6, fig.cap="**Figure 12:** PCA (left) and _t_-SNE plots (right) using only the expression values for significantly correlated HVGs in the HSC data set. Cells are coloured according to the level of _H2-Aa_ expression."----
out.pca <- plotPCA(sce, exprs_values="exprs", feature_set=chosen, colour_by="H2-Aa") + fontsize
set.seed(100)
out.tsne <- plotTSNE(sce, exprs_values="exprs", feature_set=chosen, colour_by="H2-Aa") + fontsize
jpeg("pca-tSNE-plots.jpg")
multiplot(out.pca, out.tsne, cols=2)
dev.off()

## ------------------------------------------------------------------------
saveRDS(file="hsc-data.rds", sce)
