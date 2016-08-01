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

sce <- readRDS("hsc-data.rds")

deconv.sf <- sizeFactors(sce)
sce <- computeSpikeFactors(sce)

## ----normplotspikehsc, Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the HSC data set. Axes are shown on a log-scale."----
jpeg("normalization-vs-deconvolution.jpg")
plot(sizeFactors(sce), deconv.sf, pch=16, log="xy", xlab="Size factor (spike-in)",
     ylab="Size factor (deconvolution)")
dev.off()

## ------------------------------------------------------------------------
incoming <- read.csv("nbt.3102-S7.xlsx", sheet=1, rowNames=TRUE)
incoming <- incoming[,!duplicated(colnames(incoming))] # Remove duplicated genes.
sce <- newSCESet(exprsData=t(incoming))

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(org.Mm.eg.db)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
set.seed(100)

## ----phaseplotth2, message=FALSE, fig.cap="**Figure 25:** Cell cycle phase scores from applying the pair-based classifier on the T~H~2 data set, where each point represents a cell."----
anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
keep <- !is.na(ensembl) # Remove genes without ENSEMBL IDs.
assignments <- cyclone(exprs(sce)[keep,], mm.pairs, gene.names=ensembl[keep])
jpeg("phase-scores.jpg")
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()
## ------------------------------------------------------------------------
design <- model.matrix(~ G1 + G2M, assignments$score)
fit.block <- trendVar(sce, use.spikes=NA, trend="loess", design=design)
dec.block <- decomposeVar(sce, fit.block)

## ----pcaplotth2, PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 data set. Each point represents a cell, coloured according to its G1 score. Only the top 500 HVGs were used to make each PCA plot."----
fit <- trendVar(sce, use.spikes=NA, trend="loess")
dec <- decomposeVar(sce, fit)
top.hvgs <- order(dec$bio, decreasing=TRUE)[1:500]
sce$G1score <- assignments$score$G1
out <- plotPCA(sce, select=top.hvgs, colour_by="G1score") + fontsize + ggtitle("Before removal")

top.hvgs2 <- order(dec.block$bio, decreasing=TRUE)[1:500]
corrected <- removeBatchEffect(exprs(sce), covariates=assignments$score[,c("G1", "G2M")])
sce2 <- newSCESet(exprsData=corrected, phenoData=phenoData(sce))
out2 <- plotPCA(sce2, select=top.hvgs2, colour_by="G1score") + fontsize + ggtitle("After removal")
jpeg("cell-cycle-removal.jpg")
multiplot(out, out2, cols=2)
dev.off()

## ----diffusionth2,A diffusion map for the T~H~2 data set, where each cell is coloured by its expression of _Gata3_."----
jpeg("diffusion-map.jpg")
plotDiffusionMap(sce2, colour_by="Gata3") + fontsize
dev.off()

saveRDS(sce, args[1])