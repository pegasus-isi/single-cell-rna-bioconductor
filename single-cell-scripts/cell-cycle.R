#!/usr/bin/Rscript

library(scater)
library(gplots)
library(scran)
library(org.Mm.eg.db)

args <- commandArgs(TRUE)
sce <- readRDS(args[1])

## ----phaseplothsc, message=FALSE, fig.cap="**Figure 8:** Cell cycle phase scores from applying the pair-based classifier on the HSC data set, where each point represents a cell."----
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
keep <- !is.na(ensembl)
assignments <- cyclone(sce[keep,], mm.pairs, gene.names=ensembl[keep])
jpeg(args[2])
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()

