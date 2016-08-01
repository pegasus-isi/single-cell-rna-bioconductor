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

sce <- readRDS("normalize-1.sce.RData")

g1.only <- assignments$score$G1 > 0.5
sce <- sce[,g1.only]
saveRDS(sce, "hvgs-setup-1.sce.RData")

## ------------------------------------------------------------------------
var.fit <- trendVar(sce, trend="loess", span=0.3)
saveRDS(var.fit, "var.fit.1.RData")

## ------------------------------------------------------------------------
var.out <- decomposeVar(sce, var.fit)
saveRDS(var.out, "var.out.1.RData")

## ------------------------------------------------------------------------
var.fit2 <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
var.out2 <- decomposeVar(sce, var.fit2)
saveRDS(var.out2, "var.out2.1.RData")