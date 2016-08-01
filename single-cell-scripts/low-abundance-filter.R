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

keep <- rowMeans(counts(sce)) >= 1
sce <- sce[keep,] 
sink(args[2])
sum(keep)
sink()

if(length(args) == 3){
  sce <- sce[!fData(sce)$is_feature_control_Mt,]
  saveRDS(sce, args[3])
} else{
  saveRDS(sce, args[3])
  ## ------------------------------------------------------------------------
  alt.keep <- rowSums(is_exprs(sce)) >= 10
  sink(args[4])
  sum(alt.keep)
  sink()
  
  ## ----geneplothsc, Frequency of expression against the mean expression for each gene. Circles represent endogenous genes and triangles represent spike-in transcripts or mitochondrial genes. The bars on each axis represent the location of each gene on that axis. Genes with expression frequencies higher than the dropout rate are defined as those above a non-linear trend fitted to the spike-in transcripts."----
  jpeg(args[5])
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  plotQC(sce, type = "exprs-freq-vs-mean") + fontsize
  dev.off()
}


