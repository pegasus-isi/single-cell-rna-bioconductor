#!/usr/bin/Rscript

library(scran)
library(edgeR)
library(scater)
library(Rtsne)
library(mvoutlier)
library(destiny)
library(gplots)
library(R.utils)
library(org.Mm.eg.db)
# library(TxDb.Mmusculus.UCSC.mm10.ensGene)

## Fix arguments 
args <- commandArgs(TRUE)

readFormat <- function(infile) { 
  # First column is empty.
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  # First column after row names is some useless filler.
  counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] 
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}

## ------------------------------------------------------------------------
endo.data <- readFormat("expression_mRNA_17-Aug-2014.txt")
spike.data <- readFormat("expression_spikes_17-Aug-2014.txt")
mito.data <- readFormat("expression_mito_17-Aug-2014.txt")

m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]
# 
# ## ---- echo=FALSE---------------------------------------------------------
stopifnot(identical(endo.data$metadata$cell_id, spike.data$metadata$cell_id)) # should be the same.
stopifnot(all(endo.data$metadata$cell_id==mito.data$metadata$cell_id)) # should now be the same.

## ------------------------------------------------------------------------
# all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
# metadata <- AnnotatedDataFrame(endo.data$metadata)
# sce <- newSCESet(countData=all.counts, phenoData=metadata)
# sink("sce-countLoading-dims.csv")
# dim(sce)
# sink()

# ## ------------------------------------------------------------------------
# nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
# is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
# isSpike(sce) <- is.spike
# is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
# 
# ## ------------------------------------------------------------------------
# sce <- calculateQCMetrics(sce, feature_controls=list(Spike=is.spike, Mt=is.mito)) 
# 
# saveRDS(sce, args[2])