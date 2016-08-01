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

# gunzip("GSE61533_HTSEQ_count_results.xls.gz", remove=FALSE, overwrite=TRUE)


all.counts <- read.csv('GSE61533_HTSEQ_count_results.csv', header=TRUE, row.names=1)

sce <- newSCESet(countData=all.counts)
sink("count.csv")
dim(sce)
sink()

is.spike <- grepl("^ERCC", rownames(sce))
isSpike(sce) <- is.spike
is.mito <- grepl("^mt-", rownames(sce))
sce <- calculateQCMetrics(sce, feature_controls=list(Spike=is.spike, Mt=is.mito))
saveRDS(sce, "sce.RData")

sink("column-names.txt")
head(colnames(pData(sce)))
sink()


