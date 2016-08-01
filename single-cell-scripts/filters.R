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
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

args <- commandArgs(TRUE)
sce <- readRDS(args[1])

libsize.drop <- isOutlier(sce$total_counts, n=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, n=3, type="lower")
mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, n=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_feature_controls_Spike, n=3, type="higher")

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
saveRDS(sce, args[2])

write.table(data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce)), 
           args[3])

if(length(args) == 4){
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  jpeg(args[4])
  plotPCA(sce, pca_data_input="pdata") + fontsize
  dev.off()
}

