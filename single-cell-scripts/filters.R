#!/usr/bin/Rscript

library(scran)
library(edgeR)
library(scater)
library(Rtsne)
library(mvoutlier)
library(destiny)
library(gplots)
library(R.utils)


args <- commandArgs(TRUE)
sce <- readRDS("sce.RData")
libsize.drop <- isOutlier(sce$total_counts, n=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, n=3, type="lower")
mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, n=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_feature_controls_Spike, n=3, type="higher")

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
saveRDS(sce, "sce-filters-1.RData")

write.table(data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
                       ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce)), 
            "filtered-cells-1.csv")

if(args == 4){
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  jpeg(file = "pca-filter.jpg")
  print(plotPCA(sce, pca_data_input="pdata") + fontsize)
  dev.off()
}



