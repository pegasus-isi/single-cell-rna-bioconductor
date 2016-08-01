#!/usr/bin/Rscript

args <- commandArgs(TRUE)

cluster <- readRDS("cluster.RData")
design <- readRDS("design.RData")
y <- readRDS("y.RData")
fit <- readRDS("fit.RData")
my.tree <- readRDS("my.tree.RData")
my.clusters <- readRDS("my.clusters.RData")
clust.col <- readRDS("clust.col.RData")
sce <- args[1]

result.logFC <- result.PValue <- list()
chosen.clust <- which(levels(cluster)=="2") # character, as 'cluster' is a factor.
for (clust in seq_len(nlevels(cluster))) {
  if (clust==chosen.clust) { next }
  contrast <- numeric(ncol(design))
  contrast[chosen.clust] <- 1
  contrast[clust] <- -1
  res <- glmLRT(fit, contrast=contrast)
  con.name <- paste0('vs.', levels(cluster)[clust])
  result.logFC[[con.name]] <- res$table$logFC
  result.PValue[[con.name]] <- res$table$PValue
}

## ------------------------------------------------------------------------
max.PValue <- do.call(pmax, result.PValue)
all.logFC <- do.call(cbind, result.logFC)
all.signs <- sign(all.logFC)
same.sign <- rowSums(all.signs[,1]!=all.signs)==0L
marker.set <- data.frame(Gene=rownames(y), logFC=all.logFC, 
                         PValue=max.PValue, stringsAsFactors=FALSE)
marker.set <- marker.set[same.sign,]
marker.set <- marker.set[order(marker.set$PValue),]
sink("markers.csv")
head(marker.set)
sink()

## ----heatmapmarkerbrain, Heatmap of mean-centred normalized log-expression values for the top set of markers for cluster 2 in the brain data set. Column colours represent the cluster to which each cell is assigned."----
write.table(marker.set, file="brain-marker-2.tsv", sep="\t", quote=FALSE, row.names=FALSE)
top.markers <- marker.set$Gene[1:20]
norm.exprs <- exprs(sce)[top.markers,,drop=FALSE]
heat.vals <- norm.exprs - rowMeans(norm.exprs)
jpeg("cluster2-heatmap.jpg")
heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=1,
          ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), dendrogram='none')
legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters)), pch=16)
dev.off()