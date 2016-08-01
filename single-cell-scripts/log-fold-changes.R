#!/usr/bin/Rscript

cluster <- readRDS("cluster.RData")
design <- readRDS("design.RData")
fit <- readRDS("fit.RData")
y <- readRDS("y.RData")

## ------------------------------------------------------------------------
# Automatic construction of the contrast matrix.
nclusters <- nlevels(cluster)
contrast.matrix <- matrix(0, ncol(design), nclusters) 
contrast.matrix[1,] <- -1 
diag(contrast.matrix) <- 1
contrast.matrix <- contrast.matrix[,-1]
res.any <- glmLRT(fit, contrast=contrast.matrix)

# Computing log-fold changes between each cluster and the average of the rest.
cluster.expression <- fit$coefficients[,seq_len(nclusters)] 
other.expression <- (rowSums(cluster.expression) - cluster.expression)/(nclusters-1)
log.fold.changes <- cluster.expression - other.expression
colnames(log.fold.changes) <- paste0("LogFC.for.", levels(cluster))
rownames(log.fold.changes) <- NULL

# Ordering by the likelihood ratio; p-values affected by numerical imprecision.
any.de <- data.frame(Gene=rownames(y), log.fold.changes, 
                     LR=res.any$table$LR, stringsAsFactors=FALSE)
any.de <- any.de[order(any.de$LR, decreasing=TRUE),]
sink("DE-genes.csv")
head(any.de)
sink()
