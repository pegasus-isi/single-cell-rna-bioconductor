#!/usr/bin/Rscript

# API Documentation: http://pegasus.isi.edu/documentation

if (is.element("dax3", installed.packages()[,1])) {
  require(dax3)
} else {
  stop("The R DAX Generator API is not installed.\nPlease refer to 'https://pegasus.isi.edu/documentation/dax_generator_api.php' on how to install it.")
}

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: daxgen.R DAXFILE\n")
}
daxfile <- args[1]

# Create an abstract DAG
workflow <- ADAG("advanced_bioconductor_workflow")

# Job 1- Preprocess
preprocess <- Job("e.preprocess")

# Input files
preprocess.xlsin <- File("GSE61533_HTSEQ_count_results.csv")
preprocess.xlsin <- AddPFN(preprocess.xlsin, 
                           PFN("file:///scratch/r-home/advanced_workflow/single-cell-rna-input-files/GSE61533_HTSEQ_count_results.csv", "local"))

workflow <- AddFile(workflow, preprocess.xlsin)

# Output files
preprocess.count <- File("count.csv")
preprocess.sce <- File("sce.RData")
preprocess.colNames <- File("column-names.txt")


# Specify input
preprocess <- Uses(preprocess, preprocess.xlsin, link = DAX3.Link$INPUT)

# Specify output
preprocess <- Uses(preprocess, preprocess.count, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)
preprocess <- Uses(preprocess, preprocess.sce, link = DAX3.Link$OUTPUT,
                   transfer = FALSE)
preprocess <- Uses(preprocess, preprocess.colNames, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)

# Add Job
workflow <- AddJob(workflow, preprocess)


## ---------------------------------------------------------------------------
# Job 2- Count loading

# countLoading <- Job("e.count")
# 
# # Input files
# countLoading.mRNA <- File("expression_mRNA_17-Aug-2014.txt")
# countLoading.spikes <- File("expression_spikes_17-Aug-2014.txt")
# countLoading.mito <- File("expression_mito_17-Aug-2014.txt")
# 
# countLoading.mRNA <- AddPFN(countLoading.mRNA, 
#                             PFN("file:///scratch/r-home/advanced_workflow/expression_mRNA_17-Aug-2014.txt", "local"))
# countLoading.spikes <- AddPFN(countLoading.spikes, 
#                               PFN("file:///scratch/r-home/advanced_workflow/expression_spikes_17-Aug-2014.txt", "local"))
# countLoading.mito <- AddPFN(countLoading.mito,
#                             PFN("file:///scratch/r-home/advanced_workflow/expression_mito_17-Aug-2014.txt", "local"))
# 
# workflow <- AddFile(workflow, countLoading.mRNA)
# workflow <- AddFile(workflow, countLoading.spikes)
# workflow <- AddFile(workflow, countLoading.mito)
# 
# # Specify inputs
# countLoading <- Uses(countLoading, countLoading.mRNA, link = DAX3.Link$INPUT)
# countLoading <- Uses(countLoading, countLoading.spikes, link = DAX3.Link$INPUT)
# countLoading <- Uses(countLoading, countLoading.mito, link = DAX3.Link$INPUT)
# 
# # Output files
# countLoading.scedims <- File("count-loading-sce-dimensions.csv")
# countLoading.sce <- File("sce-countLoading.RData")
# 
# # Specify outputs
# countLoading <- Uses(countLoading, countLoading.scedims, link = DAX3.Link$OUTPUT,
#                      transfer = TRUE)
# countLoading <- Uses(countLoading, countLoading.sce, link = DAX3.Link$OUTPUT,
#                      transfer = FALSE)
# 
# # Specify input arguments
# countLoading <- AddArguments(countLoading, list(countLoading.scedims, countLoading.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, countLoading)

## ---------------------------------------------------------------------------
# Job 3.1- Quality Control-1 (First run)
qualityControl.part1.1 <- Job("e.qc1.1")

# Output files
qualityControl.part1.1.libSize <- File("library-sizes-1.jpg")
qualityControl.part1.1.eg <- File("expressed-genes-1.1.jpg")

# Specify inputs
qualityControl.part1.1 <- Uses(qualityControl.part1.1, preprocess.sce, link = DAX3.Link$INPUT)

# Specify input arguments
qualityControl.part1.1 <- AddArguments(qualityControl.part1.1, list(preprocess.sce, "1e6", 
                                                                qualityControl.part1.1.libSize, 
                                                                qualityControl.part1.1.eg))

# Specify outputs
qualityControl.part1.1 <- Uses(qualityControl.part1.1, qualityControl.part1.1.libSize, 
                             link = DAX3.Link$OUTPUT, transfer = TRUE)
qualityControl.part1.1 <- Uses(qualityControl.part1.1, qualityControl.part1.1.eg, 
                             link = DAX3.Link$OUTPUT, transfer = TRUE)

# Add Job
workflow <- AddJob(workflow, qualityControl.part1.1)

## ---------------------------------------------------------------------------
# Job 3.2- Quality Control (Second run)
# qualityControl.part1.2 <- Job("e.qc1.2")
# 
# # Output files
# qualityControl.part1.2.libSize <- File("library-sizes-2.jpg")
# qualityControl.part1.2.eg <- File("expressed-genes-2.jpg")
# 
# # Specify inputs
# qualityControl.part1.2 <- Uses(qualityControl.part1.2, countLoading.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# qualityControl.part1.2 <- AddArguments(qualityControl.part1.2, list(countLoading.sce, "1e3", 
#                                                                 qualityControl.part1.2.libSize, 
#                                                                 qualityControl.part1.2.eg))
# 
# # Specify outputs
# qualityControl.part1.2 <- Uses(qualityControl.part1.2, qualityControl.part1.2.libSize, 
#                              link = DAX3.Link$OUTPUT, transfer = TRUE)
# qualityControl.part1.2 <- Uses(qualityControl.part1.2, qualityControl.part1.2.eg, 
#                              link = DAX3.Link$OUTPUT, transfer = TRUE)
# 
# # Add Job
# workflow <- AddJob(workflow, qualityControl.part1.2)

## ---------------------------------------------------------------------------
# Job 4.1- Quality Control-2 (First run)

qualityControl.part2.1 <- Job("e.qc2.1")

# Output files
qualityControl.part2.1.mito <- File("mito-proportion-1.jpg")
qualityControl.part2.1.ercc <- File("ERCC-prorportion-1.jpg")

# Specify inputs
qualityControl.part2.1 <- Uses(qualityControl.part2.1, preprocess.sce, link = DAX3.Link$INPUT)

# Specify input arguments
qualityControl.part2.1 <- AddArguments(qualityControl.part2.1, list(preprocess.sce, 
                                                                    qualityControl.part2.1.mito, 
                                                                    qualityControl.part2.1.ercc))

# Specify outputs
qualityControl.part2.1 <- Uses(qualityControl.part2.1, qualityControl.part2.1.mito, 
                               link = DAX3.Link$OUTPUT, transfer = TRUE)
qualityControl.part2.1 <- Uses(qualityControl.part2.1, qualityControl.part2.1.ercc, 
                               link = DAX3.Link$OUTPUT, transfer = TRUE)

# Add Job
workflow <- AddJob(workflow, qualityControl.part2.1)

## ---------------------------------------------------------------------------
# Job 4.2- Quality Control-2 (Second run)

# qualityControl.part2.2 <- Job("e.qc2.2")
# 
# # Output files
# qualityControl.part2.2.mito <- File("mito-proportion-2.jpg")
# qualityControl.part2.2.ercc <- File("ERCC-prorportion-2.jpg")
# 
# # Specify inputs
# qualityControl.part2.2 <- Uses(qualityControl.part2.2, countLoading.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# qualityControl.part2.2 <- AddArguments(qualityControl.part2.2, list(countLoading.sce, 
#                                                                     qualityControl.part2.2.mito, 
#                                                                     qualityControl.part2.2.ercc))
# 
# # Specify outputs
# qualityControl.part2.2 <- Uses(qualityControl.part2.2, qualityControl.part2.2.mito, 
#                                link = DAX3.Link$OUTPUT, transfer = TRUE)
# qualityControl.part2.2 <- Uses(qualityControl.part2.2, qualityControl.part2.2.ercc, 
#                                link = DAX3.Link$OUTPUT, transfer = TRUE)
# 
# # Add Job
# workflow <- AddJob(workflow, qualityControl.part2.2)


## ---------------------------------------------------------------------------
# Job 5.1- Filters (First run)

filters1.1 <- Job("e.filters1.1")

# Output files
filters1.1.filteredcells <- File("filtered-cells-1.csv")
filters1.1.sce <- File("sce-filters-1.RData")
filters1.1.pca <- File("pca-filter.jpg")

# Specify inputs
filters1.1 <- Uses(filters1.1, preprocess.sce, link = DAX3.Link$INPUT)

# Specify input arguments
filters1.1 <- AddArguments(filters1.1, list(preprocess.sce, filters1.1.sce, 
                                            filters1.1.filteredcells, filters1.1.pca))

# Specify outputs
filters1.1 <- Uses(filters1.1, filters1.1.filteredcells, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)
filters1.1 <- Uses(filters1.1, filters1.1.sce, link = DAX3.Link$OUTPUT,
                   transfer = FALSE)
filters1.1 <- Uses(filters1.1, filters1.1.pca, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)

# Add Job
workflow <- AddJob(workflow, filters1.1)

## ---------------------------------------------------------------------------
# Job 5.2- Filters (Second run)

# filters1.2 <- Job("e.filters1.2")
# 
# # Output files
# filters1.2.filteredcells <- File("filtered-cells-2.csv")
# filters1.2.sce <- File("sce-filters-2.RData")
# 
# # Specify outputs
# filters1.2 <- Uses(filters1.2, filters1.2.filteredcells, link = DAX3.Link$OUTPUT,
#                    transfer = TRUE)
# filters1.2 <- Uses(filters1.2, filters1.2.sce, link = DAX3.Link$OUTPUT,
#                    transfer = FALSE)
# 
# # Specify inputs
# filters1.2 <- Uses(filters1.2, countLoading.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# filters1.2 <- AddArguments(filters1.2, list(countLoading.sce, filters1.2.sce, 
#                                             filters1.2.filteredcells))
# 
# 
# # Add Job
# workflow <- AddJob(workflow, filters1.2)

## ---------------------------------------------------------------------------
# Job 6.1- Low abundance filter (First run)

lowAbundance1.1 <- Job("e.lowabundance1.1")

# Output files
lowAbundance1.1.keep <- File("keep-1.txt")
lowAbundance1.1.sce <- File("lowAbundance-sce-1.RData")
lowAbundance1.1.altKeep <- File("alt.keep-1.txt")
lowAbundance1.1.freqvsmean <- File("freq-vs-mean-1.jpg")

# Specify outputs
lowAbundance1.1 <- Uses(lowAbundance1.1, lowAbundance1.1.keep, link = DAX3.Link$OUTPUT,
                        transfer = TRUE)
lowAbundance1.1 <- Uses(lowAbundance1.1, lowAbundance1.1.sce, link = DAX3.Link$OUTPUT,
                        transfer = FALSE)
lowAbundance1.1 <- Uses(lowAbundance1.1, lowAbundance1.1.altKeep, link = DAX3.Link$OUTPUT,
                        transfer = TRUE)
lowAbundance1.1 <- Uses(lowAbundance1.1, lowAbundance1.1.freqvsmean, link = DAX3.Link$OUTPUT,
                        transfer = TRUE)

# Specify input
lowAbundance1.1 <- Uses(lowAbundance1.1, filters1.1.sce, link = DAX3.Link$INPUT)

# Specify input arguments
lowAbundance1.1 <- AddArguments(lowAbundance1.1, list(filters1.1.sce, lowAbundance1.1.keep,
                                                      lowAbundance1.1.sce, lowAbundance1.1.altKeep,
                                                      lowAbundance1.1.freqvsmean))

# Add Job
workflow <- AddJob(workflow, lowAbundance1.1)

## ---------------------------------------------------------------------------
# Job 6.2- Low abundance filter (Second run)

# lowAbundance1.2 <- Job("e.lowabundance1.2")
# 
# # Output files
# lowAbundance1.2.keep <- File("keep-2.txt")
# lowAbundance1.2.sce <- File("lowAbundance-sce-2.RData")
# 
# 
# # Specify outputs
# lowAbundance1.2 <- Uses(lowAbundance1.1, lowAbundance1.2.keep, link = DAX3.Link$OUTPUT,
#                         transfer = TRUE)
# lowAbundance1.2 <- Uses(lowAbundance1.1, lowAbundance1.2.sce, link = DAX3.Link$OUTPUT,
#                         transfer = FALSE)
# 
# 
# # Specify input
# lowAbundance1.2 <- Uses(lowAbundance1.2, filters1.2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# lowAbundance1.2 <- AddArguments(lowAbundance1.2, list(filters1.2.sce, lowAbundance1.2.keep,
#                                                       lowAbundance1.2.sce))
# 
# # Add job
# workflow <- AddJob(workflow, lowAbundance1.2)

## ---------------------------------------------------------------------------
# Job 7- Normalize-1

normalize1 <- Job("e.normalize.1")

# Output files
normalize1.sizeFactors <- File("size-factors.txt")
normalize1.libSize <- File("libSize-vs-sizeFactors-1.jpg")
normalize1.sce <- File("normalize-1.sce.RData")
normalize1.fontsize <- File("normalize-1.fontsize.RData")

# Specify outputs
normalize1 <- Uses(normalize1, normalize1.sizeFactors, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)
normalize1 <- Uses(normalize1, normalize1.libSize, link = DAX3.Link$OUTPUT,
                   transfer = TRUE)
normalize1 <- Uses(normalize1, normalize1.sce, link = DAX3.Link$OUTPUT,
                   transfer = FALSE)
normalize1 <- Uses(normalize1, normalize1.fontsize, link = DAX3.Link$OUTPUT,
                   transfer = FALSE)

# Specify input
normalize1 <- Uses(normalize1, lowAbundance1.1.sce, link = DAX3.Link$INPUT)

# Specify input arguments
normalize1 <- AddArguments(normalize1, list(lowAbundance1.1.sce))

# Add Job
workflow <- AddJob(workflow, normalize1)

## ---------------------------------------------------------------------------
# Job 8- Normalize-2

# normalize2 <- Job("e.normalize.2")
# 
# # Output files
# normalize2.libSize <- File("libSize-vs-sizeFactors-2.jpg")
# normalize2.sce <- File("normalize2-sce.RData")
# 
# # Specify outputs
# normalize2 <- Uses(normalize2, normalize2.libSize, link = DAX3.Link$OUTPUT,
#                    transfer = TRUE)
# normalize2 <- Uses(normalize2, normalize2.sce, link = DAX3.Link$OUTPUT,
#                    transfer = FALSE)
# 
# # Specify input
# normalize2 <- Uses(normalize2, lowAbundance1.2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# normalize2 <- AddArguments(normalize2, list(lowAbundance1.2.sce, normalize2.libSize,
#                                             normalize2.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, normalize2)

## ---------------------------------------------------------------------------
# Job 9- Reduction

reduction <- Job("e.reduction")

# Output files
reduction.pca <- File("pca-normalized.jpg")
reduction.tSNE <- File("tSNE-plots.jpg")

# Specify outputs
reduction <- Uses(reduction, reduction.pca, link = DAX3.Link$OUTPUT, transfer = TRUE)
reduction <- Uses(reduction, reduction.tSNE, link = DAX3.Link$OUTPUT, transfer = TRUE)

# Specify input 
reduction <- Uses(reduction, normalize1.sce, link = DAX3.Link$INPUT)
reduction <- Uses(reduction, normalize1.fontsize, link = DAX3.Link$INPUT)

# Add Job
workflow <- AddJob(workflow, reduction)

## ---------------------------------------------------------------------------
# Job 10.1 Cell Cycle
cellCycle1.1 <- Job("e.cellCyc1.1")

# Output file
cellCycle1.1.phase <- File("phase-assignment-1.jpg")
cellCycle1.1.assignments <- File("cell-cycle-1.assignments.RData")

# Specify output
cellCycle1.1 <- Uses(cellCycle1.1, cellCycle1.1.phase, link = DAX3.Link$OUTPUT,
                     transfer = TRUE)
cellCycle1.1 <- Uses(cellCycle1.1, cellCycle1.1.assignments, link = DAX3.Link$OUTPUT,
                     transfer = FALSE)

# Specify input
cellCycle1.1 <- Uses(cellCycle1.1, normalize1.sce, link = DAX3.Link$INPUT)

# Specify input arguments
cellCycle1.1 <- AddArguments(cellCycle1.1, list(normalize1.sce, cellCycle1.1.phase, 
                                                cellCycle1.1.assignments))

# Add Job
workflow <- AddJob(workflow, cellCycle1.1)

## ---------------------------------------------------------------------------
# Job 10.2 Cell Cycle
# cellCycle1.2 <- Job("e.cellCyc1.2")
# 
# # Output file
# cellCycle1.2.phase <- File("phase-assignment-2.jpg")
# cellCycle1.2.assignments <- File("cell-cycle-2.assignments.RData")
# # Specify output
# cellCycle1.2 <- Uses(cellCycle1.2, cellCycle1.2.phase, link = DAX3.Link$OUTPUT,
#                      transfer = TRUE)
#cellCycle1.2 <- Uses(cellCycle1.2, cellCycle1.2.assignments, link = DAX3.Link$OUTPUT,
#                     transfer = FALSE)
# 
# # Specify input
# cellCycle1.2 <- Uses(cellCycle1.2, normalize2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
#  cellCycle1.2 <- AddArguments(cellCycle1.2, list(normalize2.sce, cellCycle1.2.phase,
#                                                  cellCycle1.2.assignments))

# # Add Job
# workflow <- AddJob(workflow, cellCycle1.2)

## ---------------------------------------------------------------------------
# Job 11 HVGs setup-1

hvgs.setup.1 <- Job("e.hvgs1")

# Output files
hvgs.setup.1.sce <- File("hvgs-setup-1.sce.RData")
var.fit1 <- File("var.fit.1.RData")
var.out1 <- File("var.out.1.RData")
var.out2.1 <- File("var.out2.1.RData")


# Specify outputs
hvgs.setup.1 <- Uses(hvgs.setup.1, hvgs.setup.1.sce, link = DAX3.Link$OUTPUT,
                     transfer = FALSE)
hvgs.setup.1 <- Uses(hvgs.setup.1, var.fit1, link = DAX3.Link$OUTPUT,
                     transfer = FALSE)
hvgs.setup.1 <- Uses(hvgs.setup.1, var.out1, link = DAX3.Link$OUTPUT,
                     transfer = FALSE)
hvgs.setup.1 <- Uses(hvgs.setup.1, var.out2.1, link = DAX3.Link$OUTPUT,
                     transfer = FALSE)

# Specify input
hvgs.setup.1 <- Uses(hvgs.setup.1, normalize1.sce , link = DAX3.Link$INPUT)
hvgs.setup.1 <- Uses(hvgs.setup.1, cellCycle1.1.assignments, link = DAX3.Link$INPUT)

# Specify input arguments
hvgs.setup.1 <- AddArguments(hvgs.setup.1, list(cellCycle1.1.assignments))

# Add Job
workflow <- AddJob(workflow, hvgs.setup.1)

## ---------------------------------------------------------------------------
# Job 12 HVGs setup-2

# hvgs.setup.2 <- Job("e.hvgs2")
# 
# # Output files
# hvgs.setup.2.sce <- File("hvgs-setup-2.sce.RData")
# pca.tissue <- File("pca-tissue-sex.jpg")
# tsne.tissue <- File("tSNE-tissue-sex.jpg")
# c1.chip <- File("C1-chip-colored.jp")
# var.fit2 <- File("var.fit-2.RData")
# var.out2.2 <- File("var.out-2.RData")
# 
# # Specify outputs
# hvgs.setup.2 <- Uses(hvgs.setup.2, hvgs.setup.2.sce, link = DAX3.Link$OUTPUT,
#                      transfer = FALSE)
# hvgs.setup.2 <- Uses(hvgs.setup.2, pca.tissue, link = DAX3.Link$OUTPUT,
#                      transfer = TRUE)
# hvgs.setup.2 <- Uses(hvgs.setup.2, tsne.tissue, link = DAX3.Link$OUTPUT,
#                      transfer = TRUE)
# hvgs.setup.2 <- Uses(hvgs.setup.2, c1.chip, link = DAX3.Link$OUTPUT,
#                      transfer = TRUE)
# hvgs.setup.2 <- Uses(hvgs.setup.2, var.fit2, link = DAX3.Link$OUTPUT,
#                      transfer = FALSE)
# hvgs.setup.2 <- Uses(hvgs.setup.2, var.out2.2, link = DAX3.Link$OUTPUT,
#                      transfer = FALSE)
# 
# # Specify input
# hvgs.setup.2 <- Uses(hvgs.setup.2, normalize2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# hvgs.setup.2 <- AddArguments(hvgs.setup.2, list(normalize2.sce, hvgs.setup.2.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, hvgs.setup.2)

## ---------------------------------------------------------------------------
# Job 13 Log Expression-1

log.expression.1 <- Job("e.logExp1")

# Output files
log.expression.1.jpg <- File("log-expression-1.jpg")
log.expression.1.topData <- File("top-hvgs.RData")
log.expression.1.top <- File("top-hvgs.csv")
log.expression.1.tsv <- File("hsc-hvg.tsv")
log.expression.1.normalized <- File("normalized-log-expression.jpg")

# Specify outputs
log.expression.1 <- Uses(log.expression.1, log.expression.1.jpg, link = DAX3.Link$OUTPUT,
                         transfer = TRUE)
log.expression.1 <- Uses(log.expression.1, log.expression.1.topData, link = DAX3.Link$OUTPUT,
                         transfer = FALSE)
log.expression.1 <- Uses(log.expression.1, log.expression.1.top, link = DAX3.Link$OUTPUT,
                         transfer = TRUE)
log.expression.1 <- Uses(log.expression.1, log.expression.1.tsv, link = DAX3.Link$OUTPUT,
                         transfer = TRUE)
log.expression.1 <- Uses(log.expression.1, log.expression.1.normalized, link = DAX3.Link$OUTPUT,
                         transfer = TRUE)

# Specify inputs
log.expression.1 <- Uses(log.expression.1, var.fit1, link = DAX3.Link$INPUT)
log.expression.1 <- Uses(log.expression.1, var.out1, link = DAX3.Link$INPUT)
log.expression.1 <- Uses(log.expression.1, var.out2.1, link = DAX3.Link$INPUT)
log.expression.1 <- Uses(log.expression.1, hvgs.setup.1.sce, link = DAX3.Link$INPUT)


# Add Job
workflow <- AddJob(workflow, log.expression.1)


## ---------------------------------------------------------------------------
# Job 14 Log Expression-2

# log.expression.2 <- Job("e.logExp2")
# 
# # Output files
# log.expression.2.jpg <- File("log-expression-2.jpg")
# log.expression.2.top <- File("top-hvgs-2.csv")
# log.expression.2.hvgtsv <- File("brain-hvg.tsv")
# log.expression.2.normalized <- File("normalized-log-expression-2.jpg")
# log.expression.2.cortsv <- File("brain-cor.tsv")
# log.expression.2.genePairs <- File("gene-pairs-2.csv")
# log.expression.2.corGenes <- File("correlated-genes-2.csv")
# log.expression.2.varCor <- File("var.cor.RData")
# log.expression.2.sigCor <- File("sig.cor.RData")
# 
# # Specify outputs
# log.expression.2 <- Uses(log.expression.2, log.expression.2.jpg, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.top, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.hvgtsv, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.normalized, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.cortsv, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.genePairs, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.corGenes, link = DAX3.Link$OUTPUT,
#                          transfer = TRUE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.varCor, link = DAX3.Link$OUTPUT,
#                          transfer = FALSE)
# log.expression.2 <- Uses(log.expression.2, log.expression.2.sigCor, link = DAX3.Link$OUTPUT,
#                          transfer = FALSE)
# 
# # Specify inputs
# log.expression.2 <- Uses(log.expression.2, var.fit2, link = DAX3.Link$INPUT)
# log.expression.2 <- Uses(log.expression.2, var.out2.2, link = DAX3.Link$INPUT)
# log.expression.2 <- Uses(log.expression.2, hvgs.setup.2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# log.expression.2 <- AddArguments(log.expression.2, list(hvgs.setup.2.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, log.expression.2)

## ---------------------------------------------------------------------------
# Job 15 Correlate Pairs-1

correlate.pairs.1 <- Job("e.corPairs1")

# Output files
correlate.pairs.1.hsctsv <- File("hsc-cor.tsv")
correlate.pairs.1.genePairs <- File("gene-pairs-1.csv")
correlate.pairs.1.corPairs <- File("correlated-genes-1.csv")
correlate.pairs.1.pairCor <- File("pair-wise-correlations.csv")
correlate.pairs.1.heatmap <- File("heatmap.jpg")
correlate.pairs.1.tsne <- File("pca-tSNE-plots.jpg")
correlate.pairs.1.hsc <- File("hsc-data.rds")

# Specify outputs
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.hsctsv, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.genePairs, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.corPairs, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.pairCor, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.heatmap, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.tsne, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
correlate.pairs.1 <- Uses(correlate.pairs.1, correlate.pairs.1.hsc, link = DAX3.Link$OUTPUT,
                          transfer = TRUE)
# Specify inputs
correlate.pairs.1 <- Uses(correlate.pairs.1, hvgs.setup.1.sce, link = DAX3.Link$INPUT)
correlate.pairs.1 <- Uses(correlate.pairs.1, log.expression.1.topData, link = DAX3.Link$INPUT)

# Add Job
workflow <- AddJob(workflow, correlate.pairs.1)

## ---------------------------------------------------------------------------
# Job 16 Correlate Pairs-2

# correlate.pairs.2 <- Job("e.corPairs2")
# 
# # Output files
# correlate.pairs.2.brain <- File("brain.heat.pdf")
# correlate.pairs.2.dispersion <- File("dispersion-summary.csv")
# correlate.pairs.2.cluster <- File("cluster.Rdata")
# correlate.pairs.2.design <- File("design.Rdata")
# correlate.pairs.2.y <- File("y.Rdata")
# correlate.pairs.2.fit <- File("fit.Rdata")
# correlate.pairs.2.tree <- File("my.tree.Rdata")
# correlate.pairs.2.my.clusters <- File("my.clusters.Rdata")
# correlate.pairs.2.clustCol <- File("clust.col.Rdata")
# 
# # Specify outputs
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.brain, link = DAX3.Link$OUTPUT,
#                           transfer = TRUE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.dispersion, link = DAX3.Link$OUTPUT,
#                           transfer = TRUE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.cluster, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.design, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.y, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.fit, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.tree, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.my.clusters, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# correlate.pairs.2 <- Uses(correlate.pairs.2, correlate.pairs.2.clustCol, link = DAX3.Link$OUTPUT,
#                           transfer = FALSE)
# 
# # Specify inputs
# correlate.pairs.2 <- Uses(correlate.pairs.2, log.expression.2.varCor, link = DAX3.Link$INPUT)
# correlate.pairs.2 <- Uses(correlate.pairs.2, log.expression.2.sigCor, link = DAX3.Link$INPUT)
# correlate.pairs.2 <- Uses(correlate.pairs.2, hvgs.setup.2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# correlate.pairs.2 <- AddArguments(correlate.pairs.2, list(hvgs.setup.2.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, correlate.pairs.2)

## ---------------------------------------------------------------------------
# Job 17 Cell cycle blocking

cell.cycle.blocking <- Job("e.cellBlock")

# Output files
cell.block.normalization <- File("normalization-vs-deconvolution.jpg")
cell.block.phase <- File("phase-scores.jpg")
cell.block.removal <- File("cell-cycle-removal.jpg")
cell.block.diffusion <- File("diffusion-map.jpg")
cell.block.sce <- File("cell.cycle.blocking.sce.RData")

# Specify outputs
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.block.normalization, link = DAX3.Link$OUTPUT,
                            transfer = TRUE)
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.block.phase, link = DAX3.Link$OUTPUT,
                            transfer = TRUE)
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.block.removal, link = DAX3.Link$OUTPUT,
                            transfer = TRUE)
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.block.diffusion, link = DAX3.Link$OUTPUT,
                            transfer = TRUE)
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.block.sce, link = DAX3.Link$OUTPUT,
                            transfer = FALSE)

# Input file
cell.cycle.csv <- File("nbt.3102-S7.csv")
cell.cycle.csv <- AddPFN(cell.cycle.csv, 
                         PFN("file:///scratch/r-home/advanced_workflow/single-cell-rna-input-files/nbt.3102-S7.csv", "local"))
workflow <- AddFile(workflow, cell.cycle.csv)

# Specify input
cell.cycle.blocking <- Uses(cell.cycle.blocking, correlate.pairs.1.hsc, link = DAX3.Link$INPUT)
cell.cycle.blocking <- Uses(cell.cycle.blocking, cell.cycle.csv, link = DAX3.Link$INPUT)

# Specify input arguments
cell.cycle.blocking <- AddArguments(cell.cycle.blocking, list(cell.block.sce))

# Add Job
workflow <- AddJob(workflow, cell.cycle.blocking)

## ---------------------------------------------------------------------------
# Job 18 Extract Annotation

extract <- Job("e.extract")

# Output files
extract.ENSEMBL <- File("ENSEMBL-identifiers.csv")
extract.scedims <- File("sce-dimensions.csv")

# Specify outputs
extract <- Uses(extract, extract.ENSEMBL, link = DAX3.Link$OUTPUT,
                transfer = TRUE)
extract <- Uses(extract, extract.scedims, link = DAX3.Link$OUTPUT,
                transfer = TRUE)

# Input file
extract.count <- File("counttable_es.csv")
extract.count <- AddPFN(extract.count, 
                        PFN("file:///scratch/r-home/advanced_workflow/counttable_es.csv", "local"))
workflow <- AddFile(workflow, extract.count)

# Specify inputs
extract <- Uses(extract, extract.count, link = DAX3.Link$INPUT)
extract <- Uses(extract, cell.block.sce, link = DAX3.Link$INPUT)

# Specify input arguments
extract <- AddArguments(extract, list(cell.block.sce))

# Add Job
workflow <- AddJob(workflow, extract)

## ---------------------------------------------------------------------------
# Job 19 LRT

# lrt <- Job("e.lrt")
# 
# # Output files
# lrt.markers <- File("markers.csv")
# lrt.heatmap <- File("cluster2-heatmap.jpg")
# lrt.brain <- File("brain-marker-2.tsv")
# 
# # Specify outputs
# lrt <- Uses(lrt, lrt.markers, link = DAX3.Link$OUTPUT, transfer = TRUE)
# lrt <- Uses(lrt, lrt.heatmap, link = DAX3.Link$OUTPUT, transfer = TRUE)
# lrt <- Uses(lrt, lrt.brain, link = DAX3.Link$OUTPUT, transfer = TRUE)
# 
# # Specify inputs
# lrt <- Uses(lrt, correlate.pairs.2.cluster, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.tree, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.design, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.y, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.my.clusters, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.fit, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, correlate.pairs.2.clustCol, link = DAX3.Link$INPUT)
# lrt <- Uses(lrt, hvgs.setup.2.sce, link = DAX3.Link$INPUT)
# 
# # Specify input arguments
# lrt <- AddArguments(lrt, list(hvgs.setup.2.sce))
# 
# # Add Job
# workflow <- AddJob(workflow, lrt)

## ---------------------------------------------------------------------------
# Job 20 Log fold

# log <- Job("e.log")
# 
# # Output files
# log.genes <- File("DE-genes.csv")
# 
# # Specify output
# log <- Uses(log, log.genes, link = DAX3.Link$OUTPUT, transfer = FALSE)
# 
# # Specify inputs
# log <- Uses(log, correlate.pairs.2.cluster, link = DAX3.Link$INPUT)
# log <- Uses(log, correlate.pairs.2.design, link = DAX3.Link$INPUT)
# log <- Uses(log, correlate.pairs.2.y, link = DAX3.Link$INPUT)
# log <- Uses(log, correlate.pairs.2.fit, link = DAX3.Link$INPUT)
# 
# # Add Job
# workflow <- AddJob(workflow, log)


# Dependencies
workflow <- Depends(workflow, parent = preprocess, child = qualityControl.part1.1)
#workflow <- Depends(workflow, parent = countLoading, child = qualityControl.part1.2)
workflow <- Depends(workflow, parent = preprocess, child = qualityControl.part2.1)
#workflow <- Depends(workflow, parent = countLoading, child = qualityControl.part2.2)
workflow <- Depends(workflow, parent = preprocess, child = filters1.1)
#workflow <- Depends(workflow, parent = countLoading, child = filters1.2)
workflow <- Depends(workflow, parent = filters1.1, child = lowAbundance1.1)
#workflow <- Depends(workflow, parent = filters1.2, child = lowAbundance1.2)
workflow <- Depends(workflow, parent = lowAbundance1.1, child = normalize1)
#workflow <- Depends(workflow, parent = lowAbundance1.2, child = normalize2)
workflow <- Depends(workflow, parent = normalize1, child = reduction)
workflow <- Depends(workflow, parent = normalize1, child = cellCycle1.1)
workflow <- Depends(workflow, parent = cellCycle1.1, child = hvgs.setup.1)
workflow <- Depends(workflow, parent = normalize1, child = hvgs.setup.1)
#workflow <- Depends(workflow, parent = normalize2, child = cellCycle1.2)
#workflow <- Depends(workflow, parent = normalize2, child = hvgs.setup.2)
#workflow <- Depends(workflow, parent = cellCycle1.2, child = hvgs.setup.2)
workflow <- Depends(workflow, parent = hvgs.setup.1, child = log.expression.1)
#workflow <- Depends(workflow, parent = hvgs.setup.2, child = log.expression.2)
workflow <- Depends(workflow, parent = log.expression.1, child = correlate.pairs.1)
#workflow <- Depends(workflow, parent = log.expression.2, child = correlate.pairs.2)
workflow <- Depends(workflow, parent = correlate.pairs.1, child = cell.cycle.blocking)
workflow <- Depends(workflow, parent = cell.cycle.blocking, child = extract)
#workflow <- Depends(workflow, parent = correlate.pairs.2, child = lrt)
#workflow <- Depends(workflow, parent = correlate.pairs.2, child = log)


# Write the DAX to file
WriteXML(workflow, daxfile)