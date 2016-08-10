# advanced_bioconductor_workflow Workflow

Generating a Workflow
---------------------
Run the generate_dax.sh script.

Running a Workflow
-------------------
Run the plan_dax.sh script.

Libraries Needed:
------------------
library(scran)
library(edgeR)
library(scater)
library(Rtsne)
library(mvoutlier)
library(destiny)
library(gplots)
library(gdata)
library(R.utils)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

How to get all the input files:
--------------------------------
# Run the following R Script in a R environment.

all.urls <- c("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61533&format=file&file=GSE61533%5FHTSEQ%5Fcount%5Fresults%2Exls%2Egz", 
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mito_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt",
"http://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv", 
"http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx")
all.basenames <- basename(all.urls)
all.basenames[1] <- "GSE61533_HTSEQ_count_results.xls.gz"
all.modes <- rep("w", length(all.urls))
all.modes[c(1, 6)] <- "wb"
for (x in seq_along(all.urls)) { 
    download.file(all.urls[x], all.basenames[x], mode=all.modes[x])
}

The workflow scripts will not run with xlsx files as inputs so the scripts have been modified
to use csv files instead. Once this script is finished, the user needs to open all the xlsx
files in excel or a similar program and export the first sheet of each xlsx file as a csv file.

Workflow Description:
----------------------


