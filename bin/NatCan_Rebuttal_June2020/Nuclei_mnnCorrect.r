##############################################################
#             Integrate nuclei data with mnnCorrect()        #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html#3_mutual_nearest_neighbors
### Working dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Nuclei_integration/mnn

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load GSC + tumour + nuclei dataset
### 2) Run mnnCorrect from Batchelor package
### 3) Recore corrected matrix with Dev and IR
### 4) Save results
##############################################################


##############################################################
# 1) Install + Load packages
##############################################################
#library(devtools)
#devtools::install_github("RcppCore/Rcpp")
#devtools::install_github("satijalab/seurat-wrappers")
#devtools::install_github("hms-dbmi/conos")
#devtools::install_github('MacoskoLab/liger')
library(Seurat) #v3.1.5
library(batchelor) #v1.2.4 (used in fastMNN seurat wrapper)
library(scater)


##############################################################
# 2) Load and re-normalize data (log2)
##############################################################

dat <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/SU2C_LiveNuc_GSC_GBM_Malignant_Seurat.rds")

### dd seq tech to meta.data
technology <- dat@meta.data$orig.ident
technology <- gsub("BTSC", "cells", technology)
technology <- gsub("GBM", "cells", technology)
technology <- ifelse(technology == "cells" ,print("cells"), print("nuclei"))
dat@meta.data$technology <- technology

### isolate cell counts
cell.counts <- dat@raw.data[, rownames(dat@meta.data[dat@meta.data$technology == "cells", ])]
dim(cell.counts)

### isoalte nuclei counts
nuclei.counts <- dat@raw.data[, rownames(dat@meta.data[dat@meta.data$technology == "nuclei", ])]
dim(nuclei.counts)

### log normalize data with scater
### for now, do this without size factors
sce_nuclei <- SingleCellExperiment(list(counts=nuclei.counts))
sce_nuclei <- logNormCounts(sce_nuclei)
saveRDS(sce_nuclei, file = "snRNA_LogNormCounts.rds")

sce_cells <- SingleCellExperiment(list(counts=cell.counts))
sce_cells <- logNormCounts(sce_cells)
saveRDS(sce_cells, file = "scRNA_LogNormCounts.rds")


##############################################################
# 3) Run Mutual Nearest Neighbours
##############################################################
### cos.norm.out = When ‘cos.norm.out=TRUE’, cosine normalization is performed
### on the matrix of values used to calculate correction vectors
### (and on which those vectors are applied). This can be turned
### off to obtain corrected values on the log-scale, similar to
### the input data.
###
### Interesting thread: http://supportupgrade.bioconductor.org/p/120471/#120512
###
##############################################################
### EXAMPLE EXECUTION ON H4H
### #!/bin/bash
### #SBATCH -t 24:00:00
### #SBATCH --mem=150G
### #SBATCH -p veryhimem
### #SBATCH -c 30
### #SBATCH -N 1
### #SBATCH --account=pughlab
###
### module load R/3.6.1
###
### Rscript mnn_nuclei.r
##############################################################

#library(Seurat) #v3.1.5
library(batchelor) #v1.2.4 (used in fastMNN seurat wrapper)
library(scater) #v1.14.6
library(parallel)
library(BiocParallel)

print("Configure Parallelization with ")
print(Sys.time())
options(MulticoreParam=quote(MulticoreParam(workers=multicoreWorkers())))

print("Loading Data")
print(Sys.time())
sce_cells <- readRDS("scRNA_LogNormCounts.rds")
sce_nuclei <- readRDS("snRNA_LogNormCounts.rds")

print("Run Mnn")
print(Sys.time())
### with more cells, a larger k should be used to achieve better merging
### in non-orthogonality. To achieve this, we can set prop.k, which
### allows k to adapt to batch size
### A numeric scalar in (0, 1) specifying the proportion of cells
### in each dataset to use for mutual nearest neighbor searching.
###  If set, ‘k’ for the search in each batch is redefined as
### ‘max(k, prop.k*N)’ where ‘N’ is the number of cells in that batch.
mnn_out <- mnnCorrect(sce_cells,
                      sce_nuclei,
                      prop.k = 0.8,
                      BPPARAM = MulticoreParam(),
                      correct.all = TRUE,
                      cos.norm.out = FALSE #output corrected log-scale values
                    )

print("Saving Data")
print(Sys.time())
saveRDS(mnn_out, file = "Nuclei_liveCell_mnn.rds")
