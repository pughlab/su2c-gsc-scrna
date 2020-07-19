##############################################################
#             Integrate nuclei data with fastMNN()           #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html#3_mutual_nearest_neighbors
### Working dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Nuclei_integration/fastmnn

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load GSC + tumour + nuclei dataset
### 2) Run fastMNN on ALL genes
### 3) Recore corrected matrix with Dev and IR
### 4) Save results
##############################################################


##############################################################
# 1) Install + Load packages
##############################################################
library(Seurat) #v3.1.5
library(AUCell) #v1.8.0
library(BBmisc)
library(SeuratWrappers)


##############################################################
# 2) Load and re-normalize data (log2)
##############################################################

dat <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/SU2C_LiveNuc_GSC_GBM_Malignant_Seurat.rds")

### dd seq tech to meta.data
technology <- dat@meta.data$orig.ident
technology <- gsub("BTSC", "cells", technology)
technology <- gsub("GBM", "cells", technology)
technology <- ifelse(technology == "cells", print("cells"), print("nuclei"))
dat@meta.data$technology <- technology
dat <- UpdateSeuratObject(dat)

##############################################################
# 3) Run fastMNN with seurat wrapper
##############################################################
ribo.genes <- c(rownames(dat@assays$RNA@data)[grep("^RP[1:9]", rownames(dat@assays$RNA@data))],
              rownames(dat@assays$RNA@data)[grep("^RP[L,S]", rownames(dat@assays$RNA@data))]
             ) #667 genes
genes <- rownames(dat@assays$RNA@data)[!rownames(dat@assays$RNA@data) %in% ribo.genes] #19486 genes

dat_fMNN <- RunFastMNN(object.list = SplitObject(dat, split.by = "technology"), features = genes)

### isolate corrected expression Matrix
exprMatrix <- data.matrix(dat_fMNN@tools$RunFastMNN@assays@data@listData$reconstructed)
dim(exprMatrix)
saveRDS(exprMatrix, file = "GSCs_Tumour_LiveNuclei_fastMNN_correctedExpMat_ALLGENES.rds")

## save fastMNN corrected seurat object
saveRDS(dat_fastMNN, file = "GSCs_Tumour_LiveNuclei_fastMNN_Seurat.rds")

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
