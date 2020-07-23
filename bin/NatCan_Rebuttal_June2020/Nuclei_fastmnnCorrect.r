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

####################
# Save data
####################

### isolate corrected expression Matrix
exprMatrix <- data.matrix(dat_fMNN@tools$RunFastMNN@assays@data@listData$reconstructed)
dim(exprMatrix)
saveRDS(exprMatrix, file = "GSCs_Tumour_LiveNuclei_fastMNN_correctedExpMat_ALLGENES.rds")

## save fastMNN corrected seurat object
saveRDS(dat_fMNN, file = "GSCs_Tumour_LiveNuclei_fastMNN_Seurat.rds")

### isolate and save corrected PCA
mnn <- data.matrix(dat_fMNN@reductions$mnn@cell.embeddings)
saveRDS(mnn, file = "GSCs_Tumour_LiveNuclei_fastMNN_mnn_reducedDims.rds")


##############################################################
# 4) Score corrected gene matrix with Dev and IR signatures
##############################################################
### load corrected fastMNN expression matrix
exprMatrix <- readRDS("GSCs_Tumour_LiveNuclei_fastMNN_correctedExpMat_ALLGENES.rds")

### load gene signautres
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AUCell_Signatures_Hypoxia.Rdata")
sigs <- sigs[c("Developmental_GSC", "InjuryResponse_GSC")]

### run AUCell
cells_rankings <- AUCell_buildRankings(exprMatrix,
                                       nCores=30,
                                       plotStats=FALSE)

cells_AUC <- AUCell_calcAUC(sigs,
                            cells_rankings)

cells_assignment <- AUCell_exploreThresholds(cells_AUC,
                                             plotHist=FALSE,
                                             assign=TRUE
                                            )

AUC <- t(as.data.frame(cells_AUC@assays@data$AUC))
colnames(AUC) <- paste0(colnames(AUC), "_AUC")
AUC <- data.frame(AUC)
print(AUC[1:2, 1:2])

saveRDS(AUC, file = "AUCell_CorrectedMat_fastMNN_Nuclei.rds")
