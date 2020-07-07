##############################################################
#                                 #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/fastMNN

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load fastMNN corrected dataset
### 2) Score with Developmental and Injury Response gene signatures
### 3) Save Data
##############################################################

library(Seurat) #v3.1.5
library(AUCell) #v1.8.0
library(BBmisc)
library(SeuratWrappers)


##############################################################
### 1) Load data
##############################################################
### fastMNN corrected data
### this is kind of messed up because its only the top 2000 most variable genes in the reconctructed matrix
#test <- readRDS("Global_SU2C_GSCs_Seurat_fastMNN.rds")
#exprMatrix <- data.matrix(test@tools$RunFastMNN@assays@data@listData$reconstructed)
#dim(exprMatrix)

### load global GSC dataset
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
BTSC <- UpdateSeuratObject(BTSC)

### load and subset gene sigs
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AUCell_Signatures_Hypoxia.Rdata")
sigs <- sigs[c("Developmental_GSC", "InjuryResponse_GSC")]



##############################################################
### 2) Run fastMNN on all genes
##############################################################
genes <- rownames(BTSC@assays$RNA@data)
BTSC_fMNN <- RunFastMNN(object.list = SplitObject(BTSC, split.by = "SampleID"), features = genes)
#BTSC_fMNN <- RunUMAP(BTSC_fMNN, reduction = "mnn", dims = 1:30)
#BTSC_fMNN <- FindNeighbors(BTSC_fMNN, reduction = "mnn", dims = 1:30)
#BTSC_fMNN <- FindClusters(BTSC_fMNN, resolution = 2)
#saveRDS(BTSC_fMNN, file = "Global_SU2C_GSCs_Seurat_fastMNN_ALLGENES.rds")
rm(BTSC_fMNN)

### isolate corrected expression Matrix
exprMatrix <- data.matrix(BTSC_fMNN@tools$RunFastMNN@assays@data@listData$reconstructed)
dim(exprMatrix)
saveRDS(exprMatrix, file = "Global_GSCs_fastMNN_correctedExpMat_ALLGENES.rds")



##############################################################
### 3) Run AUCell with Dev and IR signatures
##############################################################
#exprMatrix <- readRDS("Global_GSCs_fastMNN_correctedExpMat_ALLGENES.rds")
#load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AUCell_Signatures_Hypoxia.Rdata")
#sigs <- sigs[c("Developmental_GSC", "InjuryResponse_GSC")]

cells_rankings <- AUCell_buildRankings(exprMatrix,
                                       #nCores=1,
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

x <- strsplit(rownames(AUC), "_")
AUC$SampleID <- sapply( x, "[", 1)
AUC$SampleID <- paste0(AUC$SampleID, "_L")

##############################################################
### 3) Corrlelate scores to original
##############################################################




##############################################################
### 4) Save data
##############################################################

saveRDS(AUC, file = "fastMNN_DevIR_AUCell_GSCs.rds")
### save AUCell scores + metadata
#BTSC_AUCell <- cbind(BTSC@meta.data, AUC)
#saveRDS(BTSC_AUCell, file = "fastMNN_DevIR_AUCell_meta_GSCs.rds")
