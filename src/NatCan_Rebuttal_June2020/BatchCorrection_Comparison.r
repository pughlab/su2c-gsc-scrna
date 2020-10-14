##############################################################
#             Compre Batch Correction Tools on GSCs          #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: https://github.com/satijalab/seurat-wrappers
### Working dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/Comparison


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Caclulate # clusters / method
### 2) Calculate # samples / cluster
### 3) Calculate mean silhouette width / cluster
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
library(ggplot2)
library(cluster)



##############################################################
# 2) Load Data
##############################################################

meta <- readRDS("Global_GSC_BatchCorrection_metadata.rds")

load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
original <- UpdateSeuratObject(BTSC)
saveRDS(original, file = "../original_clustering/Global_SU2C_GSCs_Seurat_Original.rds")

conos <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/CONOS/Global_SU2C_GSCs_Seurat.CONOS.rds")

liger <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/Liger/Global_SU2C_GSCs_Seurat_LIGER.rds")

fastMNN <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/fastMNN/Global_SU2C_GSCs_Seurat_fastMNN.rds")



##############################################################
# 3) Caclulate shiluoette width
##############################################################

### 3.1) Original
pca.dist <- as.matrix(original@reductions$pca@cell.embeddings[ ,1:30])
pca.dist <- dist(pca.dist)
cl <- original@meta.data[ ,"res.2"]
cl <- as.integer(gsub("C", "", cl))
table(cl)
original_sil <- cluster::silhouette(x=cl,
                          dist = pca.dist,
                          do.clus.stat = TRUE
                        )
saveRDS(sil, file = "Original_silWidth.rds")


### 3.2) Conos (doesnt have PCA)
pca.dist <- as.matrix(conos@reductions$pca@cell.embeddings[ ,1:30])
pca.dist <- dist(pca.dist)
cl <- conos@meta.data[ ,"Original_clusters"]
table(cl)
#calc sil
sil <- cluster::silhouette(as.integer(cl), pca.dist)
saveRDS(sil, file = "Conos_silWidth.rds")



### 3.3) Liger
pca.dist <- as.matrix(liger@reductions$pca@cell.embeddings[ ,1:30])
pca.dist <- dist(pca.dist)
cl <- as.integer(liger@meta.data[ ,"seurat_clusters"])
table(cl)
#calc sil
sil <- cluster::silhouette(cl, pca.dist)
saveRDS(sil, file = "Liger_silWidth.rds")
