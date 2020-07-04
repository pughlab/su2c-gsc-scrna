##############################################################
#               Scale Data within each sample first          #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: https://github.com/JinmiaoChenLab/Rphenograph
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load GSC + tumour dataset
### 2) Split by sample
### 3) Scale each sample
### 4) Re-run PCA + clustering etc
### 5) Save + Plot
##############################################################

library(Seurat)

##############################################################
### 1) Load data
##############################################################

dat <- readRDS("../data/SU2C_Live_GSC_GBM_Malignant_Seurat.rds")
pc.use <-
genes.use <-
dat <- UpdateSeuratObject(dat) #79862 cells


##############################################################
### 2) Split and Scale by Sample
##############################################################
### get a list of scaled data frames for each sample
split <- ScaleData(dat,
                  split.by = "SampleID",
                  do.center = FALSE,
                  vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                  verbose = T
                    )

##############################################################
### 3) Re-run PCA in same way as manuscript (Figure 5A)
##############################################################

ribo.genes <- c(rownames(split@assays$RNA@data)[grep("^RP[1:9]",
                rownames(split@assays$RNA@data))],
                rownames(split@assays$RNA@data)[grep("^RP[L,S]", rownames(split@assays$RNA@data))]
               ) #2006

split <- RunPCA(split,
               pc.genes = rownames(split@assays$RNA@data)[!rownames(split@assays$RNA@data) %in% ribo.genes],
               pcs.compute = 20,
               genes.print = 10
              )

##############################################################
### 4) Save data
##############################################################
meta <- split@meta.data
pca <- split@reductions$pca
saveRDS(meta, file = "ScaledBySample_GCC_GBM_meta.rds")
saveRDS(pca, file = "ScaledBySample_GCC_GBM_pca.rds")
saveRDS(split, file = "ScaledBySample_GCC_GBM_Seurat.rds")

meta_orig <- dat@meta.data
pca_orig <- dat@reductions$pca
saveRDS(meta_orig, file = "Original_GCC_GBM_meta.rds")
saveRDS(pca_orig, file = "Original_GCC_GBM_pca.rds")
