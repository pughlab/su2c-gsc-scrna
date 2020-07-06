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
######## USER DEFINED PARAMS AND FUNCTIONS
#global GSCs
dat.file <-
outFileName <-

#GSC+Tumour
dat.file <- "../../data/SU2C_Live_GSC_GBM_Malignant_Seurat.rds"
outFileName <- "GSCs_GBM"

##loads an RData file, and returns it with new name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}
##############################################################


##############################################################
### 1) Load data
##############################################################
###load data
if(grepl(".rds$", dat.file)){
  dat <- readRDS(dat.file)
} else if (grepl(".data$", dat.file)){
  dat <- loadRData(dat.file)
}
### extract processing information from original object
pca.genes <- dat@calc.params$RunPCA$pc.genes
pc.use <- 1:30
res.use <-  2
dat <- UpdateSeuratObject(dat) #79862 cells


##############################################################
### 2) Split and Scale by Sample
##############################################################
### get a list of scaled data frames for each sample
split <- ScaleData(dat,
                  split.by = "SampleID",
                  do.center = TRUE,
                  vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                  verbose = T
                    )


##############################################################
### 3) Re-run PCA +clustering in same way as manuscript
##############################################################
split <- RunPCA(split,
               pc.genes = pca.genes,
               pcs.compute = 50,
               genes.print = 0
              )
split <- FindNeighbors(split)
split <- FindClusters(split, resolution = res.use)
split <- RunUMAP(split, dims = pc.use)


##############################################################
### 4) Correlate PCA results
##############################################################



##############################################################
### 5) Save data
##############################################################

### combine meta.data

### combine pca.data


meta <- split@meta.data
pca <- split@reductions$pca
saveRDS(meta, file = "ScaledBySample_GCC_GBM_meta.rds")
saveRDS(pca, file = "ScaledBySample_GCC_GBM_pca.rds")


saveRDS(split, file = "ScaledBySample_GCC_GBM_Seurat.rds")

meta_orig <- dat@meta.data
pca_orig <- dat@reductions$pca
saveRDS(meta_orig, file = "Original_GCC_GBM_meta.rds")
saveRDS(pca_orig, file = "Original_GCC_GBM_pca.rds")


##############################################################
### 5) Repeat steps above for GSC only dataset
##############################################################

### load data

pc.use <- 1:30
dat <- UpdateSeuratObject(dat) #69k cells

### scale by sample
split <- ScaleData(dat,
                  split.by = "SampleID",
                  do.center = TRUE,
                  vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                  verbose = T
                  )
###  run PCA
ribo.genes <- c(rownames(split@assays$RNA@data)[grep("^RP[1:9]",
                rownames(split@assays$RNA@data))],
                rownames(split@assays$RNA@data)[grep("^RP[L,S]", rownames(split@assays$RNA@data))]
                                   ) #2006
split <- RunPCA(split,
                pc.genes = rownames(split@assays$RNA@data)[!rownames(split@assays$RNA@data) %in% ribo.genes],
                pcs.compute = 20,
                genes.print = 10
                )

### Find clusters (res 2)
split <- FindNeighbors(split)
split <- FindClusters(split, resolution = 2)

### Run UMAP
split <- RunUMAP(split, dims = pc.use)

### save data
meta <- cbind(split@meta.data, split@reductions$umap@cell.embeddings)
pca <- split@reductions$pca
saveRDS(meta, file = "ScaledBySample_GCCs_meta.rds")
saveRDS(pca, file = "ScaledBySample_GCCs_pca.rds")
saveRDS(split, file = "ScaledBySample_GCCs_Seurat.rds")

meta_orig <- cbind(dat@meta.data, dat@reductions$umap@cell.embeddings)
pca_orig <- dat@reductions$pca
saveRDS(meta_orig, file = "Original_GCC_meta.rds")
saveRDS(pca_orig, file = "Original_GCC_pca.rds")


##############################################################
### 6) (GSCs+Tumour) Compare original and scaled
##############################################################

### 6.1) Correlate cell embeddings between processing methods

orig_pca <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling/GSC_GBM/Original_GCC_GBM_pca.rds")
orig_pca <- orig_pca@cell.embeddings
colnames(orig_pca) <- paste0("Original.", colnames(orig_pca))
orig_pca <- orig_pca[ , 1:50]

scaled_pca <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling/GSC_GBM/ScaledBySample_GCC_GBM_pca.rds")
scaled_pca <- scaled_pca@cell.embeddings
colnames(scaled_pca) <- paste0("Scaled.", colnames(scaled_pca))

#correlate PCs between method
combined <- cbind(orig_pca, scaled_pca)
corr.mat <- cor(combined, method = "pearson")
corr.mat <- corr.mat[grep("Scaled", rownames(corr.mat)), grep("Original", rownames(corr.mat))]
saveRDS(corr.mat, file = "PearsonCorrelation_CellEmbeddings_GSC_GBM.rds")


### 6.2) Correlate gene loadings between processing methods
orig_pca <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling/GSC_GBM/Original_GCC_GBM_pca.rds")
orig_pca <- orig_pca@feature.loadings
colnames(orig_pca) <- paste0("Original.", colnames(orig_pca))
orig_pca <- orig_pca[ , 1:50]

scaled_pca <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling/GSC_GBM/ScaledBySample_GCC_GBM_pca.rds")
scaled_pca <- scaled_pca@feature.loadings
colnames(scaled_pca) <- paste0("Scaled.", colnames(scaled_pca))

combined <- cbind(orig_pca, scaled_pca)


corr.mat <- cor(combined, method = "pearson")
corr.mat <- corr.mat[grep("Scaled", rownames(corr.mat)), grep("Original", rownames(corr.mat))]
saveRDS(corr.mat, file = "PearsonCorrelation_CellEmbeddings_GSC_GBM.rds")



### 6.3) Overlap top
