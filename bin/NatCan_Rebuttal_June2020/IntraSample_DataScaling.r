##############################################################
#               Scale Data within each sample first          #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: https://github.com/JinmiaoChenLab/Rphenograph
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/PCA_SampleScaling


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load GSC + tumour and GSC alone dataset
### 2) Split by sample
### 3) Scale each sample
### 4) Re-run PCA + clustering etc
### 5) Save
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
                  verbose = T,
                  features = rownames(split@assays$RNA@scale.data)
                  )

##############################################################
### 3) Re-run PCA + clustering in same way as manuscript
##############################################################
split <- RunPCA(split,
               features = pca.genes,
               npcs = 50,
               verbose = FALSE,
               reduction.key = "Scaled.PC_"
              )
split <- FindNeighbors(split)
split <- FindClusters(split, resolution = res.use)
split <- RunUMAP(split, dims = pc.use)


##############################################################
### 4) Correlate PCA results
##############################################################
orig_pca <- dat@reductions$pca
scaled_pca <- split@reductions$pca

### 4.1) Correlate PC cell embeddings
pca <- cbind(orig_pca@cell.embeddings[, 1:50], scaled_pca@cell.embeddings)
pca <- cor(pca, method = "pearson")
pca <- pca[grep("^Scaled", colnames(pca)), grep("^PC", colnames(pca)) ]
pc.corr.file <- paste0(outFileName, "_cellEmbeddings_correlation.rds")
saveRDS(pca, file = pc.corr.file)

### 4.2) Correlate PC gene loadings
genes <- cbind(orig_pca@feature.loadings[, 1:50], scaled_pca@feature.loadings)
genes <- cor(genes, method = "pearson")
genes <- genes[grep("^Scaled", colnames(genes), grep("^PC", colnames(genes)) ]
genes.corr.file <- paste0(outFileName, "_geneLoadings_correlation.rds")
saveRDS(genes, file = genes.corr.file)

### 4.3) Overlap (Jarccard) top 50 and bottom 50 genes / PC




##############################################################
### 5) Save
##############################################################

### 5.1) meta data


### 5.3) PCA data


### 5.4) Seurat Objects

meta <- split@meta.data
pca <- split@reductions$pca
saveRDS(meta, file = "ScaledBySample_GCC_GBM_meta.rds")
saveRDS(pca, file = "ScaledBySample_GCC_GBM_pca.rds")


saveRDS(split, file = "ScaledBySample_GCC_GBM_Seurat.rds")

meta_orig <- dat@meta.data
pca_orig <- dat@reductions$pca
saveRDS(meta_orig, file = "Original_GCC_GBM_meta.rds")
saveRDS(pca_orig, file = "Original_GCC_GBM_pca.rds")
