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

##############################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 96:00:00
## #SBATCH --mem=100G
## #SBATCH -p veryhimem
## #SBATCH -c 1
## #SBATCH -N 1
## #SBATCH --account=pughlab
##
## module load R/3.6.1
##
## Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/IntraSample_DataScaling.r --dat.file /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata \
##  --outFileName GSC
##############################################################

### GSCs Only:
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata

### GSCs + Tumour:
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/SU2C_Live_GSC_GBM_Malignant_Seurat.rds

library(Seurat)
library(optparse)

##loads an RData file, and returns it with new name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

##############################################################
######## USER DEFINED PARAMS AND FUNCTIONS
#global GSCs
#dat.file <-
#outFileName <-
#GSC+Tumour
#dat.file <- "../../data/SU2C_Live_GSC_GBM_Malignant_Seurat.rds"
#outFileName <- "GSCs_GBM"
##############################################################

##############################################################
### PARSE OPTIONS
##############################################################
print("")
print("********************")
print("Parse options")
print(Sys.time())
option_list <- list(make_option("--dat.file",
                                type = "character",
                                default = NULL,
                                help = "path to Seurat object",
                                metavar= "character"
                               ),
                     make_option("--outFileName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files",
                                metavar= "character"
                              )
                            )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dat.file <- opt$dat.file
outFileName <- opt$outFileName

print(paste0("Dataset: ",dat.file))
print(paste0("Output file prefix: ", outFileName))



##############################################################
### 1) Load data
##############################################################
print("")
print("********************")
print("Loading Data")
print(Sys.time())

if(grepl(".rds$", dat.file)){
  dat <- readRDS(dat.file)
} else if (grepl(".data$", dat.file)){
  dat <- loadRData(dat.file)
}
### extract processing information from original object
pca.genes <- dat@calc.params$RunPCA$pc.genes
pc.use <- 1:30
res.use <-  2
dat <- UpdateSeuratObject(dat)


##############################################################
### 2) Split and Scale by Sample
##############################################################
print("")
print("********************")
print("Scale data within samples")
print(Sys.time())

split <- ScaleData(dat,
                  split.by = "SampleID",
                  do.center = TRUE,
                  vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                  verbose = T,
                  features = rownames(dat@assays$RNA@scale.data)
                  )
print(dim(split@assays$RNA@scale.data))
split.file <- paste0(outFileName, "_SplitData_unprocessed.rds")
print("Saving split intermediate file that takes way to long to generate.....")
print(split.file)
saveRDS(split, file = split.file)



##############################################################
### 3) Re-run PCA + clustering in same way as manuscript
##############################################################
print("")
print("********************")
print("Run PCA and Clustering")
print(Sys.time())

split <- RunPCA(split,
               features = pca.genes,
               npcs = 50,
               verbose = FALSE,
               reduction.key = "ScaledPC_"
              )
split <- FindNeighbors(split)
split <- FindClusters(split, resolution = res.use)
split <- RunUMAP(split, dims = pc.use)



##############################################################
### 5) Save Results
##############################################################
print("")
print("********************")
print("Save results")
print(Sys.time())

### 5.1) Seurat Objects
orig_seurat.file <- paste0(outFileName, "_Original_Seurat.rds")
print(orig_seurat.file)
saveRDS(dat, file = orig_seurat.file)

split_seurat.file <- paste0(outFileName, "_SplitScaled_Seurat.rds")
print(split_seurat.file)
saveRDS(split, file = split_seurat.file)


### 5.2) PCA runs
pca_orig <- dat@reductions$pca
pca_orig.file <- paste0(outFileName, "_Original_PCA.rds")
print(pca_orig.file)
saveRDS(pca_orig, file = pca_orig.file)

pca_split <- split@reductions$pca
pca_split.file <- paste0(outFileName, "_SplitScaled_PCA.rds")
print(pca_split.file)
saveRDS(pca_split, file = pca_split.file)


### 5.3) Metadata with UMAP coords
meta_orig <- dat@meta.data
meta_orig <- cbind(meta_orig, split@reductions$umap@cell.embeddings)
meta_orig.file <- paste0(outFileName, "_Original_meta.rds")
print(meta_orig.file)
saveRDS(meta_orig, file = meta_orig.file)

meta_split <- split@meta.data
meta_split <- cbind(meta_split, split@reductions$umap@cell.embeddings)
meta_split.file <- paste0(outFileName, "_SplitScaled_meta.rds")
print(meta_split.file)
saveRDS(meta_split, file = meta_split.file)


##############################################################
print("")
print("********************")
print("END OF SCRIPT")
print(Sys.time())
