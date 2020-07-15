##############################################################
#          CNV associations in public scRNA datasets         #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Split public data into normal brain and tumour
### 2) Run InferCNV using panel of normal + each sample
### 3) Plot CNV heatmap across cohort
### 4) Save per cell inferCNV scores
##############################################################

setwd("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA")


##############################################################
# 1) Load packages + data
##############################################################

library(Seurat)

##############################################################
# 1) Load packages + data
##############################################################

### Owen provided 2 public datasets:
### Bhaduri et al. --> cant use because no normal cells for refernece
### Wang et al. snRNA --> preference not to use because nuclei
### Wang et al. scRNA --> seems suitable

### Also downloaded 2 more public datasets:
### Neftel et al. -->
### Darmanis et al. -->


seurat.obj <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/public_data/Wang_scRNA_full_seurat.rds"

dat <- readRDS(seurat.obj)
