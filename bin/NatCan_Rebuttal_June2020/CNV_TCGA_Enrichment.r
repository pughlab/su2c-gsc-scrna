##############################################################
#              CNV associations in TCGA datasets             #
#                         L.Richards                         #
#                         June 2020                          #
##############################################################
### Data downloaded from cBioportal:
### https://www.cbioportal.org/study/summary?id=gbm_tcga_pan_can_atlas_2018
### https://bit.ly/2YPL2SN
### 148 samples have both RNAseq (RSEM) and CNV (GISTIC)
###
### Pub Ref: https://pubmed.ncbi.nlm.nih.gov/29625048/



##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Define intersect of GBM samples with RNAseq and CNVs
### 2) Score GBMs with Dev and IR gene signatures with GSVA
### 3) Define which samples have -1, 0, 1 for chr arms
### 4) Compare Dev and IR scores for samples with those alterations
##############################################################



##############################################################
# 1) Install + Load packages
##############################################################
#library(devtools)
#devtools::install_github("RcppCore/Rcpp")
#devtools::install_github("satijalab/seurat-wrappers")
#devtools::install_github("hms-dbmi/conos")
#devtools::install_github('MacoskoLab/liger')
