##############################################################
#             Run Batch Correction Tools on GSCs             #
#                         L.Richards                         #
#                         June 2020                          #
##############################################################
### Reference: https://github.com/satijalab/seurat-wrappers
### Working dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load Global GSC data (corresponds to Ext Data Fig 4.)
### 2) Run CONOS wrapper in Seurat
### 3) Run LIGER in Seurat
### 4) Run fastMNN in Seurat
##############################################################
