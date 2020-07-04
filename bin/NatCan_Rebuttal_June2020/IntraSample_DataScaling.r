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
