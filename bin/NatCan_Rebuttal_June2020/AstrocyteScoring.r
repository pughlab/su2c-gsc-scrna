##############################################################
#              Score GSCs for astrocyte signatres            #
#                         L.Richards                         #
#                         July 2020                     #
##############################################################
### Gene Sigs: https://docs.google.com/spreadsheets/d/1Ou8R_HvJyV8dFgXbK4v_3M10yGqW7nfJM3fvgPPzReI/edit?usp=sharing
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load Global GSC data (corresponds to Ext Data Fig 4.)
### 2) Run AUCell to score GSCs with astrocyte signatures
### 3) Correlate gene signatures to Dev and IR scores
### 4) Plot astrocyte gene sig activity bewteen Dev and IR ()
### 5) Plot Big heatmap to see which scores are most correlated ()
##############################################################
library(Seurat) #v3.1.5
library(AUCell) #v1.8.0
library(BBmisc)

##############################################################
# 1) Format gene signatures
##############################################################

### load (astrocyte) gene signatures (n=13)
astro <- read.csv("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AstrocyteSignatures_NatCan_July2020.csv")
astro <- convertColsToList(astro)
removeEmpty <- function(i){return(i[i != ""])}
removeDash <- function(i){return(i[i != "--"])}
astro <- lapply(astro, removeEmpty)
astro <- lapply(astro, removeDash)
str(astro)

### load global set of gene signatures (n=116)
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AUCell_Signatures_Hypoxia.Rdata")
sigs <- append(sigs, astro) #append astro

### save gene signatures
saveRDS(sigs, file = "./input_data/Gene_Signatures_wAstrocytes_July2020.rds")


##############################################################
# 2) Run AUCell
##############################################################

##############################################################
### EXAMPLE EXECUTION ON H4H
###
### #!/bin/bash
### #SBATCH -t 5:00:00
### #SBATCH --mem=100G
### #SBATCH -p veryhimem
### #SBATCH -c 30
### #SBATCH -N 1
### #SBATCH --account=pughlab
###
### module load R/3.6.1
###
### Rscript AstrocyteScoring.R
##############################################################

### load global GSC data
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
exprMatrix <- as.matrix(BTSC@data)
exprMatrix[1:10, 1:10]
dim(exprMatrix)

### load gene signatures
sigs <- readRDS("./input_data/Gene_Signatures_wAstrocytes_July2020.rds")

### run AUCell
cells_rankings <- AUCell_buildRankings(exprMatrix,
                                       nCores=30,
                                       plotStats=FALSE)

cells_AUC <- AUCell_calcAUC(sigs,
                            cells_rankings)

cells_assignment <- AUCell_exploreThresholds(cells_AUC,
                                             plotHist=FALSE,
                                             assign=TRUE
                                            )

AUC <- t(as.data.frame(cells_AUC@assays$data$AUC))
colnames(AUC) <- paste0(colnames(AUC), "_AUC")
AUC <- data.frame(AUC)
print(AUC[1:2, 1:2])

### merge AUCell scores to seuart meta.data
BTSC <- AddMetaData(BTSC, metadata = AUC)

### save AUCell scores + metadata
BTSC_AUCell <- BTSC@meta.data
saveRDS(BTSC_AUCell , file = "Astrocyte_AUCell_GSCs.rds")

#save the Seurat Object
save(BTSC, file = "Astrocyte_AUCell_GSCs_SeuratObj.rd")



##############################################################
# 3) Run AUCell
##############################################################
