##############################################################
#          CNV associations in public scRNA datasets         #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Working Dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Split public data into normal brain and tumour
### 2) Run InferCNV using panel of normal + each sample
### 3) Plot CNV heatmap across cohort
### 4) Save per cell inferCNV scores
##############################################################

### Owen provided 2 public datasets:
### Bhaduri et al. --> cant use because no normal cells for refernece
### Wang et al. snRNA --> prefer not to use because nuclei
### Wang et al. scRNA --> seems suitable
### I dont want to use these until Owen updates his tumour/normal methods to match the paper

### I have already processed 2 public datasets as part of OICR_Brain_NucSeq SingleR reference creation
### Move forward with these two datasets for now..
### Darmanis et al --> /cluster/projects/pughlab/projects/OICR_Brain_NucSeq/GBM/analysis/SingleR/references/scRNA/GBM_Darmanis_CellRep_2017
### Neftel et al --> /cluster/projects/pughlab/projects/OICR_Brain_NucSeq/GBM/analysis/SingleR/references/scRNA/GBM_Neftel_Cell_2019

setwd("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA")


##############################################################
# 1) Load packages
##############################################################
library(Seurat)
#library(infercnv) #v1.4 (cant re-install old version)
library(AUCell)
library(scater)
library(data.table)

##############################################################
# 2) Extract tumour and normal log matrices
##############################################################
### These will be input for inferCNV

############################
# 2.1) Daramanis et al.,
############################
### Singlecellexperiment object
sce <- readRDS("/cluster/projects/pughlab/projects/OICR_Brain_NucSeq/GBM/analysis/SingleR/references/scRNA/GBM_Darmanis_CellRep_2017/GBM_Darmanis_CellRep_2017_Ref.rds")

### convert to seurat object
### subset to malignant or reference and output matriced
seurat <- as.Seurat(sce,
                    counts = "counts",
                    data = "logcounts"
                  )
cells <- c(rownames(seurat@meta.data[seurat@meta.data$label == "Malignant_GBM_Adult", ]),
rownames(seurat@meta.data[seurat@meta.data$label == "Oligodendrocyte", ]),
rownames(seurat@meta.data[seurat@meta.data$label == "OPC", ])
)
data <- seurat@assays$RNA@data[ ,cells] ### 23368 genes x 1582 cells

### used for inferCNV input:
### (1) lognormal matrix of tumour+reference cells
DGE.name <- "Darmanis_TumourOligo_DGE.txt"
write.table(as.matrix(data),
        file = DGE.name,
        sep = "\t",
        col.names = T,
        row.names = T,
        quote = F
       )

### (2) comman delimited refernece cell barcode file
ref.bc <- c(rownames(seurat@meta.data[seurat@meta.data$label == "Oligodendrocyte", ]),
rownames(seurat@meta.data[seurat@meta.data$label == "OPC", ]))
write.table(matrix(as.character(ref.bc),nrow=1),
            file = "Darmanis_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

### (3) Save tumour only cells for AUCell
cells <- c(rownames(seurat@meta.data[seurat@meta.data$label == "Malignant_GBM_Adult", ])
)
data <- seurat@assays$RNA@data[ ,cells] ### 23368 genes x 1582 cells
meta <- seurat@meta.data[cells, ]
saveRDS(data, file = "Darmanis_TumourOnly_LogCounts.rds")
saveRDS(meta, file = "Darmanis_TumourOnly_meta.rds")

########################
# 2.2) Neftel et al.,
########################
### Singlecellexperiment object
sce <- readRDS("/cluster/projects/pughlab/projects/OICR_Brain_NucSeq/GBM/analysis/SingleR/references/scRNA/GBM_Neftel_Cell_2019/GBM_Neftel_Cell_2019_Ref.rds")
logcounts(sce) <- as.matrix(logcounts(sce))

### convert to seurat object
### subset to malignant or reference and output matriced
seurat <- as.Seurat(sce,
                    counts = "logcounts",
                    data = "logcounts"
                    )
### subset to only ADULT tumours (no peds)
patients <- names(table(seurat@meta.data[seurat@meta.data$label == "Malignant_GBM_Adult", ]$Sample))
meta.sub <- seurat@meta.data[seurat@meta.data$Sample %in% patients, ]

cells <- c(rownames(meta.sub[meta.sub$label == "Malignant_GBM_Adult", ]),
        rownames(meta.sub[meta.sub$label == "Oligodendrocyte", ]) ###only want normal from matched tumours
           )

data <- seurat@assays$RNA@data[ ,cells] ### 23686 genes x 5126 cells

### used for inferCNV input:
### (1) lognormal matrix of tumour+reference cells
DGE.name <- "Neftel_TumourOligo_DGE.txt"
write.table(as.matrix(data),
                   file = DGE.name,
                   sep = "\t",
                   col.names = T,
                   row.names = T,
                   quote = F
                  )

### (2) comman delimited refernece cell barcode file
ref.bc <- c(rownames(meta.sub[meta.sub$label == "Oligodendrocyte", ]))
write.table(matrix(as.character(ref.bc),nrow=1),
                       file = "Neftel_Reference_Barcodes.txt",
                       sep=",",
                       row.names=FALSE,
                       col.names=FALSE,
                       quote = F
                      )

### (3) save tumour only expression matrix for AUCell
cells <- c(rownames(meta.sub[meta.sub$label == "Malignant_GBM_Adult", ]))
data <- seurat@assays$RNA@data[ ,cells]
meta <- meta.sub[cells, ]
saveRDS(data, file = "Neftel_TumourOnly_LogCounts.rds")
saveRDS(meta, file = "Neftel_TumourOnly_meta.rds")


########################
# 2.3) Wang et al.,
########################
seurat <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/public_data/Wang/Wang_scRNA_full_seurat.rds")

### add labels from Owens analysis to data
seurat@meta.data$cluster <- gsub(0, "Tumour", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(1, "Tumour", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(2, "Tumour", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(3, "Macrophage", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(4, "Oligodendrocytes", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(5, "Tumour", seurat@meta.data$cluster)
seurat@meta.data$cluster <- gsub(6, "Endothelial", seurat@meta.data$cluster)

cells <- c(rownames(seurat@meta.data[seurat@meta.data$cluster == "Tumour", ]),
        rownames(seurat@meta.data[seurat@meta.data$cluster == "Oligodendrocytes", ]) ###only want normal from matched tumours
           )

data <- seurat@data[ ,cells] ### 15794 genes x 15827 cells

### used for inferCNV input:
### (1) lognormal matrix of tumour+reference cells
DGE.name <- "Wang_scRNA_TumourOligo_DGE.txt"
write.table(as.matrix(data),
                   file = DGE.name,
                   sep = "\t",
                   col.names = T,
                   row.names = T,
                   quote = F
                  )

### (2) comman delimited refernece cell barcode file
ref.bc <- c(rownames(seurat@meta.data[seurat@meta.data$cluster == "Oligodendrocytes", ]))
write.table(matrix(as.character(ref.bc),nrow=1),
            file = "Wang_scRNA_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
            )

### (3) save tumour only expression matrix for AUCell
cells <- c(rownames(seurat@meta.data[seurat@meta.data$cluster == "Tumour", ]))
data <- seurat@data[ ,cells]
meta <- seurat@meta.data[cells, ]
saveRDS(data, file = "Wang_TumourOnly_LogCounts.rds")
saveRDS(meta, file = "Wang_TumourOnly_meta.rds")


########################
# 2.4) Richards et al.,
########################

loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

seurat <- loadRData("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/data/SU2C_Live_GBM_AllCells_Seurat_Oct2019.Rdata")

cells <- c(rownames(seurat@meta.data[seurat@meta.data$CellType == "Tumour", ]),
        rownames(seurat@meta.data[seurat@meta.data$CellType == "NormalBrain", ]) ###only want normal from matched tumours
           )

data <- seurat@data[ ,cells] ### 17647 x 17982 cells

### used for inferCNV input:
### (1) lognormal matrix of tumour+reference cells
DGE.name <- "Richards_scRNA_TumourOligo_DGE.txt"
write.table(as.matrix(data),
                   file = DGE.name,
                   sep = "\t",
                   col.names = T,
                   row.names = T,
                   quote = F
                  )

### (2) comman delimited refernece cell barcode file
ref.bc <- c(rownames(seurat@meta.data[seurat@meta.data$CellType == "NormalBrain", ]))
write.table(matrix(as.character(ref.bc),nrow=1),
            file = "Richards_scRNA_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
            )

### (3) save tumour only expression matrix for AUCell
cells <- c(rownames(seurat@meta.data[seurat@meta.data$CellType == "Tumour", ]))
data <- seurat@data[ ,cells]
meta <- seurat@meta.data[cells, ]
saveRDS(data, file = "Richards_TumourOnly_LogCounts.rds")
saveRDS(meta, file = "Richards_TumourOnly_meta.rds")


##############################################################
# 3) Run InferCNV (v0.3) on Samwise (bash scripts below)
##############################################################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA
### Have to run on Mordor to use same version of InferCNV
### Cutoff of 1 works well for smartsseq2 data
### Transferred back to H4H

### module load R/3.2.2

#######################
# 3.1) Neftel et al.
#######################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 1 \
### --noise_filter 0.1 \
### --output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel \
### --vis_bound_threshold " -1,1" \
### --log_file Neftel.log.txt \
### --ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_Reference_Barcodes.txt \
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt

#######################
# 3.3) Wang
#######################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
### --noise_filter 0.1 \
### --output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Wang \
### --vis_bound_threshold " -1,1" \
### --log_file Wang_scRNA.log.txt \
### --ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Wang/Wang_scRNA_Reference_Barcodes.txt \
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Wang/Wang_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt

#######################
# 3.4) Richards
#######################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
### --noise_filter 0.1 \
### --output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Richards \
### --vis_bound_threshold " -1,1" \
### --log_file Richards.log.txt \
### --ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Richards/Richards_scRNA_Reference_Barcodes.txt \
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Richards/Richards_scRNA_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt

#######################
# 3.2) Darmanis
#######################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 1 \
### --noise_filter 0.1 \
### --output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis \
### --vis_bound_threshold " -1,1" \
### --log_file Darmanis.log.txt \
### --ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_Reference_Barcodes.txt \
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt


##############################################################
# 4) Score public datasets with AUCell
##############################################################

### load gene signatures
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/AstrocyteScoring/input_data/AUCell_Signatures_Hypoxia.Rdata")
sigs <- sigs[c("Developmental_GSC", "InjuryResponse_GSC")]

### load public datasets
darmanis <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_TumourOnly_LogCounts.rds")
darmanis_meta <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_TumourOnly_meta.rds")

neftel <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_TumourOnly_LogCounts.rds")
neftel_meta <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_TumourOnly_meta.rds")

wang <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Wang/Wang_TumourOnly_LogCounts.rds")
wang_meta <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Wang/Wang_TumourOnly_meta.rds")

richards <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Richards/Richards_TumourOnly_LogCounts.rds")
richards_meta <- readRDS("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Richards/Richards_TumourOnly_meta.rds")

### run AUCell
name <- "Wang"
meta <- wang_meta
exprMatrix <- as.matrix(wang)
#name <- "Darmanis"
#meta <- darmanis_meta
#exprMatrix <- as.matrix(darmanis)
exprMatrix[1:10, 1:10]
dim(exprMatrix)
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(sigs,cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)

AUC <- t(as.data.frame(cells_AUC@assays@data$AUC))
colnames(AUC) <- paste0(colnames(AUC), "_AUC")
AUC <- data.frame(AUC)
print(AUC[1:2, 1:2])
meta <- cbind(meta, AUC)

### merge AUCell scores to seuart meta.data
save.file <- paste0("./", name, "/", name, "_AUCell_Scores.rds")
saveRDS(meta, file = save.file)


##############################################################
# 5) Clean up CNV calls
##############################################################

### Order the genes by genome position and order
genePos <- read.table("GenePos_GRCh38.txt")
head(genePos)

### load in cnv data for cutoff=1
#name <- "Darmanis"
name <- "Richards"
#input.file <- paste0("./", name, "/Cutoff_1/expression_post_viz_transform.txt")
input.file <- paste0("./", name, "/cutoff_0.1/expression_post_viz_transform.txt")
obs <- data.table::fread(input.file)
obs <- data.frame(obs)
rownames(obs) <- obs$V1
obs <- obs[ , -1]
obs[1:5, 1:5]
dim(obs)
max(obs)
min(obs)

#get all genes in CNV plot from genePos
CNV.genes <- genePos[genePos$V1 %in% rownames(obs), ]
CNV.genes <- CNV.genes[ ,1:2]
rownames(CNV.genes) <- CNV.genes[,1]
CNV.genes[,1] <- NULL
colnames(CNV.genes) <- c("Chromosome")

order <- c(grep("^1$", CNV.genes$Chromosome),
grep("^2$", CNV.genes$Chromosome),
grep("^3$", CNV.genes$Chromosome),
grep("^4$", CNV.genes$Chromosome),
grep("^5$", CNV.genes$Chromosome),
grep("^6$", CNV.genes$Chromosome),
grep("^7$", CNV.genes$Chromosome),
grep("^8$", CNV.genes$Chromosome),
grep("^9$", CNV.genes$Chromosome),
grep("^10$", CNV.genes$Chromosome),
grep("^11$", CNV.genes$Chromosome),
grep("^12$", CNV.genes$Chromosome),
grep("^13$", CNV.genes$Chromosome),
grep("^14$", CNV.genes$Chromosome),
grep("^15$", CNV.genes$Chromosome),
grep("^16$", CNV.genes$Chromosome),
grep("^17$", CNV.genes$Chromosome),
grep("^18$", CNV.genes$Chromosome),
grep("^19$", CNV.genes$Chromosome),
grep("^20$", CNV.genes$Chromosome),
grep("^21$", CNV.genes$Chromosome),
grep("^22$", CNV.genes$Chromosome)
           )

CNV.genes$Number <- 1:length(CNV.genes$Chromosome)
CNV.genes <- CNV.genes[match(order, CNV.genes$Number),]
CNV.genes$Number <- NULL
head(CNV.genes)

##order expression matrix by this
order <- rownames(CNV.genes)
obs <- obs[order, ]
#save.file <- paste0("./", name, "/Cutoff_1/", name, "_Ordered_CNV_matrix.rds")
save.file <- paste0("./", name, "/cutoff_0.1/", name, "_Ordered_CNV_matrix.rds")
saveRDS(obs, file = save.file)

#save.file2 <- paste0("./", name, "/Cutoff_1/", name, "_CNVgenes.rds")
save.file2 <- paste0("./", name, "/cutoff_0.1/", name, "_CNVgenes.rds")
saveRDS(CNV.genes, file = save.file2)
