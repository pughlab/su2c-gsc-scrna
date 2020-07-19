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
library(infercnv) #v1.4 (cant re-install old version)

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


##############################################################
# 3) Run InferCNV (v0.3) on Samwise (bash scripts below)
##############################################################
### /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA
### Have to run on Mordor to use same version of InferCNV
### Cutoff of 1 works well for smartsseq2 data

module load R/3.2.2

#######################
# 3.1) Neftel et al.
#######################
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
--noise_filter 0.1 \
--output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel \
--vis_bound_threshold " -1,1" \
--log_file Neftel.log.txt \
--ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_Reference_Barcodes.txt \
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Neftel/Neftel_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt

#######################
# 3.2) Darmanis
#######################
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
--noise_filter 0.1 \
--output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis \
--vis_bound_threshold " -1,1" \
--log_file Darmanis.log.txt \
--ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_Reference_Barcodes.txt \
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/Darmanis/Darmanis_TumourOligo_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/GenePos_GRCh38.txt
