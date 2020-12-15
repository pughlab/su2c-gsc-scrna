
----
# Inferring scRNA CNVs with CONICS
----
L.Richards  

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/CONICS/BrainCells

## 1.0 Run CONICS on Brain Cells Only


```R
### barcodes of cells with brain ID

Immune <- c("C45", "C50", "C24", "C10", "C33", "C40", "C15", "C47", "C26")
Endothelial <- "C42"

#dat$res.1.8 <- paste0("C", dat$res.1.8)

brain.cells <- rownames(dat[dat$res.1.8 %nin% c(Immune, Endothelial), ])
expr.dat <- as.matrix(merged@data)
expr.dat <- expr.dat[ ,brain.cells]
dim(expr.dat)

saveRDS(expr.dat, file = "BrainCells_47224_LogNorm_CONICS_input.rds")
```


```R
###write out brain cell type matrix logNormalized
########################################

library(CONICSmat)
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/")
                )
#USD
file.prefix <- "SU2C_NucSeq_BrainCells_LogNormal_CONICS"

########################################

print("###########")
print("Load Data")
print(Sys.time())
print("")

    #merged <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/Nuclei_Tumours_Merged_noG800_AUCell_Seurat.rds")
    

print("###########")
print("Extract LogNormalized Counts")
print(Sys.time())
print("")

    expr.dat <- readRDS("BrainCells_47224_LogNorm_CONICS_input.rds")


print("###########")
print("Run CONICS")
print(Sys.time())
print("")

    #make a patient string across all the cells
    patients <- gsub('.{17}$', '', colnames(expr.dat))
    print(unique(patients))

    #download and read in chr regions fro GRCH38 
    regions <- read.table("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/CNVs/chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
    head(regions, n=5)

    #get chromosome positions of gene in expression matrix
    #dont need additional gene filtering - we already did it
    gene_pos <- getGenePositions(rownames(expr.dat))
    head(gene_pos)

    #calculate normalization factors
    normFactor <- calcNormFactors(expr.dat)

    l <- plotAll(expr.dat,
             normFactor,
             regions,
             gene_pos,
             file.prefix
                )

    #plot heamtap
    
    #hi <- plotHistogram(l,
    #                    expr.dat,
    #                    clusters =1 , #pick two cluster model to separate tumour and normal cells (we expect two cell types)
    #                    zscoreThreshold = 4,
    #                    patients
    #                   )

        
        
print("###########")
print("Save Data")
print(Sys.time())
print("")
        
    save.file <- paste0(file.prefix, ".rds")
    print(save.file)
    saveRDS(l, file = save.file)

print("###########")
print("FIN")
print(Sys.time())
print("")
```

----
## 2.0 Explore results
----

Cell stem cell paper used the following cutoffs:
> Nominate chrs with somatic CNVsas those with (1) a significant deviation of the log-likelihood of the model compared to a one-component model (likelihood ratio test < 0.001) and (2) a difference in Bayesian Inference Criterion (BIC) > 300      
>  
> For somatic chrs, use a cutoff of posterior probability (pp > 0.8) to infer the presence/absence of the respective CNV in a cell   
>  
> Cell with CNV alterations were classified as tumor cells, whereas cells with statistical probabilities suggesting no CNV alterations were classified as normal cells.
>  
> Cells that could not be clearly assigned to a genotype (e.g., 0.2 < pp < 0.8) remained unclassified. 





```R
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggrepel)
library(dplyr)
library(Hmisc)

setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/CONICS/BrainCells")
```
logNormal <- read.table("SU2C_NucSeq_BrainCells_LogNormal_CONICS_BIC_LR.txt",
                    sep="\t",
                    header=T,
                    row.names=1,
                    check.names=F
                   )

colnames(logNormal) <- gsub(" ", "_", colnames(logNormal))
colnames(logNormal) <- gsub("-", "_", colnames(logNormal))

logNormal_pp <- readRDS("SU2C_NucSeq_BrainCells_LogNormal_CONICS.rds")
colnames(logNormal_pp) <- paste0("chr", colnames(logNormal_pp))
cells <- rownames(logNormal_pp)

## subset to only the 

CONICS_UMAP <- readRDS("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/Nuclei_Tumours_Merged_noG800_tSNE_UMAP_AUCell.rds")
head(CONICS_UMAP)

meta <- CONICS_UMAP[cells, c(colnames(CONICS_UMAP[1:15]), "UMAP1", "UMAP2") ]
meta <- cbind(meta, logNormal_pp)

saveRDS(meta, file = "CONICS_meta_BrainCells.rds")

```R
CONICS_UMAP <- readRDS("CONICS_meta_BrainCells.rds")
CONICS_UMAP$res.1.8 <- paste0("C", CONICS_UMAP$res.1.8)
```


```R
Oligodendrocytes <- c("C7", "C27", "C21")
Neurons <- c("C32", "C29")

#dat$res.1.8 <- paste0("C", dat$res.1.8)


cell.type <- c()
cell.type[CONICS_UMAP$res.1.8 %in% Oligodendrocytes] <- "Oligodendrocytes"
cell.type[CONICS_UMAP$res.1.8 %in% Neurons] <- "Neurons"
cell.type[CONICS_UMAP$res.1.8 %nin% c(Oligodendrocytes, Neurons)] <- "Tumour"
CONICS_UMAP$CellType <- cell.type
table(CONICS_UMAP$CellType)
```


    
             Neurons Oligodendrocytes           Tumour 
                1377             3513            42334 



```R
pdf("~/Desktop/Nuclei_CONICS_mean_BRAINONLY.pdf", width = 15, height = 6)

regions <- c("chr7p", "chr7q", "chr10p", "chr10q", "chr9p", "chr17q")


for (i in 1:length(regions)){

region <-  regions[i]
mid <- mean(CONICS_UMAP[,region])
mid
    
oligo <- CONICS_UMAP[CONICS_UMAP$CellType == "Oligodendrocytes", ]
neuron <- CONICS_UMAP[CONICS_UMAP$CellType == "Neurons", ]
tumour <- CONICS_UMAP[CONICS_UMAP$CellType == "Tumour", ]

oligo_UMAP <- ggplot(oligo, aes(x=UMAP1, y=UMAP2, color=oligo[,region])) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab", name = "Posterior Probability \n Z-score") + 
                ggtitle(paste0("Oligodendrocytes: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))


neuron_UMAP <- ggplot(neuron, aes(x=UMAP1, y=UMAP2, color=neuron[,region])) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                     high="red", space ="Lab", name = "Posterior Probability \n Z-score") + 
                ggtitle(paste0("Neurons: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))




tumour_UMAP <- ggplot(tumour, aes(x=UMAP1, y=UMAP2, color=tumour[,region])) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                     high="red", space ="Lab", name = "Posterior Probability \n Z-score") + 
                ggtitle(paste0("Tumour: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))


print(ggarrange(oligo_UMAP, 
          neuron_UMAP, 
          tumour_UMAP, 
          ncol = 3, 
          common.legend = TRUE
         ))


}

dev.off()
```


<strong>pdf:</strong> 2



```R
pdf("~/Desktop/Nuclei_CONICS_mean_BRAINONLY_pp0.8.pdf", width = 15, height = 6)

regions <- c("chr7p", "chr7q", "chr10p", "chr10q", "chr9p", "chr17q")


for (i in 1:length(regions)){

region <-  regions[i]
mid <- mean(CONICS_UMAP[,region])
mid
    
oligo <- CONICS_UMAP[CONICS_UMAP$CellType == "Oligodendrocytes", ]
neuron <- CONICS_UMAP[CONICS_UMAP$CellType == "Neurons", ]
tumour <- CONICS_UMAP[CONICS_UMAP$CellType == "Tumour", ]

oligo_UMAP <- ggplot(oligo, aes(x=UMAP1, y=UMAP2, color=oligo[,region] >= 0.8)) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_manual(values=c("grey", "black")) +
                ggtitle(paste0("Oligodendrocytes: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))


neuron_UMAP <- ggplot(neuron, aes(x=UMAP1, y=UMAP2, color=neuron[,region] >= 0.8)) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_manual(values=c("grey", "black")) +
                ggtitle(paste0("Neurons: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))




tumour_UMAP <- ggplot(tumour, aes(x=UMAP1, y=UMAP2, color=tumour[,region] >= 0.8)) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_manual(values=c("grey", "black")) +
                ggtitle(paste0("Tumour: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))


print(ggarrange(oligo_UMAP, 
          neuron_UMAP, 
          tumour_UMAP, 
          ncol = 3, 
          common.legend = TRUE
         ))


}

dev.off()
```


<strong>pdf:</strong> 2



```R


tumour_UMAP <- ggplot(tumour, aes(x=UMAP1, y=UMAP2, color=tumour[,"chr10q"] <= 0.4)) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_manual(values=c("grey", "black")) +
                ggtitle(paste0("Tumour: ", region)) + theme_bw() + 
                 theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        xlim(c(min(CONICS_UMAP$UMAP1), max(CONICS_UMAP$UMAP1))) +
                        ylim(c(min(CONICS_UMAP$UMAP2), max(CONICS_UMAP$UMAP2)))
tumour_UMAP
```




![png](output_11_1.png)


### Compare CNV signal across cells


```R
CONICS_UMAP$TumourSample <- ifelse(CONICS_UMAP$CellType == "Tumour", 
                                   as.character(CONICS_UMAP$SampleID), 
                                   CONICS_UMAP$CellType
                                  )
```


```R
par(mfrow=c(1,2))
boxplot(CONICS_UMAP$chr10p ~ CONICS_UMAP$TumourSample, las = 2, main = "chr10p")
boxplot(CONICS_UMAP$chr7q ~ CONICS_UMAP$TumourSample, las = 2, main = "chr7q")
boxplot(CONICS_UMAP$chr17q ~ CONICS_UMAP$TumourSample, las = 2, main = "chr17q")
```


![png](output_14_0.png)



![png](output_14_1.png)


----
## Run InferCNV
----

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/inferCNV

### 1.0 Prepare input files


```R
### prepare input files

suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))

## load in data
merged <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/Nuclei_Tumours_Merged_noG800_AUCell_Seurat.rds")
meta <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/Nuclei_Tumours_Merged_noG800_tSNE_UMAP_AUCell_CellType.rds")

merged <- AddMetaData(merged, 
                      metadata = meta[ ,c("UMAP1", "UMAP2", "CellType")]
                     )

## extract out the raw counts for each cell type

oligo_raw <- merged@raw.data[ ,rownames(merged@meta.data[merged@meta.data$CellType == "Oligodendrocytes", ])]
dim(oligo_raw)

tumour_raw <- merged@raw.data[ ,rownames(merged@meta.data[merged@meta.data$CellType == "Tumour", ])]
dim(tumour_raw)

neuron_raw <- merged@raw.data[ ,rownames(merged@meta.data[merged@meta.data$CellType == "Neurons", ])]
dim(neuron_raw)


## combine subsets 
## 47,224 cells

expr.mat <- cbind(oligo_raw,
                  neuron_raw,
                  tumour_raw
                 )

cnv.dat <- LogNormalize(data = expr.mat,
                        scale.factor = 100000
                        )

#save input files
saveRDS(cnv.dat, file = "SU2C_Nuclei_TumourOligoNeuron_cnvinput.rds") #normlized data as R object

write.table(as.matrix(cnv.dat), 
            file = "SU2C_Nuclei_TumourOligoNeuron_cnvinput.txt", 
            sep = "\t", 
            col.names = T, 
            row.names = T, 
            quote = F
           )

#---Write out comma delinited reference cell barcode file---#

oligos <- colnames(oligo_raw)
neurons <- colnames(neuron_raw)
oligo_neuron <- c(oligos, neurons)

write.table(matrix(as.character(oligos),nrow=1), 
            file = "Oligo_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

write.table(matrix(as.character(neurons),nrow=1), 
            file = "Neurons_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

write.table(matrix(as.character(oligo_neuron),nrow=1), 
            file = "OligoNeuron_Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

```

### 2.0 Run inferCNV
#!/bin/bash
#
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 48:00:00

###################################
### Run 'InferCNV' on snRNA-seq
### L.Richards
### February 2020
###################################

### Example Usage:
### sbatch -p himem --mem 30G runInferCNV.sh REFERENCE_NAME WORK_DIR
### sbatch -p himem --mem 30G runInferCNV.sh Oligo /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/inferCNV

cd $2

module load R/3.2.2

mkdir $1_reference

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 --noise_filter 0.1 --output_dir $2/$1_reference --vis_bound_threshold " -1,1" --log_file SU2C_NucSeq_$1_Reference_log.txt --ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/inferCNV/$1_Reference_Barcodes.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NucSeq_Integration/merge/malignant_classification/inferCNV/SU2C_Nuclei_TumourOligoNeuron_cnvinput.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/tumours/input/GenePos_GRCh38.txt