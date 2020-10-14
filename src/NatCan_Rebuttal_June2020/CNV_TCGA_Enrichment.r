##############################################################
#              CNV associations in TCGA datasets             #
#                         L.Richards                         #
#                         June 2020                          #
##############################################################
### Data downloaded from cBioportal (hg19):
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

setwd("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_TCGA/")



##############################################################
# 1) Load packages
##############################################################

library(GSVA)
library(CONICSmat)



##############################################################
# 2) Preprocess Data
##############################################################

### 2.1) Read in RNAseq data in the form of RSEM
### 162 samples, 20501 genes with hugo symbols
rnaseq <- read.table("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_TCGA/data/gbm_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt", sep = "\t", header = T)
rnaseq <- rnaseq[!rnaseq$Hugo_Symbol == "", ] #remove empty hugp symbols
rnaseq <- rnaseq[!duplicated(rnaseq$Hugo_Symbol), ] #remove duplicated genes (n=1)
rownames(rnaseq) <- rnaseq$Hugo_Symbol
rnaseq$Hugo_Symbol <- NULL
rnaseq$Entrez_Gene_Id <- NULL

### 2.2) Read in cnv data in form of GISTIC gene-level calls
###
cnvs <- read.table("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_TCGA/data/gbm_tcga_pan_can_atlas_2018/data_CNA.txt", sep = "\t", header = T)
cnvs <- cnvs[!cnvs$Hugo_Symbol == "", ] #remove empty hugp symbols
cnvs <- cnvs[!duplicated(cnvs$Hugo_Symbol), ]
rownames(cnvs) <- cnvs$Hugo_Symbol
cnvs$Hugo_Symbol <- NULL
cnvs$Entrez_Gene_Id <- NULL

### 2.3) Create genomic posistion file
### Need to average all gene CNV scores across arm
### hg19 cytobands downloaded from: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
### Used esembl to retireve gene positions for all cnv genes
### ^http://useast.ensembl.org/biomart/martview/fdcf4f31143ae0c80b08e411e55a5250

### need to run this step on the build parition (requires internet)
### run in home directory on H4H
### cnvs <- read.table("data_CNA.txt", sep = "\t", header = T)
### gene.pos <- getGenePositions(gene_names = as.character(cnvs$Hugo_Symbol))
### saveRDS(gene.pos, file = "TCGA_GBM_CNV_GenePos.rds")
### Downloaded chr arm position file from CONICS github repo
genepos <- readRDS("./data/TCGA_GBM_CNV_GenePos.rds")
chrarms <- read.table("./data/chr_arms.txt", header = T)
genepos <- genepos[!genepos$chromosome_name == 0 ,] #remove genes mapped to chr0
genepos <- genepos[!genepos$chromosome_name == 23 ,]
#map gene ends to chr arm

gene.arm <- c()

for (i in 1:nrow(genepos)){

    arm.subset <- chrarms[chrarms$Chrom == genepos$chromosome_name[i] ,]
    gene.arm[i] <- ifelse(genepos$end_position[i] > arm.subset[1 ,"End"],
           as.character(arm.subset$Idf[2]), #greater than p arm end = q
           as.character(arm.subset$Idf[1]) #NOT greater than p arm end = p
         )
}
genepos$arm <- gene.arm
saveRDS(genepos, file = "./data/TCGA_GBM_CNV_GenePos_arm.rds")

### 2.4) Subset rnaseq and cnv dfs to common samples
common <- intersect(colnames(rnaseq), colnames(cnvs)) #148 = matched cbio
cnvs <- cnvs[ ,common]
dim(cnvs)
rnaseq <- rnaseq[ ,common]
dim(rnaseq)

saveRDS(cnvs, file = "./data/TCGA_GBM_CNV_processed.rds")
saveRDS(rnaseq, file = "./data/TCGA_GBM_RNA_processed.rds")

rm(list = ls()) #clean up environ



##############################################################
# 3) Score RNAseq data with Dev and IR gene signatures
##############################################################

rnaseq <- readRDS("./data/TCGA_GBM_RNA_processed.rds")
rnaseq <- as.matrix(rnaseq)
load("./data/AUCell_Signatures_Hypoxia.Rdata") #object called sigs
sigs <- sigs[c("Developmental_GSC", "InjuryResponse_GSC")]

gsva.scores <- gsva(expr = rnaseq,
                    gset.idx.list = sigs,
                    method = "gsva"
                  )
gsva.scores <- data.frame(t(gsva.scores))
saveRDS(gsva.scores, file = "TCGA_GBM_RNA_GSVA.rds")

## plot results
pdf("TCGA_GBM_RNA_GSVA.pdf", width = 6, height = 6)
plot(gsva.scores$Developmental_GSC,
     gsva.scores$InjuryResponse_GSC,
     xlab = "Developmental Score",
     ylab = "Injury Response Score",
     main = "TCGA GBM (n=148)"
   )
abline(h=0, v=0, lty =2, col = "red")
dev.off()



##############################################################
# 4) Classify samples by chr arm
##############################################################
### For each chr arm, we want to bin samples into loss, neutral and gain
### Do we average or take median across genes

cnvs <- readRDS("./data/TCGA_GBM_CNV_processed.rds")
arm <- readRDS("./data/TCGA_GBM_CNV_GenePos_arm.rds")

### 4.1) Classify
### Neutral = 0
### Gain = >=1
### Loss = <=-1
cnvs[cnvs > 0] <- 1 #treat high and low level the same
cnvs[cnvs < 0]<- (-1)


arm_avg <- list() #average scores across arm
arm_med <- list()#take median score across arm
arm_mod <- list() #most common value = score
arm_50 <-  list() #classify based on 50% of calls (ie. 50% del or 50% gain to be classified)
chr.counts <- list() #counts of 0, 1, -1 per sample per arm

### function to get mode from a numeric vector
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

chr.arms <- unique(arm$arm)

for(i in 1:length(chr.arms)){

    print(chr.arms[i])
    genes <- arm[arm$arm == chr.arms[i], ]
    genes <- as.character(genes$hgnc_symbol) #get genes on chr arm

    name <- paste0("chr", chr.arms[i])
    cnv.subset <- cnvs[genes, ] #subset cnv data to chr arm
    arm_avg[[name]] <- colMeans(cnv.subset) #mean
    arm_med[[name]] <- apply(cnv.subset, 2, FUN = median) #median
    arm_mod[[name]] <- apply(cnv.subset, 2, FUN = getmode) #mode

    freq <- list()
    for(j in 1:ncol(cnv.subset)){ #count values for each sample

       counts <- c(as.numeric(table(cnv.subset[ ,j] < 0)["TRUE"]),
                    as.numeric(table(cnv.subset[ ,j] == 0)["TRUE"]),
                    as.numeric(table(cnv.subset[ ,j] > 0)["TRUE"])
                  )
       counts[is.na(counts)] <- 0
       names(counts) <- c("Deletion", "Neutral", "Gain")
       freq[[colnames(cnv.subset)[j]]] <- c(counts)

    }

    chr.counts[[name]] <- freq

}

master.counts <- list()
for(x in 2:length(names(chr.counts))){

    chr.name <- names(chr.counts)[x]
    print(chr.name)
    subset <- chr.counts[[x]]
    subset <- do.call(rbind, subset)
    subset <- data.frame(subset)
    subset$Chr <- as.character(names(chr.counts)[x])
    subset$Sample <- row.names(subset)
    master.counts[[chr.name]] <- subset

}
master.counts <- do.call(rbind, master.counts)
arm_avg <- do.call(rbind, arm_avg)
arm_med <- do.call(rbind, arm_med)
arm_mod <- do.call(rbind, arm_mod)

saveRDS(master.counts, file = "TCGA_GBM_CNV_Counts.rds")
saveRDS(arm_avg, file = "TCGA_GBM_CNV_armAvg.rds")
saveRDS(arm_med, file = "TCGA_GBM_CNV_armMed.rds")
saveRDS(arm_mod, file = "TCGA_GBM_CNV_armMode.rds")

### 4.2) Determine optimal CNV arm classification
### see jupyter notebook TCGA_CNV_Plotting.ipynb

### merge GSVA results with CNV results
avg <- t(avg)
avg <- round(avg)
