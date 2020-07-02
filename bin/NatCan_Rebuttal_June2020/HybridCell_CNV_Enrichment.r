##############################################################
#              CNV associations in scRNA with hybrids        #
#                         L.Richards                         #
#                         June 2020                          #
##############################################################
### Ref Jupyter Notebooks:
### Identify_C1C2_CNVs_G800Lremoved_Nov2019
### NoG800_chisquared_subtypeCNV_Dec2019



##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Classify hybrid cells as Dev or IR
### 2)
##############################################################

library(Seurat)
library(effsize) #v0.8.0
library(ggplot2)
library(reshape)
library(ggpubr)

setwd("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_HybridCells")

##############################################################
# 1) Load data and classify hybrid cells
##############################################################
### Contains 3 objects:
### BTSC.CNVs = cnvs at sc level
### avg.cnv.df = cluster cnvs
### CNV.genes = gene ordering
### meta = AUCell scores and metadata
load("./data/GlobalBTSC_CNVs_NoG800L.Rdata")
### Load CNV data averged by chr arm per cell
arm.CNVs <- readRDS("./data/Cell_armlevelCNVs_Dec222019.rds")
arm.CNVs <- arm.CNVs[ ,-c(1:6)]
colnames(arm.CNVs) <- paste0("chr", colnames(arm.CNVs))

### Classify cells
dat <- meta[ ,c("orig.ident", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")]
colnames(dat) <- c("Sample", "C1_AUC", "C2_AUC")

dat$Dev <- dat$C1_AUC > 0.11 #AUCell cutoffs
dat$IR <- dat$C2_AUC > 0.2 #AUCell cutodff
dat$ID <- paste0(dat$Dev, dat$IR)

dat$ID <- gsub("TRUEFALSE", "Dev", dat$ID)
dat$ID <- gsub("FALSETRUE", "IR", dat$ID)
dat$ID <- gsub("FALSEFALSE", "LowLow", dat$ID)
dat$ID <- gsub("TRUETRUE", "HiHi", dat$ID)
head(dat)

### Re-classify hybrid cells using highest score
### New cell numbers = 25,295 Dev + 40,360 IR
for (i in 1:nrow(dat)){
    cell <- dat$ID[i]
    if(dat$ID[i] %in% c("HiHi" , "LowLow")){
       dat$ID[i] <- ifelse(dat$C1_AUC[i] > dat$C2_AUC[i],
                          print("Dev"),
                          print("IR")
                          )
    }
}



##############################################################
# 2) Compare CNV signal between bins
##############################################################

chrs <- colnames(arm.CNVs)[grep("chr", colnames(arm.CNVs))]
boxplot.dat <- cbind(dat, arm.CNVs)
effect.size <- c()

pdf("With_Hybrids_chrArm_comparison.pdf", width = 4, height = 4)

for(i in 1:length(chrs)){

  chr <- chrs[i]
  print(chr)

  Dev <- boxplot.dat[grep("Dev", boxplot.dat$ID), chr]
  IR <- boxplot.dat[grep("IR", boxplot.dat$ID), chr]

  b <- cohen.d(Dev, #effect size
        IR,
        hedges.correction=TRUE
       )

       p <- ggviolin(boxplot.dat,
              x = "ID",
              y = chr,
              color = "ID",
              palette = c("red", "black"),
              add = "mean_sd",
              main = paste0("Hedges g: ", round(b$estimate, 3), "\nMagnitude: ", b$magnitude),
              legend = "none"
              ) + stat_compare_means() #  Add p-value
              effect.size[i] <- b$estimate
              print(p)
}

dev.off()



##############################################################
# 3) Get Effect sizes for analysis WITHOUT hybrid
##############################################################

dat2 <- meta[ ,c("orig.ident", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")]
colnames(dat2) <- c("Sample", "C1_AUC", "C2_AUC")

dat2$Dev <- dat2$C1_AUC > 0.11 #AUCell cutoffs
dat2$IR <- dat2$C2_AUC > 0.2 #AUCell cutodff
dat2$ID <- paste0(dat2$Dev, dat2$IR)

dat2$ID <- gsub("TRUEFALSE", "Dev", dat2$ID)
dat2$ID <- gsub("FALSETRUE", "IR", dat2$ID)
dat2$ID <- gsub("FALSEFALSE", "LowLow", dat2$ID)
dat2$ID <- gsub("TRUETRUE", "HiHi", dat2$ID)
head(dat2)

#now remove hybrid cells
dat2 <- dat2[!dat2$ID == "HiHi", ]
dat2 <- dat2[!dat2$ID == "LowLow", ]
table(dat2$ID)

#filter arm.cnvs
arm.CNVs_2 <- arm.CNVs[rownames(dat2), ]
dim(arm.CNVs_2)

chrs <- colnames(arm.CNVs_2)[grep("chr", colnames(arm.CNVs_2))]
chrs
boxplot.dat2 <- cbind(dat2, arm.CNVs_2)
effect.size2 <- c()

#pdf("With_Hybrids_chrArm_comparison.pdf", width = 4, height = 4)

for(i in 1:length(chrs)){

  chr <- chrs[i]
  print(chr)

  Dev <- boxplot.dat2[grep("Dev", boxplot.dat2$ID), chr]
  IR <- boxplot.dat2[grep("IR", boxplot.dat2$ID), chr]

  b <- cohen.d(Dev, #effect size
        IR,
        hedges.correction=TRUE
       )

       p <- ggviolin(boxplot.dat2,
              x = "ID",
              y = chr,
              color = "ID",
              palette = c("red", "black"),
              add = "mean_sd",
              main = paste0("Hedges g: ", round(b$estimate, 3), "\nMagnitude: ", b$magnitude),
              legend = "none"
              ) + stat_compare_means() #  Add p-value
              effect.size2[i] <- b$estimate
              #print(p)
}

#dev.off()
##############################################################
# 3) Plot Effect Size Across Chromosomes
##############################################################

pdf("HYbrid_CNV_Enrichment_EffectSize.pdf", height = 5, width = 8)

cols <-rep("grey", length(chrs))
cols[c(12,17,28)] <- "red"

names(effect.size) <- chrs
p <- barplot(effect.size,
        las = 2,
        ylim = c(-1,2),
        ylab = "Effect Size (Hedge's g)",
        col = cols
       )
abline(h=0.8, lty = 2, col = "black")
abline(h=-0.8,  lty = 2, col = "black")
abline(h=0)
points(x = p, y = effect.size2 ) #add points from experiment without hybrid cells

dev.off()


### plot out ones with large effect size
pdf("Hybrid_enrichedCNVs_largeEffectSize.pdf", width = 8, height = 3)
ggviolin(boxplot.dat,
         x = "ID",
         y = c("chr6q", "chr9p", "chr19p"),
         combine = TRUE,
         color = "ID",
         palette = c("red", "black"),
         ylab = "InferCNV Score \n(chr arm mean)",
         add = "mean_sd",
         legend = "none",
         xlab = ""
        )
dev.off()