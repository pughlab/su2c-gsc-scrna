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
### 2) Look for enriched CNVs
##############################################################

library(Seurat)
library(effsize) #v0.8.0
library(ggplot2)
library(reshape)
library(ggpubr)

setwd("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/EnrichmentAnalysis/effect_size")

##############################################################
# 1) Load data and classify tumour cells as IR or Dev
##############################################################
### Use AUCell cutoffs from GSCs
### Dev = 0.11
### IR = 0.2
dataset <- "Darmanis"
dat <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_public_scRNA/EnrichmentAnalysis/Darmanis_ChrArm_plotting.rds"

dat <- readRDS(dat)

dat$Dev <- dat$Developmental_GSC_AUC > 0.11 #AUCell cutoffs
dat$IR <- dat$InjuryResponse_GSC_AUC > 0.2 #AUCell cutodff
dat$ID <- paste0(dat$Dev, dat$IR)

dat$ID <- gsub("TRUEFALSE", "Dev", dat$ID)
dat$ID <- gsub("FALSETRUE", "IR", dat$ID)
dat$ID <- gsub("FALSEFALSE", "LowLow", dat$ID)
dat$ID <- gsub("TRUETRUE", "HiHi", dat$ID)

### get rid of hybrids
dat <- dat[!dat$ID == "HiHi", ]
dat <- dat[!dat$ID == "LowLow", ]

### Re-classify hybrid cells using highest score
### New cell numbers = 25,295 Dev + 40,360 IR
#for (i in 1:nrow(dat)){
#    cell <- dat$ID[i]
#    if(dat$ID[i] %in% c("HiHi" , "LowLow")){
#       dat$ID[i] <- ifelse(dat$C1_AUC[i] > dat$C2_AUC[i],
#                          print("Dev"),
#                          print("IR")
#                          )
#    }
#}


##############################################################
# 2) Compare CNV signal between bins
##############################################################

chrs <- gsub("_Avg", "", colnames(dat)[grep("Avg", colnames(dat))])
effect.size <- c()

plot.file <- paste0(dataset, "_chrArm_comparison.pdf")
pdf(plot.file, width = 4, height = 4)

for(i in 1:length(chrs)){

  chr <- chrs[i]
  print(chr)

  Dev <- dat[grep("Dev", dat$ID), paste0(chr, "_Avg")]
  IR <- dat[grep("IR", dat$ID), paste0(chr, "_Avg")]

  b <- cohen.d(Dev, #effect size
        IR,
        hedges.correction=TRUE
       )

       p <- ggviolin(dat,
              x = "ID",
              y = paste0(chr, "_Avg"),
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
# 2) Plot Effect Size Across Chromosomes
##############################################################

plot.file2 <- paste0(dataset, "_EffectSize.pdf")
pdf(plot.file2, height = 5, width = 8)

#cols <-rep("grey", length(chrs))
cols <- ifelse(abs(effect.size) > 0.8, "red", "grey")

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
#points(x = p, y = effect.size2 ) #add points from experiment without hybrid cells

dev.off()


### plot out ones with large effect size
plot.file3 <- paste0(dataset, "_GSC_chrArm_comparison.pdf")
pdf(plot.file3, width = 8, height = 3)
ggviolin(dat,
         x = "ID",
         y = c("chr6q_Avg", "chr9p_Avg", "chr19p_Avg"),
         combine = TRUE,
         color = "ID",
         palette = c("red", "black"),
         ylab = "InferCNV Score \n(chr arm mean)",
         add = "mean_sd",
         legend = "none",
         xlab = ""
        )
dev.off()
