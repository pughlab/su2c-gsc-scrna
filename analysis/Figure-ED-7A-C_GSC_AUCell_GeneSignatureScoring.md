
---
# Score GSCs with Developmental and Injury Response gene signatures with AUCell
---
L.Richards  
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/ScoreBTSCs

C1 = Developmental  
C2 = Injury Response  

----
# 1.0 Score with AUCell
----
#run on samwise screen -r 5943.pts-4.samwise

module load igenome-human/hg38
module load R/3.5.0

Rscript AUcell_BTSC.R > AUCell_log.txtlibrary("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
       )
library(AUCell)
library(RColorBrewer)
#library(NMF)

print("load gene signature for scoring")
print(Sys.time())

load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/AUCell/TumourOnly/StemSig/GeneSignatures.RData")
str(sigs)

print("load global BTSC object")
print(Sys.time())

load("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs.Rdata")

print("Isolate expression matrix'")
print(Sys.time())

exprMatrix <- as.matrix(BTSC@data)
exprMatrix[1:10, 1:10]
dim(exprMatrix)

print("AUCell: Build cell rankings")
print(Sys.time())

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)

print("AUCell: Calculate enrichment for gene signatures (AUC)")
print(Sys.time())

cells_AUC <- AUCell_calcAUC(sigs, cells_rankings)
save(cells_AUC, file="BTSCs_AUC_C1C2Stem.RData")

print("Output plots with AUC thresholds")
print(Sys.time())

set.seed(123)
par(mfrow=c(2,3)) 

pdf("BTSCs_Genesignatures_AUCell.pdf", width = 11, height = 9)

cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                             plotHist=TRUE, 
                                             assign=TRUE) 
#dev.off()


tsne_coords <- BTSC@dr$tsne@cell.embeddings
umap_coords <- BTSC@dr$umap@cell.embeddings
selectedThresholds <- getThresholdSelected(cells_assignment)

AUCell_plotTSNE(tSNE = tsne_coords, 
                exprMat = exprMatrix, 
                cellsAUC = cells_AUC, 
                thresholds = selectedThresholds
               )

AUCell_plotTSNE(tSNE = umap_coords, 
                exprMat = exprMatrix, 
                cellsAUC = cells_AUC, 
                thresholds = selectedThresholds
               )

dev.off()

AUC <- t(as.data.frame(cells_AUC@assays$data$AUC))
colnames(AUC) <- paste0(colnames(AUC), "_AUC")
head(AUC)

AUC <- data.frame(AUC)
head(AUC)

print("Save data")
print(Sys.time())

save(AUC, file = "BTSC_AUC_output.RData")

BTSC <- AddMetaData(BTSC,
                          metadata = AUC)

save(BTSC, file = "SU2C_BTSCs_Seurat_AUC.RData")
---
# Visualize Distribution of Scores across Samples
---

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/ScoreBTSCs/Cells


```R
library(ggplot2)
```

    Warning message:
    “package ‘ggplot2’ was built under R version 3.4.4”


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/ScoreBTSCs/Cells/BTSC_AUC_output.RData")
```

#### Plot distribution of AUCell scores



```R
C1 <- data.frame(AUC[ ,c("RNA.GSC.c1_AUC")])
rownames(C1) <- rownames(AUC)
C1$Sample <- gsub('.{17}$', '', rownames(C1))
colnames(C1) <- c("Score", "Sample")
C1$Cluster <- "C1"
C1$Sample <- with(C1, reorder(Sample, Score, median))

C2 <- data.frame(AUC[ ,c("RNA.GSC.c2_AUC")])
rownames(C2) <- rownames(AUC)
C2$Sample <- gsub('.{17}$', '', rownames(C2))
colnames(C2) <- c("Score", "Sample")
C2$Cluster <- "C2"
C2$Sample <- with(C1, reorder(Sample, Score, median))

dat <- rbind(C1, C2)
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>Score</th><th scope=col>Sample</th><th scope=col>Cluster</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>0.1845547</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>0.1326318</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>0.1362233</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>0.1508534</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>0.1330394</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>0.1601819</td><td>BT127_L  </td><td>C1       </td></tr>
</tbody>
</table>




```R
C1_z <- data.frame(AUC[ ,c("RNA.GSC.c1_AUC")])
rownames(C1_z) <- rownames(AUC)
C1_z$Sample <- gsub('.{17}$', '', rownames(C1_z))
colnames(C1_z) <- c("Score", "Sample")
C1_z$Cluster <- "C1"
C1_z$Score <- scale(C1_z$Score, center = TRUE, scale = TRUE)
C1_z$Sample <- with(C1_z, reorder(Sample, Score, median))

C2_z <- data.frame(AUC[ ,c("RNA.GSC.c2_AUC")])
rownames(C2_z) <- rownames(AUC)
C2_z$Sample <- gsub('.{17}$', '', rownames(C2_z))
colnames(C2_z) <- c("Score", "Sample")
C2_z$Cluster <- "C2"
C2_z$Score <- scale(C2_z$Score, center = TRUE, scale = TRUE)

C2_z$Sample <- with(C1_z, reorder(Sample, Score, median))

dat_z <- rbind(C1_z, C2_z)
head(dat_z)
```


<table>
<thead><tr><th></th><th scope=col>Score</th><th scope=col>Sample</th><th scope=col>Cluster</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>1.9046355</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>0.6596044</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>0.7457241</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>1.0965305</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>0.6693790</td><td>BT127_L  </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>1.3202153</td><td>BT127_L  </td><td>C1       </td></tr>
</tbody>
</table>




```R
pdf("~/Desktop/NewCohort_C1.C2_gradient.pdf", width = 10, height = 3)

p <- ggplot(dat, aes(x=Sample, y=Score, fill = Cluster)) +theme_classic() + 
      geom_violin()  + ylab("Raw AUC Score") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
      scale_fill_manual(values=c("red", "black")) + scale_x_discrete()
      
p


q <- ggplot(dat_z, aes(x=Sample, y=Score, fill = Cluster)) +theme_classic() + 
      geom_violin()  + ylab("AUC Z-score") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
      scale_fill_manual(values=c("red", "black")) + scale_x_discrete()
      
q

dev.off()
```






<strong>pdf:</strong> 2


---
# Visualize Classification of Cells Across Samples
---


```R
dat <- data.frame(AUC[ ,c("RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")])
colnames(dat) <- c("C1_AUC", "C2_AUC")

dat$Sample <- gsub('.{17}$', '', rownames(dat))


dat$C1 <- dat$C1_AUC > 0.11
dat$C2 <- dat$C2_AUC > 0.2
dat$ID <- paste0(dat$C1, dat$C2)

dat$ID <- gsub("TRUEFALSE", "C1", dat$ID)
dat$ID <- gsub("TRUETRUE", "C1+C2", dat$ID)
dat$ID <- gsub("FALSETRUE", "C2", dat$ID)
dat$ID <- gsub("FALSEFALSE", "Unclassified", dat$ID)

head(dat)
table(dat$ID)

a <- prop.table(table(dat$ID, dat$Sample),2)
head(a)
```


<table>
<thead><tr><th></th><th scope=col>C1_AUC</th><th scope=col>C2_AUC</th><th scope=col>Sample</th><th scope=col>C1</th><th scope=col>C2</th><th scope=col>ID</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>0.1845547</td><td>0.1712612</td><td>BT127_L  </td><td>TRUE     </td><td>FALSE    </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>0.1326318</td><td>0.1620893</td><td>BT127_L  </td><td>TRUE     </td><td>FALSE    </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>0.1362233</td><td>0.1541338</td><td>BT127_L  </td><td>TRUE     </td><td>FALSE    </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>0.1508534</td><td>0.1656422</td><td>BT127_L  </td><td>TRUE     </td><td>FALSE    </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>0.1330394</td><td>0.1608234</td><td>BT127_L  </td><td>TRUE     </td><td>FALSE    </td><td>C1       </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>0.1601819</td><td>0.2073333</td><td>BT127_L  </td><td>TRUE     </td><td> TRUE    </td><td>C1+C2    </td></tr>
</tbody>
</table>




    
              C1        C1+C2           C2 Unclassified 
           25336         2170        40254         1633 



                  
                        BT127_L      BT147_L       BT48_L       BT67_L       BT73_L
      C1           0.7825484765 0.9976798144 0.9979550102 0.9758735441 0.6707095209
      C1+C2        0.2029085873 0.0000000000 0.0013633265 0.0166389351 0.0084899939
      C2           0.0055401662 0.0000000000 0.0006816633 0.0024958403 0.0030321407
      Unclassified 0.0090027701 0.0023201856 0.0000000000 0.0049916805 0.3177683445
                  
                         BT84_L       BT89_L       BT94_L       G523_L       G549_L
      C1           0.9699331849 0.9966405375 0.9991823385 0.9540720961 0.0049504950
      C1+C2        0.0250556793 0.0000000000 0.0008176615 0.0072096128 0.0043680839
      C2           0.0005567929 0.0000000000 0.0000000000 0.0253671562 0.9868957484
      Unclassified 0.0044543430 0.0033594625 0.0000000000 0.0133511348 0.0037856727
                  
                         G564_L       G566_L       G583_L       G620_L       G637_L
      C1           0.0026975991 0.0000000000 0.0099942890 0.9921085859 0.9342030630
      C1+C2        0.0097113569 0.0004668534 0.0122786979 0.0066287879 0.0391378332
      C2           0.9867817642 0.9995331466 0.8520845231 0.0006313131 0.0175836642
      Unclassified 0.0008092797 0.0000000000 0.1256424900 0.0006313131 0.0090754396
                  
                         G729_L       G797_L       G799_L       G800_L       G837_L
      C1           0.0000000000 0.0051724138 0.0875706215 0.0125735688 0.0008534244
      C1+C2        0.0010413051 0.0336206897 0.1704331450 0.2851792402 0.0000000000
      C2           0.9989586949 0.9586206897 0.6911487759 0.7011771001 0.9790911031
      Unclassified 0.0000000000 0.0025862069 0.0508474576 0.0010700910 0.0200554726
                  
                         G851_L       G876_L       G885_L       G895_L     G945-I_L
      C1           0.0766155896 0.0827338129 0.9370229008 0.9281345566 0.0000000000
      C1+C2        0.0146568954 0.1187050360 0.0038167939 0.0275229358 0.0003824823
      C2           0.7768154564 0.7913669065 0.0028625954 0.0099388379 0.9992350354
      Unclassified 0.1319120586 0.0071942446 0.0562977099 0.0344036697 0.0003824823
                  
                       G945-J_L     G945-K_L     G946-J_L     G946-K_L
      C1           0.0000000000 0.0000000000 0.9006505027 1.0000000000
      C1+C2        0.0055887454 0.0005433306 0.0183323477 0.0000000000
      C2           0.9934476778 0.9956533551 0.0496747487 0.0000000000
      Unclassified 0.0009635768 0.0038033143 0.0313424009 0.0000000000



```R
a <- as.data.frame.matrix(a)
order <- c("C1", "C2", "C1+C2", "Unclassified")
col.order <- names(sort(a[1,]))
col.order
a <- as.matrix(a[order, col.order])
a
```


<ol class=list-inline>
	<li>'G566_L'</li>
	<li>'G729_L'</li>
	<li>'G945-I_L'</li>
	<li>'G945-J_L'</li>
	<li>'G945-K_L'</li>
	<li>'G837_L'</li>
	<li>'G564_L'</li>
	<li>'G549_L'</li>
	<li>'G797_L'</li>
	<li>'G583_L'</li>
	<li>'G800_L'</li>
	<li>'G851_L'</li>
	<li>'G876_L'</li>
	<li>'G799_L'</li>
	<li>'BT73_L'</li>
	<li>'BT127_L'</li>
	<li>'G946-J_L'</li>
	<li>'G895_L'</li>
	<li>'G637_L'</li>
	<li>'G885_L'</li>
	<li>'G523_L'</li>
	<li>'BT84_L'</li>
	<li>'BT67_L'</li>
	<li>'G620_L'</li>
	<li>'BT89_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
	<li>'BT94_L'</li>
	<li>'G946-K_L'</li>
</ol>




<table>
<thead><tr><th></th><th scope=col>G566_L</th><th scope=col>G729_L</th><th scope=col>G945-I_L</th><th scope=col>G945-J_L</th><th scope=col>G945-K_L</th><th scope=col>G837_L</th><th scope=col>G564_L</th><th scope=col>G549_L</th><th scope=col>G797_L</th><th scope=col>G583_L</th><th scope=col>⋯</th><th scope=col>G885_L</th><th scope=col>G523_L</th><th scope=col>BT84_L</th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>BT89_L</th><th scope=col>BT147_L</th><th scope=col>BT48_L</th><th scope=col>BT94_L</th><th scope=col>G946-K_L</th></tr></thead>
<tbody>
	<tr><th scope=row>C1</th><td>0.0000000000</td><td>0.000000000 </td><td>0.0000000000</td><td>0.0000000000</td><td>0.0000000000</td><td>0.0008534244</td><td>0.0026975991</td><td>0.004950495 </td><td>0.005172414 </td><td>0.009994289 </td><td>⋯           </td><td>0.937022901 </td><td>0.954072096 </td><td>0.9699331849</td><td>0.975873544 </td><td>0.9921085859</td><td>0.996640538 </td><td>0.997679814 </td><td>0.9979550102</td><td>0.9991823385</td><td>1           </td></tr>
	<tr><th scope=row>C2</th><td>0.9995331466</td><td>0.998958695 </td><td>0.9992350354</td><td>0.9934476778</td><td>0.9956533551</td><td>0.9790911031</td><td>0.9867817642</td><td>0.986895748 </td><td>0.958620690 </td><td>0.852084523 </td><td>⋯           </td><td>0.002862595 </td><td>0.025367156 </td><td>0.0005567929</td><td>0.002495840 </td><td>0.0006313131</td><td>0.000000000 </td><td>0.000000000 </td><td>0.0006816633</td><td>0.0000000000</td><td>0           </td></tr>
	<tr><th scope=row>C1+C2</th><td>0.0004668534</td><td>0.001041305 </td><td>0.0003824823</td><td>0.0055887454</td><td>0.0005433306</td><td>0.0000000000</td><td>0.0097113569</td><td>0.004368084 </td><td>0.033620690 </td><td>0.012278698 </td><td>⋯           </td><td>0.003816794 </td><td>0.007209613 </td><td>0.0250556793</td><td>0.016638935 </td><td>0.0066287879</td><td>0.000000000 </td><td>0.000000000 </td><td>0.0013633265</td><td>0.0008176615</td><td>0           </td></tr>
	<tr><th scope=row>Unclassified</th><td>0.0000000000</td><td>0.000000000 </td><td>0.0003824823</td><td>0.0009635768</td><td>0.0038033143</td><td>0.0200554726</td><td>0.0008092797</td><td>0.003785673 </td><td>0.002586207 </td><td>0.125642490 </td><td>⋯           </td><td>0.056297710 </td><td>0.013351135 </td><td>0.0044543430</td><td>0.004991681 </td><td>0.0006313131</td><td>0.003359462 </td><td>0.002320186 </td><td>0.0000000000</td><td>0.0000000000</td><td>0           </td></tr>
</tbody>
</table>




```R
#pdf("~/Desktop/NewCohort_PropC1.C2.Clasificaton_AUCell.pdf")

par(xpd = T, mar = par()$mar + c(0,0,7,0))
barplot(a,
        col = c("red", "black", "lightgrey", "white"),
        las = 2,
        ylab = "Proportion Cells"
       )

legend(0, 1.3, 
       legend = c("Proneural (C1)", "Immunomesenchymal (C2)", "Intermediate High State (C1+C2+)", "Intermediate Low State (C1-C2-)"), 
       fill = c("red", "black", "grey", "white"),
       bty = 'n'
      )

#dev.off()
```


![png](output_14_0.png)



```R
rbPal <- colorRampPalette(c("black", "red"))


#This adds a column of color values
# based on the y values
dat$Col<- rbPal(100)[as.numeric(cut(dat$C1_AUC, breaks = 100))]
unique(dat$Col)
```


<ol class=list-inline>
	<li>'#880000'</li>
	<li>'#570000'</li>
	<li>'#5A0000'</li>
	<li>'#690000'</li>
	<li>'#710000'</li>
	<li>'#6E0000'</li>
	<li>'#900000'</li>
	<li>'#4D0000'</li>
	<li>'#760000'</li>
	<li>'#7B0000'</li>
	<li>'#610000'</li>
	<li>'#7E0000'</li>
	<li>'#730000'</li>
	<li>'#520000'</li>
	<li>'#850000'</li>
	<li>'#B10000'</li>
	<li>'#670000'</li>
	<li>'#5F0000'</li>
	<li>'#8D0000'</li>
	<li>'#800000'</li>
	<li>'#640000'</li>
	<li>'#5C0000'</li>
	<li>'#4A0000'</li>
	<li>'#6C0000'</li>
	<li>'#830000'</li>
	<li>'#8B0000'</li>
	<li>'#920000'</li>
	<li>'#550000'</li>
	<li>'#790000'</li>
	<li>'#3D0000'</li>
	<li>'#380000'</li>
	<li>'#9A0000'</li>
	<li>'#AC0000'</li>
	<li>'#950000'</li>
	<li>'#A70000'</li>
	<li>'#970000'</li>
	<li>'#9D0000'</li>
	<li>'#AA0000'</li>
	<li>'#A20000'</li>
	<li>'#450000'</li>
	<li>'#420000'</li>
	<li>'#480000'</li>
	<li>'#9F0000'</li>
	<li>'#4F0000'</li>
	<li>'#400000'</li>
	<li>'#3B0000'</li>
	<li>'#B40000'</li>
	<li>'#A40000'</li>
	<li>'#C10000'</li>
	<li>'#AF0000'</li>
	<li>'#D00000'</li>
	<li>'#B60000'</li>
	<li>'#B90000'</li>
	<li>'#BC0000'</li>
	<li>'#C80000'</li>
	<li>'#C60000'</li>
	<li>'#BE0000'</li>
	<li>'#300000'</li>
	<li>'#210000'</li>
	<li>'#330000'</li>
	<li>'#360000'</li>
	<li>'#2E0000'</li>
	<li>'#2B0000'</li>
	<li>'#290000'</li>
	<li>'#D50000'</li>
	<li>'#E50000'</li>
	<li>'#C30000'</li>
	<li>'#D80000'</li>
	<li>'#CB0000'</li>
	<li>'#260000'</li>
	<li>'#240000'</li>
	<li>'#190000'</li>
	<li>'#1C0000'</li>
	<li>'#1E0000'</li>
	<li>'#170000'</li>
	<li>'#140000'</li>
	<li>'#120000'</li>
	<li>'#0F0000'</li>
	<li>'#0A0000'</li>
	<li>'#0C0000'</li>
	<li>'#050000'</li>
	<li>'#070000'</li>
	<li>'#000000'</li>
	<li>'#020000'</li>
	<li>'#D30000'</li>
	<li>'#CE0000'</li>
	<li>'#E20000'</li>
	<li>'#FF0000'</li>
</ol>




```R
samples <- unique(dat$Sample)

pdf("~/Desktop/NewCohortAUC_cell_sample.pdf")

for(i in 1:length(samples)){

sample <- samples[i]
print(sample)
sample.dat <- dat[dat$Sample == sample, ]
print(dim(sample.dat))



plot(sample.dat$C1_AUC, 
     sample.dat$C2_AUC,
     xlim = c(min(dat$C1_AUC), max(dat$C1_AUC)),
     ylim = c(min(dat$C2_AUC), max(dat$C2_AUC)),
     cex = 0.8,
     pch = 20,
     ylab = "Immunomesenchymal Program (C2)",
     xlab = "Proneural Program (C1)",
     col = sample.dat$Col
    )
abline(h=0.2, v = 0.11, lty = 2, col = "black", lwd = 2)
legend("topright",
       legend = paste0(sample, "\n", nrow(sample.dat), " cells" ),
      bty = 'n',
       xjust = 1
      )
    
   }

dev.off()
```

    [1] "BT127_L"
    [1] 1444    7
    [1] "BT147_L"
    [1] 862   7
    [1] "BT48_L"
    [1] 1467    7
    [1] "BT67_L"
    [1] 1202    7
    [1] "BT73_L"
    [1] 1649    7
    [1] "BT84_L"
    [1] 1796    7
    [1] "BT89_L"
    [1] 893   7
    [1] "BT94_L"
    [1] 1223    7
    [1] "G523_L"
    [1] 3745    7
    [1] "G549_L"
    [1] 3434    7
    [1] "G564_L"
    [1] 3707    7
    [1] "G566_L"
    [1] 2142    7
    [1] "G583_L"
    [1] 3502    7
    [1] "G620_L"
    [1] 3168    7
    [1] "G637_L"
    [1] 3526    7
    [1] "G729_L"
    [1] 2881    7
    [1] "G797_L"
    [1] 1160    7
    [1] "G799_L"
    [1] 1062    7
    [1] "G800_L"
    [1] 3738    7
    [1] "G837_L"
    [1] 4687    7
    [1] "G851_L"
    [1] 1501    7
    [1] "G876_L"
    [1] 834   7
    [1] "G885_L"
    [1] 1048    7
    [1] "G895_L"
    [1] 1308    7
    [1] "G945-I_L"
    [1] 5229    7
    [1] "G945-J_L"
    [1] 5189    7
    [1] "G945-K_L"
    [1] 3681    7
    [1] "G946-J_L"
    [1] 1691    7
    [1] "G946-K_L"
    [1] 1624    7



<strong>pdf:</strong> 2



```R

```
