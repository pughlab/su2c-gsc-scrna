
----
# Plotting malignant cell classification for GBM cells
---
L.Richards  
~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/plotting/

----
## 1.0 Quanitify malignant cell number across patients
----


```R
###load in malignant cell data frame
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/plotting/SU2C_TumourCellsOnly_tSNE_UMAP.Rdata")
```


```R
dimReduc.dat$PatientID <- paste0(dimReduc.dat$PatientID, "_T")
```


```R
patient.cols <- c("#1b9e77", 
                  "#d95f02",
                  "#7570b3",
                  "#e7298a",
                  "#66a61e",
                  "#e6ab02",
                  "#666666")

pdf("~/Desktop/TumourCellCounts.pdf", height = 5, width = 3)

xx <- barplot(table(dimReduc.dat$PatientID), 
        las = 2,
        col = patient.cols,
        ylim = c(0, 8000),
        #horiz = T,
              ylab = "Cells"
       )
text(x = xx, 
     y = table(dimReduc.dat$PatientID), 
               label = table(dimReduc.dat$PatientID), 
                             pos = 3, 
                             cex = 0.8, 
                             col = "black"
              )

dev.off()
```


<strong>pdf:</strong> 2


---
## 2.0 Set up plotting object
---
library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
       )

##load in all tumour data

load("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/plotting/SU2C_GBM_AllCells_Seurat.RData")

##marker gene list

markers <- c("PTPRC", #immune cells
             "MOG", #normal brain cells
             "MAG", #normal brain marker
             "EGFR", #tumour cell markers
             "CD3D", 
             "CD2",
             "ITGAM",
             "FCGR3A",
             "CD14",
             "TMEM119"
            )

##add UMAP coords to the metadata

merged <- AddMetaData(merged,
                      metadata = data.frame(merged@dr$umap@cell.embeddings)
                     
                     )

##add genes to meta.data
subset <- merged@data[markers, ]
subset <- t(data.frame(as.matrix(subset)))
rownames(subset) <- gsub("\\.", "-", rownames(subset))

merged <- AddMetaData(merged,
                      metadata = subset,
                      col.name = colnames(subset)
                     )

#### Add brain and immune scores to data

suppressPackageStartupMessages({
    #library(Seurat)
    library(AUCell)
    library(RColorBrewer)
    library(NMF)
})

exprMatrix <- as.matrix(merged@data)
dat[1:10, 1:10]
dim(exprMatrix)

load("~/pughlab/projects/Multiregional_GBM_scRNAseq/Manuscript/cell_annotation/GeneSignatures.RData")
gene.sets <- gene.signatures
str(gene.sets)

#build cell rankings
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file="SU2C_GBM_cells_rankings.RData")

#calculate enrichment for gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(gene.sets, cells_rankings)
save(cells_AUC, file="SU2C_GBM_cells_AUC.RData")

set.seed(123)
par(mfrow=c(2,3)) 

pdf("SU2C_GBM_signature_AUCell.pdf", width = 11, height = 9)

cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                             plotHist=TRUE, 
                                             assign=TRUE) 
#dev.off()

tsne_coords <- merged@dr$tsne@cell.embeddings
umap_coords <- merged@dr$umap@cell.embeddings
selectedThresholds <- getThresholdSelected(cells_assignment)


par(mfrow=c(2,3))

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

###add AUC data to metadat

merged <- AddMetaData(merged, metadata = AUC, col.name = colnames(AUC))

##add classification data

brain.cells <- cells_assignment$Brain.GTEx$assignment
immune.cells <- cells_assignment$ESTIMATE.Immune$assignment

merged@meta.data$BrainCells <- rownames(merged@meta.data) %in% brain.cells
merged@meta.data$ImmuneCells <- rownames(merged@meta.data) %in% immune.cells

merged@meta.data$Classification <- paste(merged@meta.data$BrainCells, merged@meta.data$ImmuneCells, sep = "_")
merged@meta.data$Classification <- gsub("FALSE_TRUE", "Immune", merged@meta.data$Classification)
merged@meta.data$Classification <- gsub("TRUE_TRUE", "Unclassified", merged@meta.data$Classification)
merged@meta.data$Classification <- gsub("FALSE_FALSE", "Unclassified", merged@meta.data$Classification)
merged@meta.data$Classification <- gsub("TRUE_FALSE", "Brain", merged@meta.data$Classification)
table(merged@meta.data$Classification)

##add in cutoffs for classification histogram
merged@meta.data$BrainCutoff <- 0.389
merged@meta.data$ImmuneCutoff <- 0.0692## save seruat object

meta <- data.frame(merged@meta.data)
save(meta, file = "GBM_AllCells_annotation_meta_Oct2019.Rdata")

save(merged, file = "GBM_AllCells_annotation_meta_Seurat_Oct2019.Rdata")###load in CONICS results and merge 
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/CNVs/CONICS/SU2C_GBM_CONICS_Brain.Rdata")

##add CONICS scores to metadata
cnv.meta <- CONICS_UMAP[ ,15:ncol(CONICS_UMAP)]
save(cnv.meta, file = "CNV_meta_GBM_CONICS.Rdata")
---
# 3.0 Plot Marker genes on UMAP
---



```R
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_annotation/plotting/")

```


```R
load("CNV_meta_GBM_CONICS.Rdata")
load("GBM_AllCells_annotation_meta_Oct2019.Rdata")

head(meta)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>SampleID</th><th scope=col>PatientID</th><th scope=col>SampleType</th><th scope=col>Pathology</th><th scope=col>Stage</th><th scope=col>percent.mito</th><th scope=col>res.1.5</th><th scope=col>⋯</th><th scope=col>FCGR3A</th><th scope=col>CD14</th><th scope=col>TMEM119</th><th scope=col>Brain.GTEx_AUC</th><th scope=col>ESTIMATE.Immune_AUC</th><th scope=col>BrainCells</th><th scope=col>ImmuneCells</th><th scope=col>Classification</th><th scope=col>BrainCutoff</th><th scope=col>ImmuneCutoff</th></tr></thead>
<tbody>
	<tr><th scope=row>G1003-A_T_AAACCTGAGTATCGAA</th><td>2356                   </td><td> 6656                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.03202526             </td><td>6                      </td><td>⋯                      </td><td>0                      </td><td>0.0000000              </td><td>0                      </td><td>0.3903429              </td><td>0.007161754            </td><td> TRUE                  </td><td>FALSE                  </td><td>Brain                  </td><td>0.389                  </td><td>0.0692                 </td></tr>
	<tr><th scope=row>G1003-A_T_AAACCTGCAGTCCTTC</th><td>3621                   </td><td>14575                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.04136036             </td><td>11                     </td><td>⋯                      </td><td>0                      </td><td>0.5229809              </td><td>0                      </td><td>0.4086430              </td><td>0.002351621            </td><td> TRUE                  </td><td>FALSE                  </td><td>Brain                  </td><td>0.389                  </td><td>0.0692                 </td></tr>
	<tr><th scope=row>G1003-A_T_AAACCTGTCAGTTTGG</th><td>3906                   </td><td>13102                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.09149278             </td><td>26                     </td><td>⋯                      </td><td>0                      </td><td>0.0000000              </td><td>0                      </td><td>0.3711021              </td><td>0.008631517            </td><td>FALSE                  </td><td>FALSE                  </td><td>Unclassified           </td><td>0.389                  </td><td>0.0692                 </td></tr>
	<tr><th scope=row>G1003-A_T_AAACCTGTCGGAATCT</th><td>3106                   </td><td> 8885                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.10908065             </td><td>7                      </td><td>⋯                      </td><td>0                      </td><td>0.0000000              </td><td>0                      </td><td>0.4240318              </td><td>0.006867802            </td><td> TRUE                  </td><td>FALSE                  </td><td>Brain                  </td><td>0.389                  </td><td>0.0692                 </td></tr>
	<tr><th scope=row>G1003-A_T_AAACCTGTCGTTTAGG</th><td>3217                   </td><td> 9028                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.08781083             </td><td>7                      </td><td>⋯                      </td><td>0                      </td><td>0.0000000              </td><td>0                      </td><td>0.4354693              </td><td>0.013976110            </td><td> TRUE                  </td><td>FALSE                  </td><td>Brain                  </td><td>0.389                  </td><td>0.0692                 </td></tr>
	<tr><th scope=row>G1003-A_T_AAACGGGAGAAACCGC</th><td>2502                   </td><td> 5496                  </td><td>G1003-A_T              </td><td>G1003-A_T              </td><td>G1003                  </td><td>Tumour                 </td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY                </td><td>0.06689756             </td><td>7                      </td><td>⋯                      </td><td>0                      </td><td>0.0000000              </td><td>0                      </td><td>0.3989979              </td><td>0.009638082            </td><td> TRUE                  </td><td>FALSE                  </td><td>Brain                  </td><td>0.389                  </td><td>0.0692                 </td></tr>
</tbody>
</table>




```R
patient.cols <- c("#1b9e77", 
                  "#d95f02",
                  "#7570b3",
                  "#e7298a",
                  "#66a61e",
                  "#e6ab02",
                  "#666666"
                 )



markers <- c("PTPRC", #immune cells
             "MOG", #normal brain cells
             "MAG", #normal brain marker
             "EGFR", #tumour cell markers
             "CD3D", 
             "CD2",
             "ITGAM",
             "FCGR3A",
             "CD14",
             "TMEM119"
            )

plist <- list()

for(i in 1:length(markers)){
    
gene <- markers[i]
print(gene)

p <- ggplot(meta, aes_string(x="UMAP1", y="UMAP2", color=gene)) + 
                   geom_point(alpha = 0.8, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                  scale_color_gradient(low="grey", high="blue") +  
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        labs(color = "Expression") + ggtitle(gene) +
                        theme(legend.position = "none")
                    
plist[[gene]] <- p
    
}
                    



```

    [1] "PTPRC"
    [1] "MOG"
    [1] "MAG"
    [1] "EGFR"
    [1] "CD3D"
    [1] "CD2"
    [1] "ITGAM"
    [1] "FCGR3A"
    [1] "CD14"
    [1] "TMEM119"



```R
pdf("~/Desktop/GeneMarkers_nolegend.pdf", width = 15, height = 6)

do.call("grid.arrange", c(plist, ncol=5))

dev.off()
```


<strong>pdf:</strong> 2


---
# 4.0 Plot Brain and Immune Gene Signatures From AUCelll
---


```R



        

brain_auc <- ggplot(meta, aes_string(x="UMAP1", y="UMAP2", color="Brain.GTEx_AUC")) + 
                   geom_point(alpha = 0.8, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                  scale_color_gradient(low="grey", high="blue") +  
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        labs(color = "AUC") + ggtitle("GTEx Brain") +
                        theme(legend.position = "none")



immune_auc <- ggplot(meta, aes_string(x="UMAP1", y="UMAP2", color="ESTIMATE.Immune_AUC")) + 
                   geom_point(alpha = 0.8, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                  scale_color_gradient(low="grey", high="blue") +  
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        labs(color = "AUC") + ggtitle("ESTIMATE Immune") +
                        theme(legend.position = "none")

brain_class <- ggplot(meta, aes_string(x="UMAP1", y="UMAP2", color="BrainCells")) + 
                   geom_point(alpha = 0.8, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                  scale_color_manual(values = c("grey", "black")) +  
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        labs(color = "Classification") + ggtitle("Brain") +  theme(legend.position='none')

immune_class <- ggplot(meta, aes_string(x="UMAP1", y="UMAP2", color="ImmuneCells")) + 
                   geom_point(alpha = 0.8, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                  scale_color_manual(values = c("grey", "black")) +  
                    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) + 
                        labs(color = "Classification") + ggtitle("Immune") + theme(legend.position='none')

```


```R
brain_hist <- qplot(meta$Brain.GTEx_AUC,
      geom = "histogram",
      bins = 60,
      xlab = "AUC",
      main = "GTEx Brain Gene Signature",
      fill = I("white"),
      col = I("black"),
      ylab = "Cell Frequency"
      ) + theme_classic() + geom_vline(xintercept = meta$BrainCutoff, col = "red")

immune_hist <- qplot(meta$ESTIMATE.Immune_AUC,
      geom = "histogram",
      bins = 60,
      xlab = "AUC",
      main = "ESTIMATE Immune Gene Signature",
      fill = I("white"),
      col = I("black"),
      ylab = "Cell Frequency"
      ) + theme_classic() + geom_vline(xintercept = meta$ImmuneCutoff, col = "red")


```


```R
pdf("~/Desktop/scRNA_GBM_AUCellClass_legends.pdf", width = 9, height = 6 )

grid.arrange(brain_hist, brain_auc, brain_class,
             immune_hist, immune_auc, immune_class,
             ncol = 3
            )
dev.off()


```


<strong>pdf:</strong> 2



```R

pdf("~/Desktop/scRNA_GBM_AUCellClass_NO_legends.pdf", width = 9, height = 6 )

grid.arrange(brain_hist, brain_auc, brain_class,
             immune_hist, immune_auc, immune_class,
             ncol = 3
            )
dev.off()
```


<strong>pdf:</strong> 2


---
# 5.0 Plot CONICS scores
---


```R
umaps <- meta[rownames(cnv.meta), c("UMAP1", "UMAP2") ]
cnv.meta <- cbind(cnv.meta, umaps)
head(cnv.meta)
```


<table>
<thead><tr><th></th><th scope=col>chr1p</th><th scope=col>chr1q</th><th scope=col>chr2p</th><th scope=col>chr2q</th><th scope=col>chr3p</th><th scope=col>chr3q</th><th scope=col>chr4p</th><th scope=col>chr4q</th><th scope=col>chr5p</th><th scope=col>chr5q</th><th scope=col>⋯</th><th scope=col>chr18q</th><th scope=col>chr19p</th><th scope=col>chr19q</th><th scope=col>chr20p</th><th scope=col>chr20q</th><th scope=col>chr21q</th><th scope=col>chr22q</th><th scope=col>hi</th><th scope=col>UMAP1</th><th scope=col>UMAP2</th></tr></thead>
<tbody>
	<tr><th scope=row>G1003-A_T_AACTCAGGTTTACTCT</th><td>0.006577591 </td><td>0.16988500  </td><td>0.92337107  </td><td>0.12749888  </td><td>0.003298319 </td><td>0.02425839  </td><td>0.3224717   </td><td>0.05034575  </td><td>0.2558471   </td><td>0.02837204  </td><td>⋯           </td><td>0.4444247   </td><td>0.9999769   </td><td>0.9165722   </td><td>0.002879431 </td><td>1.471056e-02</td><td>0.8640733   </td><td>0.9998397   </td><td>Normal      </td><td>7.978218    </td><td>-6.295468   </td></tr>
	<tr><th scope=row>G1003-A_T_AACTTTCCAGCTGTTA</th><td>0.005797757 </td><td>0.02621779  </td><td>0.55643371  </td><td>0.35206516  </td><td>0.006746463 </td><td>0.02650921  </td><td>0.4385753   </td><td>0.05156368  </td><td>0.1489874   </td><td>0.02853041  </td><td>⋯           </td><td>0.2571512   </td><td>0.9997533   </td><td>0.9253219   </td><td>0.018380466 </td><td>8.036185e-04</td><td>0.8133084   </td><td>0.9998852   </td><td>Normal      </td><td>9.953024    </td><td>-7.508868   </td></tr>
	<tr><th scope=row>G1003-A_T_ACGGAGAAGACAGACC</th><td>0.012871498 </td><td>0.03455282  </td><td>0.94086827  </td><td>0.08988776  </td><td>0.002682724 </td><td>0.03659838  </td><td>0.3468459   </td><td>0.06348949  </td><td>0.1547147   </td><td>0.04240630  </td><td>⋯           </td><td>0.3954680   </td><td>0.9999975   </td><td>0.8925683   </td><td>0.238281615 </td><td>1.323097e-06</td><td>0.7735523   </td><td>0.9928854   </td><td>Normal      </td><td>6.945198    </td><td>-8.325158   </td></tr>
	<tr><th scope=row>G1003-A_T_ACTATCTTCCAAAGTC</th><td>0.006254103 </td><td>0.02336432  </td><td>0.12071298  </td><td>0.09246091  </td><td>0.024273248 </td><td>0.02442434  </td><td>0.3788937   </td><td>0.05693365  </td><td>0.2509516   </td><td>0.02776765  </td><td>⋯           </td><td>0.3664499   </td><td>0.9998572   </td><td>0.8046951   </td><td>0.023671348 </td><td>1.076633e-03</td><td>0.8258113   </td><td>0.9937123   </td><td>Normal      </td><td>6.356769    </td><td>-7.776512   </td></tr>
	<tr><th scope=row>G1003-A_T_ACTGAGTCAGATGGCA</th><td>0.007125522 </td><td>0.02258658  </td><td>0.04465254  </td><td>0.09551925  </td><td>0.002704177 </td><td>0.02488746  </td><td>0.6710757   </td><td>0.12031854  </td><td>0.1939507   </td><td>0.02993228  </td><td>⋯           </td><td>0.9354430   </td><td>0.9999849   </td><td>0.9063837   </td><td>0.288770254 </td><td>4.874944e-03</td><td>0.8557956   </td><td>0.9991503   </td><td>Normal      </td><td>6.851467    </td><td>-8.251102   </td></tr>
	<tr><th scope=row>G1003-A_T_AGGCCACCATTTGCTT</th><td>0.006274028 </td><td>0.02243835  </td><td>0.57148966  </td><td>0.11547646  </td><td>0.002571897 </td><td>0.02539707  </td><td>0.3806867   </td><td>0.10870137  </td><td>0.1514356   </td><td>0.02932271  </td><td>⋯           </td><td>0.8248902   </td><td>0.9947717   </td><td>0.9259768   </td><td>0.417673785 </td><td>1.953605e-02</td><td>0.8565818   </td><td>0.9997476   </td><td>Normal      </td><td>6.983804    </td><td>-8.346859   </td></tr>
</tbody>
</table>




```R
regions <- c("chr7p", "chr7q", "chr10p", "chr10q")

cnv.list <- list()

for(i in 1:length(regions)){
    
    region <- regions[i]
    print(region)
    mid <- mean(cnv.meta[,region])

cnv_UMAP <- ggplot(cnv.meta, aes_string(x="UMAP1", y="UMAP2", color=region)) + 
                         geom_point(alpha = 0.6, size = 0.8, pch = 16) +  
                       labs(x = "UMAP 1", y = "UMAP 2") + 
                scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                     high="red", space ="Lab", name = "Posterior \nProbability \n Z-score") + 
                ggtitle(region) + theme_bw() + 
                 theme(legend.title=element_text(size=8),
                       axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank()) 
cnv.list[[region]]   <- cnv_UMAP
    }

cnv.list[[1]]
```

    [1] "chr7p"
    [1] "chr7q"
    [1] "chr10p"
    [1] "chr10q"





![png](output_22_2.png)



```R
pdf("~/Desktop/SU2C_GBM_CONICS_nolegend.pdf", width = 8, height = 8)

do.call("grid.arrange", c(cnv.list, ncol=2))

dev.off()
```


<strong>pdf:</strong> 2

