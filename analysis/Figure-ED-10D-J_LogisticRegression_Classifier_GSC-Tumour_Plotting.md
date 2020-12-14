
----
# Plot GSC-like tumour cells from classifier
---

L.Richards

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen

File from Owen:  
/cluster/projects/su2c_csc/scRNAseq_analysis/scRNA_seq_GBM_tissue/owen_scRNA_GSCs_Tumour_combined/logistic_regression_T_v_L_seurat.rds


```R
library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
       )

```

---
## 1.0 Format the data for plotting
---
dat <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_seurat.rds")
str(dat)

pca <- data.frame(dat@dr$pca@cell.embeddings)
pca <- pca[, 1:5]

meta <- data.frame(dat@meta.data)
meta <- cbind(meta, pca)
head(meta)

saveRDS(meta, file = "logistic_regression_T_v_L_metadata.rds")
---
## 2.0 Differential gene expression on cell types
---


```R
library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
       )

dat <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_seurat.rds")


#give cells new ID based on results of classifier

dat@meta.data$Class <- paste0(dat@meta.data$pred_class, dat@meta.data$misclassified)
table(dat@meta.data$Class)


dat@meta.data$Class <- gsub("TumourFALSE", "Tumour", dat@meta.data$Class)
dat@meta.data$Class <- gsub("LineFALSE", "GSC", dat@meta.data$Class)
dat@meta.data$Class <- gsub("LineTRUE", "GSC-like_Tumour", dat@meta.data$Class)
dat@meta.data$Class <- gsub("TumourTRUE", "Tumour-like_GSC", dat@meta.data$Class)
table(dat@meta.data$Class)

dat <- SetIdent(dat, ident.use = factor(dat@meta.data$Class))
table(dat@ident)

#saveRDS(dat, file = "logistic_regression_T_v_L_seurat_classification.rds")



## Do DE in Seruat with Wilcoxon rank sum test

### 1. GSC-like tumour cells versus tumour "GSC.like.T_markers"
### 2. Tumour-like GSCs versus tumour
### 3.0 All GSCs versus all tumour (not based on linear regression results methods)

GSC.like.T_markers <- FindMarkers(dat, 
                                  ident.1 = "GSC-like_Tumour", 
                                  ident.2 = "Tumour", 
                                  test.use = "wilcox"
                                 )
saveRDS(GSC.like.T_markers, 
        file = "GSC.like.T_markers.rds"
       )

write.csv(GSC.like.T_markers, 
          file = "GSC.like.T_markers.csv",
          quote = FALSE
          )

Tumour.like.GSCs_markers <- FindMarkers(dat, 
                                        ident.1 = "Tumour-like_GSC", 
                                        ident.2 = "GSC", 
                                        test.use = "wilcox"
                                       )

saveRDS(Tumour.like.GSCs_markers, 
        file = "Tumour-like_GSC_markers.rds"
       )

write.csv(Tumour.like.GSCs_markers, 
          file = "Tumour-like_GSC_markers.csv",
          quote = FALSE
          )
          
          
#### Add marker genes from DE to meta.data for downstream plotting

genes <- c('SUSD2', 
           'SAA1',
           'CAPS',
           'CHI3L1',
           'AQP4',
           'GFAP',
           'ZFP36',
           'TPPP3',
           'NTRK2',
           'ID3',
           'HBB',
           'HILPDA',
           'CYR61',
           'CLU',
           'C9orf24',
           'HMGB2', 
           'TYMS', 
           'STMN1', 
           'MLLT11', 
           'VGF', 
           'MDM2', 
           'CCT6A', 
           'DLL3', 
           'STMN2' ,
           'TUBB' ,
           'CDK4', 
           'NUP107' ,
           'SOX11' ,
           'SOX4' ,
           'HIST1H4C', 
           'SEC61G',
           'CDK4', 
           'CDKN2A', 
           'HIST1H4C' ,
           'OCIAD2' ,
           'MDM2' ,
           'H2AFZ', 
           'NEFL' ,
           'PTMA', 
           'HMGN2' ,
           'SET' ,
           'DEK' ,
           'MDK' ,
           'DKK1' ,
           'PPP1R14B', 
           'HMGN1',
           'TTYH1', 
           'HOPX' ,
           'COL1A2', 
           'CHI3L1', 
           'GFAP' ,
           'CRYAB' ,
           'NOVA1', 
           'WSB1', 
           'NMB', 
           'RIC3' ,
           'AGT' ,
           'MEG3', 
           'CLU', 
           'APOE',
           'NEAT1' ,
           'MALAT1'
           )
genes <- as.character(unique(genes))
genes

subset <- t(data.matrix(dat@data[genes, ]))
head(subset)

pc <- dat@dr$pca@cell.embeddings

meta <- dat@meta.data
head(meta)

meta.genes <- cbind(meta, subset, pc)
head(meta.genes)

saveRDS(meta.genes, file = "logistic_regression_T_v_L_metadata_DEgenes.rds")
saveRDS(dat, file = "logistic_regression_T_v_L_seurat_classification.rds")

```

---
## 3.0 Visualize DE Genes of GSC-like tumour cells
---



```R
GSC.like.tumour <- read.csv("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/GSC.like.T_markers.csv")

GSC.like.tumour$log2FC <- log2(exp(GSC.like.tumour$avg_logFC))
GSC.like.tumour$log10.adj.p <-  -log10(GSC.like.tumour$p_val_adj)

#filter by FDR of 0.01
GSC.like.tumour <- GSC.like.tumour[GSC.like.tumour$p_val_adj <= 0.01, ]
dim(GSC.like.tumour)

head(GSC.like.tumour)
#plot(GSC.like.tumour$p_val_adj)
```


<ol class=list-inline>
	<li>896</li>
	<li>8</li>
</ol>




<table>
<thead><tr><th scope=col>X</th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FC</th><th scope=col>log10.adj.p</th></tr></thead>
<tbody>
	<tr><td>TUBB     </td><td>0        </td><td> 1.045300</td><td>0.887    </td><td>0.757    </td><td>0        </td><td> 1.508048</td><td>Inf      </td></tr>
	<tr><td>GJA1     </td><td>0        </td><td>-1.114220</td><td>0.110    </td><td>0.600    </td><td>0        </td><td>-1.607480</td><td>Inf      </td></tr>
	<tr><td>ANXA1    </td><td>0        </td><td>-1.129634</td><td>0.342    </td><td>0.788    </td><td>0        </td><td>-1.629717</td><td>Inf      </td></tr>
	<tr><td>CLU      </td><td>0        </td><td>-1.179591</td><td>0.840    </td><td>0.960    </td><td>0        </td><td>-1.701791</td><td>Inf      </td></tr>
	<tr><td>NTRK2    </td><td>0        </td><td>-1.257728</td><td>0.213    </td><td>0.725    </td><td>0        </td><td>-1.814519</td><td>Inf      </td></tr>
	<tr><td>GFAP     </td><td>0        </td><td>-1.370528</td><td>0.581    </td><td>0.862    </td><td>0        </td><td>-1.977254</td><td>Inf      </td></tr>
</tbody>
</table>




```R
test <- GSC.like.tumour[order(GSC.like.tumour$log2FC), ]
test$Index <- 1:nrow(test)

top15 <- as.character(test[1:15, ]$X)
bottom15 <- as.character(test[(nrow(test)-15):nrow(test), ]$X)

top15
bottom15

GSC.T <- ggscatter(test, 
          x = "Index", 
          y = "log2FC",
          color = "grey", 
          shape = 20,
          palette = c("grey"),
          label = "X", 
          label.select = c(top15, bottom15),
          #label.rectangle = TRUE,
          font.label = c(7, "black", "plain"),
          repel = "TRUE",
          size = 1,
          xlab = "Ranked Genes",
          ylab = "Average log2 Fold Change"
         )  +  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      geom_hline(yintercept = c(1,-1), lty = 2, col = "darkred") + ggtitle("DE GSC-like tumour vs tumour (FDR < 0.01)")
```


<ol class=list-inline>
	<li>'SUSD2'</li>
	<li>'SAA1'</li>
	<li>'CAPS'</li>
	<li>'CHI3L1'</li>
	<li>'AQP4'</li>
	<li>'GFAP'</li>
	<li>'ZFP36'</li>
	<li>'TPPP3'</li>
	<li>'NTRK2'</li>
	<li>'ID3'</li>
	<li>'HBB'</li>
	<li>'HILPDA'</li>
	<li>'CYR61'</li>
	<li>'CLU'</li>
	<li>'C9orf24'</li>
</ol>




<ol class=list-inline>
	<li>'HMGB2'</li>
	<li>'TYMS'</li>
	<li>'STMN1'</li>
	<li>'MLLT11'</li>
	<li>'VGF'</li>
	<li>'MDM2'</li>
	<li>'CCT6A'</li>
	<li>'DLL3'</li>
	<li>'STMN2'</li>
	<li>'TUBB'</li>
	<li>'CDK4'</li>
	<li>'NUP107'</li>
	<li>'SOX11'</li>
	<li>'SOX4'</li>
	<li>'HIST1H4C'</li>
	<li>'SEC61G'</li>
</ol>



---
## 4.0 Visualize DE Genes of tumour-like GSC cells
---


```R
Tumour.like.GSC <- read.csv("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/Tumour-like_GSC_markers.csv")

Tumour.like.GSC$log2FC <- log2(exp(Tumour.like.GSC$avg_logFC))
Tumour.like.GSC$log10.adj.p <-  -log10(Tumour.like.GSC$p_val_adj)

#filter by FDR of 0.01
Tumour.like.GSC <- Tumour.like.GSC[Tumour.like.GSC$p_val_adj <= 0.01, ]
dim(Tumour.like.GSC)

head(Tumour.like.GSC)
#plot(GSC.like.tumour$p_val_adj)
```


<ol class=list-inline>
	<li>1125</li>
	<li>8</li>
</ol>




<table>
<thead><tr><th scope=col>X</th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FC</th><th scope=col>log10.adj.p</th></tr></thead>
<tbody>
	<tr><td>AQP4      </td><td>0         </td><td> 0.3472748</td><td>0.214     </td><td>0.025     </td><td>0         </td><td> 0.5010117</td><td>Inf       </td></tr>
	<tr><td>TMA7      </td><td>0         </td><td>-0.7491595</td><td>0.694     </td><td>0.984     </td><td>0         </td><td>-1.0808087</td><td>Inf       </td></tr>
	<tr><td>PPIA      </td><td>0         </td><td>-0.7592584</td><td>0.797     </td><td>0.996     </td><td>0         </td><td>-1.0953784</td><td>Inf       </td></tr>
	<tr><td>SRP14     </td><td>0         </td><td>-0.7597104</td><td>0.691     </td><td>0.986     </td><td>0         </td><td>-1.0960304</td><td>Inf       </td></tr>
	<tr><td>SNRPG     </td><td>0         </td><td>-0.8033821</td><td>0.624     </td><td>0.967     </td><td>0         </td><td>-1.1590354</td><td>Inf       </td></tr>
	<tr><td>MRPL33    </td><td>0         </td><td>-0.8953037</td><td>0.566     </td><td>0.932     </td><td>0         </td><td>-1.2916502</td><td>Inf       </td></tr>
</tbody>
</table>




```R
test <- Tumour.like.GSC[order(Tumour.like.GSC$log2FC), ]
test$Index <- 1:nrow(test)

top15 <- as.character(test[1:15, ]$X)
bottom15 <- as.character(test[(nrow(test)-15):nrow(test), ]$X)

top15
bottom15

Tumour.GSC <- ggscatter(test, 
          x = "Index", 
          y = "log2FC",
          color = "grey", 
          shape = 20,
          palette = c("grey"),
          label = "X", 
          label.select = c(top15, bottom15),
          #label.rectangle = TRUE,
          font.label = c(7, "black", "plain"),
          repel = "TRUE",
          size = 1,
          xlab = "Ranked Genes",
          ylab = "Average log2 Fold Change"
         )  +  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      geom_hline(yintercept = c(1,-1), lty = 2, col = "darkred") + ggtitle("DE Tumour-like GSC vs GSC (FDR < 0.01)")
```


<ol class=list-inline>
	<li>'CDK4'</li>
	<li>'CDKN2A'</li>
	<li>'HIST1H4C'</li>
	<li>'OCIAD2'</li>
	<li>'MDM2'</li>
	<li>'H2AFZ'</li>
	<li>'NEFL'</li>
	<li>'PTMA'</li>
	<li>'HMGN2'</li>
	<li>'SET'</li>
	<li>'DEK'</li>
	<li>'MDK'</li>
	<li>'DKK1'</li>
	<li>'PPP1R14B'</li>
	<li>'HMGN1'</li>
</ol>




<ol class=list-inline>
	<li>'TTYH1'</li>
	<li>'HOPX'</li>
	<li>'COL1A2'</li>
	<li>'CHI3L1'</li>
	<li>'GFAP'</li>
	<li>'CRYAB'</li>
	<li>'NOVA1'</li>
	<li>'WSB1'</li>
	<li>'NMB'</li>
	<li>'RIC3'</li>
	<li>'AGT'</li>
	<li>'MEG3'</li>
	<li>'CLU'</li>
	<li>'APOE'</li>
	<li>'NEAT1'</li>
	<li>'MALAT1'</li>
</ol>



-----


```R
pdf("~/Desktop/DE_Classifier.pdf", height = 6, width =12)
ggarrange(GSC.T, Tumour.GSC, ncol = 2)
dev.off()
```




<strong>pdf:</strong> 2


---
## 4.0 Compare expression of cells across classification categories
---


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/")

dat <- readRDS("logistic_regression_T_v_L_metadata_DEgenes.rds")
head(dat)
colnames(dat)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640        </td><td>  875       </td><td>BTSC        </td><td>0.043428571 </td><td> 0.0957859  </td><td> 0.07303696 </td><td>S           </td><td> 0.022748941</td><td>0.10663091  </td><td>0.000000000 </td><td>⋯           </td><td>1.3949986   </td><td>-1.2290755  </td><td> 0.04078296 </td><td>-0.3940755  </td><td>-4.3912124  </td><td>-0.5042649  </td><td> 3.3774399  </td><td> 7.9740340  </td><td>-0.9635511  </td><td>-3.9862054  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036        </td><td> 2408       </td><td>BTSC        </td><td>0.002076412 </td><td> 0.0535880  </td><td> 0.30872825 </td><td>G2M         </td><td>-0.255140251</td><td>0.19311215  </td><td>0.095520701 </td><td>⋯           </td><td>0.5659140   </td><td> 1.0364168  </td><td>-1.30456924 </td><td>-5.4628116  </td><td>-1.9005555  </td><td> 1.0018027  </td><td>-2.5513941  </td><td>-0.3275569  </td><td> 4.7457853  </td><td> 3.1820851  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240        </td><td>10058       </td><td>BTSC        </td><td>0.078047326 </td><td> 0.2119064  </td><td>-0.22082262 </td><td>S           </td><td> 0.432729024</td><td>0.10956094  </td><td>0.007489049 </td><td>⋯           </td><td>2.2555504   </td><td>-1.6466283  </td><td> 1.18936923 </td><td> 2.4664831  </td><td> 1.0089357  </td><td> 0.7129536  </td><td>-1.7719795  </td><td>-1.0386423  </td><td> 1.7797592  </td><td> 0.1388518  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337        </td><td>10798       </td><td>BTSC        </td><td>0.061863308 </td><td>-0.1322671  </td><td>-0.20464285 </td><td>G1          </td><td> 0.072375715</td><td>0.07132768  </td><td>0.000000000 </td><td>⋯           </td><td>2.6829950   </td><td>-3.1942123  </td><td>-1.53382334 </td><td> 1.0749278  </td><td> 3.2954332  </td><td>-2.6045133  </td><td>-0.3074987  </td><td>-0.6009474  </td><td> 1.2862561  </td><td>-2.3595514  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140        </td><td>14601       </td><td>BTSC        </td><td>0.081501267 </td><td> 0.4082643  </td><td> 0.40187978 </td><td>S           </td><td> 0.006384479</td><td>0.25821914  </td><td>0.007630352 </td><td>⋯           </td><td>2.6098085   </td><td> 0.4529079  </td><td>-1.02563061 </td><td> 0.5009398  </td><td> 2.2269572  </td><td>-1.2858385  </td><td>-2.4643894  </td><td>-0.9815629  </td><td>-0.7160265  </td><td>-0.7608637  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543        </td><td>  820       </td><td>BTSC        </td><td>0.108536585 </td><td>-0.1085568  </td><td>-0.14419756 </td><td>G1          </td><td> 0.035640733</td><td>0.04803011  </td><td>0.000000000 </td><td>⋯           </td><td>0.8944369   </td><td>-1.6530236  </td><td>-0.08019493 </td><td> 2.2799159  </td><td> 0.4626099  </td><td>-3.2290846  </td><td> 0.7074034  </td><td> 2.3155190  </td><td>-0.3342554  </td><td> 1.4165346  </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'nGene'</li>
	<li>'nUMI'</li>
	<li>'orig.ident'</li>
	<li>'percent.mito'</li>
	<li>'S.Score'</li>
	<li>'G2M.Score'</li>
	<li>'Phase'</li>
	<li>'CC.Difference'</li>
	<li>'Zhong_NPCs_upreg_AUC'</li>
	<li>'Zhong_Excitatory_neurons_upreg_AUC'</li>
	<li>'Zhong_Interneurons_upreg_AUC'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Astrocytes_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC1_upreg_AUC'</li>
	<li>'nowakowski_nEN.early2_upreg_AUC'</li>
	<li>'nowakowski_nEN.late_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.1_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.2_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC2_upreg_AUC'</li>
	<li>'nowakowski_nEN.early1_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC3_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.3_upreg_AUC'</li>
	<li>'nowakowski_IN.STR_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE2_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN2_upreg_AUC'</li>
	<li>'nowakowski_IPC.div1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN3_upreg_AUC'</li>
	<li>'nowakowski_nIN1_upreg_AUC'</li>
	<li>'nowakowski_nIN2_upreg_AUC'</li>
	<li>'nowakowski_nIN3_upreg_AUC'</li>
	<li>'nowakowski_nIN4_upreg_AUC'</li>
	<li>'nowakowski_nIN5_upreg_AUC'</li>
	<li>'nowakowski_Endothelial_upreg_AUC'</li>
	<li>'nowakowski_U4_upreg_AUC'</li>
	<li>'nowakowski_Mural_upreg_AUC'</li>
	<li>'nowakowski_Glyc_upreg_AUC'</li>
	<li>'nowakowski_U1_upreg_AUC'</li>
	<li>'nowakowski_Microglia_upreg_AUC'</li>
	<li>'nowakowski_U3_upreg_AUC'</li>
	<li>'nowakowski_Choroid_upreg_AUC'</li>
	<li>'nowakowski_RG.early_upreg_AUC'</li>
	<li>'nowakowski_U2_upreg_AUC'</li>
	<li>'nowakowski_oRG_upreg_AUC'</li>
	<li>'nowakowski_tRG_upreg_AUC'</li>
	<li>'nowakowski_vRG_upreg_AUC'</li>
	<li>'nowakowski_RG.div1_upreg_AUC'</li>
	<li>'nowakowski_OPC_upreg_AUC'</li>
	<li>'nowakowski_Astrocyte_upreg_AUC'</li>
	<li>'nowakowski_RG.div2_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG1_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG2_upreg_AUC'</li>
	<li>'nowakowski_MGE.div_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC1_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC2_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC3_upreg_AUC'</li>
	<li>'cahoy_astrocyte_AUC'</li>
	<li>'cahoy_oligodendrocyte_AUC'</li>
	<li>'cahoy_neuron_AUC'</li>
	<li>'cahoy_astro_young_AUC'</li>
	<li>'cahoy_astro_mature_AUC'</li>
	<li>'cahoy_OPC_AUC'</li>
	<li>'cahoy_OL_myel_AUC'</li>
	<li>'cahoy_MOG_pos_AUC'</li>
	<li>'cahoy_astro_in_vivo_AUC'</li>
	<li>'cahoy_cultured_astroglia_AUC'</li>
	<li>'a1.astro_AUC'</li>
	<li>'a2.astro_AUC'</li>
	<li>'mizrak_MicrogliaA_AUC'</li>
	<li>'mizrak_MicrogliaB_1_AUC'</li>
	<li>'mizrak_MicrogliaB_2_AUC'</li>
	<li>'mizrak_Endothelial1_AUC'</li>
	<li>'mizrak_Endothelial2_AUC'</li>
	<li>'mizrak_Endothelial3_AUC'</li>
	<li>'mizrak_Pericyte_AUC'</li>
	<li>'mizrak_Fibroblast_AUC'</li>
	<li>'mizrak_vSMC_AUC'</li>
	<li>'mizrak_Neuron1_AUC'</li>
	<li>'mizrak_Neuron2_AUC'</li>
	<li>'mizrak_Neuron3_AUC'</li>
	<li>'mizrak_Astro1_AUC'</li>
	<li>'mizrak_Astro2_AUC'</li>
	<li>'mizrak_Astro3_AUC'</li>
	<li>'mizrak_Astro4_AUC'</li>
	<li>'mizrak_Astro5_AUC'</li>
	<li>'mizrak_Ependymal_AUC'</li>
	<li>'mizrak_aNSC1_AUC'</li>
	<li>'mizrak_aNSC2_TAC_AUC'</li>
	<li>'mizrak_NB_AUC'</li>
	<li>'mizrak_lateral_astro_M_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_M_upreg_AUC'</li>
	<li>'mizrak_M_astro_lateral_upreg_AUC'</li>
	<li>'mizrak_F_astro_septal_upreg_AUC'</li>
	<li>'mizrak_lateral_astro_F_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_F_upreg_AUC'</li>
	<li>'mizrak_M_astro_septal_upreg_AUC'</li>
	<li>'mizrak_F_astro_lateral_upreg_AUC'</li>
	<li>'RNA.GSC.c1_AUC'</li>
	<li>'RNA.GSC.c2_AUC'</li>
	<li>'glioma.stem.cell_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_NEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC'</li>
	<li>'Suva_Stemness_AUC'</li>
	<li>'scGBM_AUC'</li>
	<li>'scBTSC_AUC'</li>
	<li>'Neftel_MES2_AUC'</li>
	<li>'Neftel_MES1_AUC'</li>
	<li>'Neftel_AC_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_NPC2_AUC'</li>
	<li>'Neftel_G1.S_AUC'</li>
	<li>'Neftel_G2.M_AUC'</li>
	<li>'SampleID'</li>
	<li>'C1_C2_diff'</li>
	<li>'SampleType'</li>
	<li>'Lab'</li>
	<li>'culture_cond'</li>
	<li>'CultureMethod'</li>
	<li>'res.0.8'</li>
	<li>'seurat_clusters'</li>
	<li>'is_637_L_800_L'</li>
	<li>'pred_class'</li>
	<li>'misclassified'</li>
	<li>'Class'</li>
	<li>'SUSD2'</li>
	<li>'SAA1'</li>
	<li>'CAPS'</li>
	<li>'CHI3L1'</li>
	<li>'AQP4'</li>
	<li>'GFAP'</li>
	<li>'ZFP36'</li>
	<li>'TPPP3'</li>
	<li>'NTRK2'</li>
	<li>'ID3'</li>
	<li>'HBB'</li>
	<li>'HILPDA'</li>
	<li>'CYR61'</li>
	<li>'CLU'</li>
	<li>'C9orf24'</li>
	<li>'HMGB2'</li>
	<li>'TYMS'</li>
	<li>'STMN1'</li>
	<li>'MLLT11'</li>
	<li>'VGF'</li>
	<li>'MDM2'</li>
	<li>'CCT6A'</li>
	<li>'DLL3'</li>
	<li>'STMN2'</li>
	<li>'TUBB'</li>
	<li>'CDK4'</li>
	<li>'NUP107'</li>
	<li>'SOX11'</li>
	<li>'SOX4'</li>
	<li>'HIST1H4C'</li>
	<li>'SEC61G'</li>
	<li>'CDKN2A'</li>
	<li>'OCIAD2'</li>
	<li>'H2AFZ'</li>
	<li>'NEFL'</li>
	<li>'PTMA'</li>
	<li>'HMGN2'</li>
	<li>'SET'</li>
	<li>'DEK'</li>
	<li>'MDK'</li>
	<li>'DKK1'</li>
	<li>'PPP1R14B'</li>
	<li>'HMGN1'</li>
	<li>'TTYH1'</li>
	<li>'HOPX'</li>
	<li>'COL1A2'</li>
	<li>'CRYAB'</li>
	<li>'NOVA1'</li>
	<li>'WSB1'</li>
	<li>'NMB'</li>
	<li>'RIC3'</li>
	<li>'AGT'</li>
	<li>'MEG3'</li>
	<li>'APOE'</li>
	<li>'NEAT1'</li>
	<li>'MALAT1'</li>
	<li>'PC1'</li>
	<li>'PC2'</li>
	<li>'PC3'</li>
	<li>'PC4'</li>
	<li>'PC5'</li>
	<li>'PC6'</li>
	<li>'PC7'</li>
	<li>'PC8'</li>
	<li>'PC9'</li>
	<li>'PC10'</li>
	<li>'PC11'</li>
	<li>'PC12'</li>
	<li>'PC13'</li>
	<li>'PC14'</li>
	<li>'PC15'</li>
	<li>'PC16'</li>
	<li>'PC17'</li>
	<li>'PC18'</li>
	<li>'PC19'</li>
	<li>'PC20'</li>
	<li>'PC21'</li>
	<li>'PC22'</li>
	<li>'PC23'</li>
	<li>'PC24'</li>
	<li>'PC25'</li>
	<li>'PC26'</li>
	<li>'PC27'</li>
	<li>'PC28'</li>
	<li>'PC29'</li>
	<li>'PC30'</li>
	<li>'PC31'</li>
	<li>'PC32'</li>
	<li>'PC33'</li>
	<li>'PC34'</li>
	<li>'PC35'</li>
	<li>'PC36'</li>
	<li>'PC37'</li>
	<li>'PC38'</li>
	<li>'PC39'</li>
	<li>'PC40'</li>
	<li>'PC41'</li>
	<li>'PC42'</li>
	<li>'PC43'</li>
	<li>'PC44'</li>
	<li>'PC45'</li>
	<li>'PC46'</li>
	<li>'PC47'</li>
	<li>'PC48'</li>
	<li>'PC49'</li>
	<li>'PC50'</li>
	<li>'PC51'</li>
	<li>'PC52'</li>
	<li>'PC53'</li>
	<li>'PC54'</li>
	<li>'PC55'</li>
	<li>'PC56'</li>
	<li>'PC57'</li>
	<li>'PC58'</li>
	<li>'PC59'</li>
	<li>'PC60'</li>
	<li>'PC61'</li>
	<li>'PC62'</li>
	<li>'PC63'</li>
	<li>'PC64'</li>
	<li>'PC65'</li>
	<li>'PC66'</li>
	<li>'PC67'</li>
	<li>'PC68'</li>
	<li>'PC69'</li>
	<li>'PC70'</li>
	<li>'PC71'</li>
	<li>'PC72'</li>
	<li>'PC73'</li>
	<li>'PC74'</li>
	<li>'PC75'</li>
	<li>'PC76'</li>
	<li>'PC77'</li>
	<li>'PC78'</li>
	<li>'PC79'</li>
	<li>'PC80'</li>
	<li>'PC81'</li>
	<li>'PC82'</li>
	<li>'PC83'</li>
	<li>'PC84'</li>
	<li>'PC85'</li>
	<li>'PC86'</li>
	<li>'PC87'</li>
	<li>'PC88'</li>
	<li>'PC89'</li>
	<li>'PC90'</li>
	<li>'PC91'</li>
	<li>'PC92'</li>
	<li>'PC93'</li>
	<li>'PC94'</li>
	<li>'PC95'</li>
	<li>'PC96'</li>
	<li>'PC97'</li>
	<li>'PC98'</li>
	<li>'PC99'</li>
	<li>'PC100'</li>
</ol>



### 1) Stacked bar chart of cell cycle phases


```R
p4 <- ggplot() + 
       geom_bar(aes(y = , 
                    x = Class, 
                    fill = product), 
                data = charts.data,          
                stat="identity")
```


```R
counts <- t(prop.table(table(dat$Class, dat$Phase), 1)
           )
counts <- counts[c("G1", "S", "G2M"), ]
counts

pdf("~/Desktop/CellCycle_GSC_T.pdf", width = 4, height = 4)

barplot(counts,
        ylab = "Proportion Cells (%)",
        col = c("#abdda4", "#2b83ba", "#d7191c"),
        las = 2
       )

dev.off()
```


         
                 GSC GSC-like_Tumour     Tumour Tumour-like_GSC
      G1  0.41555301      0.36566440 0.84503911      0.79531657
      S   0.31873430      0.45974782 0.09855908      0.06071119
      G2M 0.26571269      0.17458778 0.05640181      0.14397225



<strong>pdf:</strong> 2


### 3) Gene Expression violins


```R
colnames(dat)
```


<ol class=list-inline>
	<li>'nGene'</li>
	<li>'nUMI'</li>
	<li>'orig.ident'</li>
	<li>'percent.mito'</li>
	<li>'S.Score'</li>
	<li>'G2M.Score'</li>
	<li>'Phase'</li>
	<li>'CC.Difference'</li>
	<li>'Zhong_NPCs_upreg_AUC'</li>
	<li>'Zhong_Excitatory_neurons_upreg_AUC'</li>
	<li>'Zhong_Interneurons_upreg_AUC'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Astrocytes_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC1_upreg_AUC'</li>
	<li>'nowakowski_nEN.early2_upreg_AUC'</li>
	<li>'nowakowski_nEN.late_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.1_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.2_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC2_upreg_AUC'</li>
	<li>'nowakowski_nEN.early1_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC3_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.3_upreg_AUC'</li>
	<li>'nowakowski_IN.STR_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE2_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN2_upreg_AUC'</li>
	<li>'nowakowski_IPC.div1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN3_upreg_AUC'</li>
	<li>'nowakowski_nIN1_upreg_AUC'</li>
	<li>'nowakowski_nIN2_upreg_AUC'</li>
	<li>'nowakowski_nIN3_upreg_AUC'</li>
	<li>'nowakowski_nIN4_upreg_AUC'</li>
	<li>'nowakowski_nIN5_upreg_AUC'</li>
	<li>'nowakowski_Endothelial_upreg_AUC'</li>
	<li>'nowakowski_U4_upreg_AUC'</li>
	<li>'nowakowski_Mural_upreg_AUC'</li>
	<li>'nowakowski_Glyc_upreg_AUC'</li>
	<li>'nowakowski_U1_upreg_AUC'</li>
	<li>'nowakowski_Microglia_upreg_AUC'</li>
	<li>'nowakowski_U3_upreg_AUC'</li>
	<li>'nowakowski_Choroid_upreg_AUC'</li>
	<li>'nowakowski_RG.early_upreg_AUC'</li>
	<li>'nowakowski_U2_upreg_AUC'</li>
	<li>'nowakowski_oRG_upreg_AUC'</li>
	<li>'nowakowski_tRG_upreg_AUC'</li>
	<li>'nowakowski_vRG_upreg_AUC'</li>
	<li>'nowakowski_RG.div1_upreg_AUC'</li>
	<li>'nowakowski_OPC_upreg_AUC'</li>
	<li>'nowakowski_Astrocyte_upreg_AUC'</li>
	<li>'nowakowski_RG.div2_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG1_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG2_upreg_AUC'</li>
	<li>'nowakowski_MGE.div_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC1_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC2_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC3_upreg_AUC'</li>
	<li>'cahoy_astrocyte_AUC'</li>
	<li>'cahoy_oligodendrocyte_AUC'</li>
	<li>'cahoy_neuron_AUC'</li>
	<li>'cahoy_astro_young_AUC'</li>
	<li>'cahoy_astro_mature_AUC'</li>
	<li>'cahoy_OPC_AUC'</li>
	<li>'cahoy_OL_myel_AUC'</li>
	<li>'cahoy_MOG_pos_AUC'</li>
	<li>'cahoy_astro_in_vivo_AUC'</li>
	<li>'cahoy_cultured_astroglia_AUC'</li>
	<li>'a1.astro_AUC'</li>
	<li>'a2.astro_AUC'</li>
	<li>'mizrak_MicrogliaA_AUC'</li>
	<li>'mizrak_MicrogliaB_1_AUC'</li>
	<li>'mizrak_MicrogliaB_2_AUC'</li>
	<li>'mizrak_Endothelial1_AUC'</li>
	<li>'mizrak_Endothelial2_AUC'</li>
	<li>'mizrak_Endothelial3_AUC'</li>
	<li>'mizrak_Pericyte_AUC'</li>
	<li>'mizrak_Fibroblast_AUC'</li>
	<li>'mizrak_vSMC_AUC'</li>
	<li>'mizrak_Neuron1_AUC'</li>
	<li>'mizrak_Neuron2_AUC'</li>
	<li>'mizrak_Neuron3_AUC'</li>
	<li>'mizrak_Astro1_AUC'</li>
	<li>'mizrak_Astro2_AUC'</li>
	<li>'mizrak_Astro3_AUC'</li>
	<li>'mizrak_Astro4_AUC'</li>
	<li>'mizrak_Astro5_AUC'</li>
	<li>'mizrak_Ependymal_AUC'</li>
	<li>'mizrak_aNSC1_AUC'</li>
	<li>'mizrak_aNSC2_TAC_AUC'</li>
	<li>'mizrak_NB_AUC'</li>
	<li>'mizrak_lateral_astro_M_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_M_upreg_AUC'</li>
	<li>'mizrak_M_astro_lateral_upreg_AUC'</li>
	<li>'mizrak_F_astro_septal_upreg_AUC'</li>
	<li>'mizrak_lateral_astro_F_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_F_upreg_AUC'</li>
	<li>'mizrak_M_astro_septal_upreg_AUC'</li>
	<li>'mizrak_F_astro_lateral_upreg_AUC'</li>
	<li>'RNA.GSC.c1_AUC'</li>
	<li>'RNA.GSC.c2_AUC'</li>
	<li>'glioma.stem.cell_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_NEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC'</li>
	<li>'Suva_Stemness_AUC'</li>
	<li>'scGBM_AUC'</li>
	<li>'scBTSC_AUC'</li>
	<li>'Neftel_MES2_AUC'</li>
	<li>'Neftel_MES1_AUC'</li>
	<li>'Neftel_AC_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_NPC2_AUC'</li>
	<li>'Neftel_G1.S_AUC'</li>
	<li>'Neftel_G2.M_AUC'</li>
	<li>'SampleID'</li>
	<li>'C1_C2_diff'</li>
	<li>'SampleType'</li>
	<li>'Lab'</li>
	<li>'culture_cond'</li>
	<li>'CultureMethod'</li>
	<li>'res.0.8'</li>
	<li>'seurat_clusters'</li>
	<li>'is_637_L_800_L'</li>
	<li>'pred_class'</li>
	<li>'misclassified'</li>
	<li>'Class'</li>
	<li>'SUSD2'</li>
	<li>'SAA1'</li>
	<li>'CAPS'</li>
	<li>'CHI3L1'</li>
	<li>'AQP4'</li>
	<li>'GFAP'</li>
	<li>'ZFP36'</li>
	<li>'TPPP3'</li>
	<li>'NTRK2'</li>
	<li>'ID3'</li>
	<li>'HBB'</li>
	<li>'HILPDA'</li>
	<li>'CYR61'</li>
	<li>'CLU'</li>
	<li>'C9orf24'</li>
	<li>'HMGB2'</li>
	<li>'TYMS'</li>
	<li>'STMN1'</li>
	<li>'MLLT11'</li>
	<li>'VGF'</li>
	<li>'MDM2'</li>
	<li>'CCT6A'</li>
	<li>'DLL3'</li>
	<li>'STMN2'</li>
	<li>'TUBB'</li>
	<li>'CDK4'</li>
	<li>'NUP107'</li>
	<li>'SOX11'</li>
	<li>'SOX4'</li>
	<li>'HIST1H4C'</li>
	<li>'SEC61G'</li>
	<li>'CDKN2A'</li>
	<li>'OCIAD2'</li>
	<li>'H2AFZ'</li>
	<li>'NEFL'</li>
	<li>'PTMA'</li>
	<li>'HMGN2'</li>
	<li>'SET'</li>
	<li>'DEK'</li>
	<li>'MDK'</li>
	<li>'DKK1'</li>
	<li>'PPP1R14B'</li>
	<li>'HMGN1'</li>
	<li>'TTYH1'</li>
	<li>'HOPX'</li>
	<li>'COL1A2'</li>
	<li>'CRYAB'</li>
	<li>'NOVA1'</li>
	<li>'WSB1'</li>
	<li>'NMB'</li>
	<li>'RIC3'</li>
	<li>'AGT'</li>
	<li>'MEG3'</li>
	<li>'APOE'</li>
	<li>'NEAT1'</li>
	<li>'MALAT1'</li>
</ol>




```R
compare_means(AQP4 ~ Class,  data = dat)


ggboxplot(dat, x = "Class", 
          y = "HMGN1",
          color = "Class", 
          palette = "jco") 
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50)     # Add global p-value
```


<table>
<thead><tr><th scope=col>.y.</th><th scope=col>group1</th><th scope=col>group2</th><th scope=col>p</th><th scope=col>p.adj</th><th scope=col>p.format</th><th scope=col>p.signif</th><th scope=col>method</th></tr></thead>
<tbody>
	<tr><td>AQP4                                                   </td><td><span style=white-space:pre-wrap>GSC            </span></td><td>Tumour-like_GSC                                        </td><td> 0.000000e+00                                          </td><td> 0.000000e+00                                          </td><td>&lt;2e-16                                              </td><td>****                                                   </td><td>Wilcoxon                                               </td></tr>
	<tr><td>AQP4                                                   </td><td><span style=white-space:pre-wrap>GSC            </span></td><td><span style=white-space:pre-wrap>Tumour         </span></td><td> 0.000000e+00                                          </td><td> 0.000000e+00                                          </td><td>&lt;2e-16                                              </td><td>****                                                   </td><td>Wilcoxon                                               </td></tr>
	<tr><td>AQP4                                                   </td><td><span style=white-space:pre-wrap>GSC            </span></td><td>GSC-like_Tumour                                        </td><td> 0.000000e+00                                          </td><td> 0.000000e+00                                          </td><td>&lt;2e-16                                              </td><td>****                                                   </td><td>Wilcoxon                                               </td></tr>
	<tr><td>AQP4                                                   </td><td>Tumour-like_GSC                                        </td><td><span style=white-space:pre-wrap>Tumour         </span></td><td>1.353945e-222                                          </td><td>2.707891e-222                                          </td><td>&lt;2e-16                                              </td><td>****                                                   </td><td>Wilcoxon                                               </td></tr>
	<tr><td>AQP4           </td><td>Tumour-like_GSC</td><td>GSC-like_Tumour</td><td> 2.079656e-01  </td><td> 2.079656e-01  </td><td>0.21           </td><td>ns             </td><td>Wilcoxon       </td></tr>
	<tr><td>AQP4                                                   </td><td><span style=white-space:pre-wrap>Tumour         </span></td><td>GSC-like_Tumour                                        </td><td> 0.000000e+00                                          </td><td> 0.000000e+00                                          </td><td>&lt;2e-16                                              </td><td>****                                                   </td><td>Wilcoxon                                               </td></tr>
</tbody>
</table>






![png](output_21_2.png)


---
### 5) Color PCA plot by gene expression
---


```R
library(RColorBrewer)
```


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/")

dat <- readRDS("logistic_regression_T_v_L_metadata_DEgenes.rds")
head(dat)
colnames(dat)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640        </td><td>  875       </td><td>BTSC        </td><td>0.043428571 </td><td> 0.0957859  </td><td> 0.07303696 </td><td>S           </td><td> 0.022748941</td><td>0.10663091  </td><td>0.000000000 </td><td>⋯           </td><td>1.3949986   </td><td>-1.2290755  </td><td> 0.04078296 </td><td>-0.3940755  </td><td>-4.3912124  </td><td>-0.5042649  </td><td> 3.3774399  </td><td> 7.9740340  </td><td>-0.9635511  </td><td>-3.9862054  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036        </td><td> 2408       </td><td>BTSC        </td><td>0.002076412 </td><td> 0.0535880  </td><td> 0.30872825 </td><td>G2M         </td><td>-0.255140251</td><td>0.19311215  </td><td>0.095520701 </td><td>⋯           </td><td>0.5659140   </td><td> 1.0364168  </td><td>-1.30456924 </td><td>-5.4628116  </td><td>-1.9005555  </td><td> 1.0018027  </td><td>-2.5513941  </td><td>-0.3275569  </td><td> 4.7457853  </td><td> 3.1820851  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240        </td><td>10058       </td><td>BTSC        </td><td>0.078047326 </td><td> 0.2119064  </td><td>-0.22082262 </td><td>S           </td><td> 0.432729024</td><td>0.10956094  </td><td>0.007489049 </td><td>⋯           </td><td>2.2555504   </td><td>-1.6466283  </td><td> 1.18936923 </td><td> 2.4664831  </td><td> 1.0089357  </td><td> 0.7129536  </td><td>-1.7719795  </td><td>-1.0386423  </td><td> 1.7797592  </td><td> 0.1388518  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337        </td><td>10798       </td><td>BTSC        </td><td>0.061863308 </td><td>-0.1322671  </td><td>-0.20464285 </td><td>G1          </td><td> 0.072375715</td><td>0.07132768  </td><td>0.000000000 </td><td>⋯           </td><td>2.6829950   </td><td>-3.1942123  </td><td>-1.53382334 </td><td> 1.0749278  </td><td> 3.2954332  </td><td>-2.6045133  </td><td>-0.3074987  </td><td>-0.6009474  </td><td> 1.2862561  </td><td>-2.3595514  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140        </td><td>14601       </td><td>BTSC        </td><td>0.081501267 </td><td> 0.4082643  </td><td> 0.40187978 </td><td>S           </td><td> 0.006384479</td><td>0.25821914  </td><td>0.007630352 </td><td>⋯           </td><td>2.6098085   </td><td> 0.4529079  </td><td>-1.02563061 </td><td> 0.5009398  </td><td> 2.2269572  </td><td>-1.2858385  </td><td>-2.4643894  </td><td>-0.9815629  </td><td>-0.7160265  </td><td>-0.7608637  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543        </td><td>  820       </td><td>BTSC        </td><td>0.108536585 </td><td>-0.1085568  </td><td>-0.14419756 </td><td>G1          </td><td> 0.035640733</td><td>0.04803011  </td><td>0.000000000 </td><td>⋯           </td><td>0.8944369   </td><td>-1.6530236  </td><td>-0.08019493 </td><td> 2.2799159  </td><td> 0.4626099  </td><td>-3.2290846  </td><td> 0.7074034  </td><td> 2.3155190  </td><td>-0.3342554  </td><td> 1.4165346  </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'nGene'</li>
	<li>'nUMI'</li>
	<li>'orig.ident'</li>
	<li>'percent.mito'</li>
	<li>'S.Score'</li>
	<li>'G2M.Score'</li>
	<li>'Phase'</li>
	<li>'CC.Difference'</li>
	<li>'Zhong_NPCs_upreg_AUC'</li>
	<li>'Zhong_Excitatory_neurons_upreg_AUC'</li>
	<li>'Zhong_Interneurons_upreg_AUC'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Astrocytes_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC1_upreg_AUC'</li>
	<li>'nowakowski_nEN.early2_upreg_AUC'</li>
	<li>'nowakowski_nEN.late_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.1_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.2_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC2_upreg_AUC'</li>
	<li>'nowakowski_nEN.early1_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC3_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.3_upreg_AUC'</li>
	<li>'nowakowski_IN.STR_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE2_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN2_upreg_AUC'</li>
	<li>'nowakowski_IPC.div1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN3_upreg_AUC'</li>
	<li>'nowakowski_nIN1_upreg_AUC'</li>
	<li>'nowakowski_nIN2_upreg_AUC'</li>
	<li>'nowakowski_nIN3_upreg_AUC'</li>
	<li>'nowakowski_nIN4_upreg_AUC'</li>
	<li>'nowakowski_nIN5_upreg_AUC'</li>
	<li>'nowakowski_Endothelial_upreg_AUC'</li>
	<li>'nowakowski_U4_upreg_AUC'</li>
	<li>'nowakowski_Mural_upreg_AUC'</li>
	<li>'nowakowski_Glyc_upreg_AUC'</li>
	<li>'nowakowski_U1_upreg_AUC'</li>
	<li>'nowakowski_Microglia_upreg_AUC'</li>
	<li>'nowakowski_U3_upreg_AUC'</li>
	<li>'nowakowski_Choroid_upreg_AUC'</li>
	<li>'nowakowski_RG.early_upreg_AUC'</li>
	<li>'nowakowski_U2_upreg_AUC'</li>
	<li>'nowakowski_oRG_upreg_AUC'</li>
	<li>'nowakowski_tRG_upreg_AUC'</li>
	<li>'nowakowski_vRG_upreg_AUC'</li>
	<li>'nowakowski_RG.div1_upreg_AUC'</li>
	<li>'nowakowski_OPC_upreg_AUC'</li>
	<li>'nowakowski_Astrocyte_upreg_AUC'</li>
	<li>'nowakowski_RG.div2_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG1_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG2_upreg_AUC'</li>
	<li>'nowakowski_MGE.div_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC1_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC2_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC3_upreg_AUC'</li>
	<li>'cahoy_astrocyte_AUC'</li>
	<li>'cahoy_oligodendrocyte_AUC'</li>
	<li>'cahoy_neuron_AUC'</li>
	<li>'cahoy_astro_young_AUC'</li>
	<li>'cahoy_astro_mature_AUC'</li>
	<li>'cahoy_OPC_AUC'</li>
	<li>'cahoy_OL_myel_AUC'</li>
	<li>'cahoy_MOG_pos_AUC'</li>
	<li>'cahoy_astro_in_vivo_AUC'</li>
	<li>'cahoy_cultured_astroglia_AUC'</li>
	<li>'a1.astro_AUC'</li>
	<li>'a2.astro_AUC'</li>
	<li>'mizrak_MicrogliaA_AUC'</li>
	<li>'mizrak_MicrogliaB_1_AUC'</li>
	<li>'mizrak_MicrogliaB_2_AUC'</li>
	<li>'mizrak_Endothelial1_AUC'</li>
	<li>'mizrak_Endothelial2_AUC'</li>
	<li>'mizrak_Endothelial3_AUC'</li>
	<li>'mizrak_Pericyte_AUC'</li>
	<li>'mizrak_Fibroblast_AUC'</li>
	<li>'mizrak_vSMC_AUC'</li>
	<li>'mizrak_Neuron1_AUC'</li>
	<li>'mizrak_Neuron2_AUC'</li>
	<li>'mizrak_Neuron3_AUC'</li>
	<li>'mizrak_Astro1_AUC'</li>
	<li>'mizrak_Astro2_AUC'</li>
	<li>'mizrak_Astro3_AUC'</li>
	<li>'mizrak_Astro4_AUC'</li>
	<li>'mizrak_Astro5_AUC'</li>
	<li>'mizrak_Ependymal_AUC'</li>
	<li>'mizrak_aNSC1_AUC'</li>
	<li>'mizrak_aNSC2_TAC_AUC'</li>
	<li>'mizrak_NB_AUC'</li>
	<li>'mizrak_lateral_astro_M_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_M_upreg_AUC'</li>
	<li>'mizrak_M_astro_lateral_upreg_AUC'</li>
	<li>'mizrak_F_astro_septal_upreg_AUC'</li>
	<li>'mizrak_lateral_astro_F_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_F_upreg_AUC'</li>
	<li>'mizrak_M_astro_septal_upreg_AUC'</li>
	<li>'mizrak_F_astro_lateral_upreg_AUC'</li>
	<li>'RNA.GSC.c1_AUC'</li>
	<li>'RNA.GSC.c2_AUC'</li>
	<li>'glioma.stem.cell_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_NEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC'</li>
	<li>'Suva_Stemness_AUC'</li>
	<li>'scGBM_AUC'</li>
	<li>'scBTSC_AUC'</li>
	<li>'Neftel_MES2_AUC'</li>
	<li>'Neftel_MES1_AUC'</li>
	<li>'Neftel_AC_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_NPC2_AUC'</li>
	<li>'Neftel_G1.S_AUC'</li>
	<li>'Neftel_G2.M_AUC'</li>
	<li>'SampleID'</li>
	<li>'C1_C2_diff'</li>
	<li>'SampleType'</li>
	<li>'Lab'</li>
	<li>'culture_cond'</li>
	<li>'CultureMethod'</li>
	<li>'res.0.8'</li>
	<li>'seurat_clusters'</li>
	<li>'is_637_L_800_L'</li>
	<li>'pred_class'</li>
	<li>'misclassified'</li>
	<li>'Class'</li>
	<li>'SUSD2'</li>
	<li>'SAA1'</li>
	<li>'CAPS'</li>
	<li>'CHI3L1'</li>
	<li>'AQP4'</li>
	<li>'GFAP'</li>
	<li>'ZFP36'</li>
	<li>'TPPP3'</li>
	<li>'NTRK2'</li>
	<li>'ID3'</li>
	<li>'HBB'</li>
	<li>'HILPDA'</li>
	<li>'CYR61'</li>
	<li>'CLU'</li>
	<li>'C9orf24'</li>
	<li>'HMGB2'</li>
	<li>'TYMS'</li>
	<li>'STMN1'</li>
	<li>'MLLT11'</li>
	<li>'VGF'</li>
	<li>'MDM2'</li>
	<li>'CCT6A'</li>
	<li>'DLL3'</li>
	<li>'STMN2'</li>
	<li>'TUBB'</li>
	<li>'CDK4'</li>
	<li>'NUP107'</li>
	<li>'SOX11'</li>
	<li>'SOX4'</li>
	<li>'HIST1H4C'</li>
	<li>'SEC61G'</li>
	<li>'CDKN2A'</li>
	<li>'OCIAD2'</li>
	<li>'H2AFZ'</li>
	<li>'NEFL'</li>
	<li>'PTMA'</li>
	<li>'HMGN2'</li>
	<li>'SET'</li>
	<li>'DEK'</li>
	<li>'MDK'</li>
	<li>'DKK1'</li>
	<li>'PPP1R14B'</li>
	<li>'HMGN1'</li>
	<li>'TTYH1'</li>
	<li>'HOPX'</li>
	<li>'COL1A2'</li>
	<li>'CRYAB'</li>
	<li>'NOVA1'</li>
	<li>'WSB1'</li>
	<li>'NMB'</li>
	<li>'RIC3'</li>
	<li>'AGT'</li>
	<li>'MEG3'</li>
	<li>'APOE'</li>
	<li>'NEAT1'</li>
	<li>'MALAT1'</li>
	<li>'PC1'</li>
	<li>'PC2'</li>
	<li>'PC3'</li>
	<li>'PC4'</li>
	<li>'PC5'</li>
	<li>'PC6'</li>
	<li>'PC7'</li>
	<li>'PC8'</li>
	<li>'PC9'</li>
	<li>'PC10'</li>
	<li>'PC11'</li>
	<li>'PC12'</li>
	<li>'PC13'</li>
	<li>'PC14'</li>
	<li>'PC15'</li>
	<li>'PC16'</li>
	<li>'PC17'</li>
	<li>'PC18'</li>
	<li>'PC19'</li>
	<li>'PC20'</li>
	<li>'PC21'</li>
	<li>'PC22'</li>
	<li>'PC23'</li>
	<li>'PC24'</li>
	<li>'PC25'</li>
	<li>'PC26'</li>
	<li>'PC27'</li>
	<li>'PC28'</li>
	<li>'PC29'</li>
	<li>'PC30'</li>
	<li>'PC31'</li>
	<li>'PC32'</li>
	<li>'PC33'</li>
	<li>'PC34'</li>
	<li>'PC35'</li>
	<li>'PC36'</li>
	<li>'PC37'</li>
	<li>'PC38'</li>
	<li>'PC39'</li>
	<li>'PC40'</li>
	<li>'PC41'</li>
	<li>'PC42'</li>
	<li>'PC43'</li>
	<li>'PC44'</li>
	<li>'PC45'</li>
	<li>'PC46'</li>
	<li>'PC47'</li>
	<li>'PC48'</li>
	<li>'PC49'</li>
	<li>'PC50'</li>
	<li>'PC51'</li>
	<li>'PC52'</li>
	<li>'PC53'</li>
	<li>'PC54'</li>
	<li>'PC55'</li>
	<li>'PC56'</li>
	<li>'PC57'</li>
	<li>'PC58'</li>
	<li>'PC59'</li>
	<li>'PC60'</li>
	<li>'PC61'</li>
	<li>'PC62'</li>
	<li>'PC63'</li>
	<li>'PC64'</li>
	<li>'PC65'</li>
	<li>'PC66'</li>
	<li>'PC67'</li>
	<li>'PC68'</li>
	<li>'PC69'</li>
	<li>'PC70'</li>
	<li>'PC71'</li>
	<li>'PC72'</li>
	<li>'PC73'</li>
	<li>'PC74'</li>
	<li>'PC75'</li>
	<li>'PC76'</li>
	<li>'PC77'</li>
	<li>'PC78'</li>
	<li>'PC79'</li>
	<li>'PC80'</li>
	<li>'PC81'</li>
	<li>'PC82'</li>
	<li>'PC83'</li>
	<li>'PC84'</li>
	<li>'PC85'</li>
	<li>'PC86'</li>
	<li>'PC87'</li>
	<li>'PC88'</li>
	<li>'PC89'</li>
	<li>'PC90'</li>
	<li>'PC91'</li>
	<li>'PC92'</li>
	<li>'PC93'</li>
	<li>'PC94'</li>
	<li>'PC95'</li>
	<li>'PC96'</li>
	<li>'PC97'</li>
	<li>'PC98'</li>
	<li>'PC99'</li>
	<li>'PC100'</li>
</ol>




```R
plot(dat$PC1, dat$PC2)
```


![png](output_25_0.png)



```R
plot_GSC_t_hex <- function(dat, x, y, sig, name){ 


plot.title <-  paste(sig, "\n", name, "  ", nrow(dat), "cells")
    
p <- ggplot(dat, aes_string(x=x, y=y, z=sig)) +
 stat_summary_hex(bins=100, fun = "median") +
 theme_classic(base_size=8) +
 scale_fill_gradientn("Median Raw \nAUC Score", colours = c(brewer.pal(n = 8, name = "YlOrRd"))) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
ggtitle(plot.title) + labs(z = "Expression")


    
    }
```


```R
x <- "PC1"
y <- "PC2"
name <- "GSCs+Tumour"
sigs <- colnames(dat)[133:(133+55)]
sigs
```


<ol class=list-inline>
	<li>'SUSD2'</li>
	<li>'SAA1'</li>
	<li>'CAPS'</li>
	<li>'CHI3L1'</li>
	<li>'AQP4'</li>
	<li>'GFAP'</li>
	<li>'ZFP36'</li>
	<li>'TPPP3'</li>
	<li>'NTRK2'</li>
	<li>'ID3'</li>
	<li>'HBB'</li>
	<li>'HILPDA'</li>
	<li>'CYR61'</li>
	<li>'CLU'</li>
	<li>'C9orf24'</li>
	<li>'HMGB2'</li>
	<li>'TYMS'</li>
	<li>'STMN1'</li>
	<li>'MLLT11'</li>
	<li>'VGF'</li>
	<li>'MDM2'</li>
	<li>'CCT6A'</li>
	<li>'DLL3'</li>
	<li>'STMN2'</li>
	<li>'TUBB'</li>
	<li>'CDK4'</li>
	<li>'NUP107'</li>
	<li>'SOX11'</li>
	<li>'SOX4'</li>
	<li>'HIST1H4C'</li>
	<li>'SEC61G'</li>
	<li>'CDKN2A'</li>
	<li>'OCIAD2'</li>
	<li>'H2AFZ'</li>
	<li>'NEFL'</li>
	<li>'PTMA'</li>
	<li>'HMGN2'</li>
	<li>'SET'</li>
	<li>'DEK'</li>
	<li>'MDK'</li>
	<li>'DKK1'</li>
	<li>'PPP1R14B'</li>
	<li>'HMGN1'</li>
	<li>'TTYH1'</li>
	<li>'HOPX'</li>
	<li>'COL1A2'</li>
	<li>'CRYAB'</li>
	<li>'NOVA1'</li>
	<li>'WSB1'</li>
	<li>'NMB'</li>
	<li>'RIC3'</li>
	<li>'AGT'</li>
	<li>'MEG3'</li>
	<li>'APOE'</li>
	<li>'NEAT1'</li>
	<li>'MALAT1'</li>
</ol>




```R
p <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
p[[i]] <- plot_GSC_t_hex(dat, x, y, sig, name)
                                        
}
```


```R
ml <- marrangeGrob(p, nrow=3, ncol=3)
ggsave("~/Desktop/GSCs_Tumours_genes_hex.pdf", ml, height = 15, width = 15, dpi = 72, limitsize = FALSE)

```

---
# 4.0 Plot classifer cells on PCA
---


```R
library(ggplot2)
library(gridExtra)
library(ggpubr)
```

    Warning message:
    “package ‘ggplot2’ was built under R version 3.4.4”Warning message:
    “package ‘ggpubr’ was built under R version 3.4.4”Loading required package: magrittr



```R
dat <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_metadata.rds")
```


```R
head(dat$pred_class) #line
head(dat$misclassified) #TRUE 

table(dat$pred_class == "Line" & dat$misclassified == TRUE)
table(dat$pred_class == "Tumour" & dat$misclassified == TRUE)
```


<ol class=list-inline>
	<li>'Line'</li>
	<li>'Line'</li>
	<li>'Line'</li>
	<li>'Line'</li>
	<li>'Line'</li>
	<li>'Tumour'</li>
</ol>




<ol class=list-inline>
	<li>FALSE</li>
	<li>FALSE</li>
	<li>FALSE</li>
	<li>FALSE</li>
	<li>FALSE</li>
	<li>TRUE</li>
</ol>




    
    FALSE  TRUE 
    77800  2062 



    
    FALSE  TRUE 
    78709  1153 



```R
p <- ggplot(data=dat) +
 geom_density_2d(aes(x=PC1, y=PC2, color=SampleType), bins = 15, binwidth = 1) + 
scale_colour_manual(values=c("#225ea8", "black", "grey")) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
theme(legend.position = "none") + theme(text = element_text(size=12))+

geom_point(data = dat[dat$pred_class == "Line" & dat$misclassified == TRUE ,], aes(x=PC1, y=PC2),
          color = "black", alpha = 0.3) + ggtitle("GSC-like Tumour Cells (n=2062 cells)")



q <- ggplot(data=dat) +
 geom_density_2d(aes(x=PC1, y=PC2, color=SampleType), bins = 15, binwidth = 1) + 
scale_colour_manual(values=c("#225ea8", "black", "grey")) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
theme(legend.position = "none") + theme(text = element_text(size=12))+

geom_point(data = dat[dat$pred_class == "Tumour" & dat$misclassified == TRUE ,], aes(x=PC1, y=PC2),
          color = "#225ea8", alpha = 0.3) + ggtitle("Tumour-like GSC Cells (n=1153 cells)")

p
q
```






![png](output_34_2.png)



![png](output_34_3.png)



```R
pdf("~/Desktop/PCA_classifer.pdf", width = 10, height = 5)
grid.arrange(p,
             q,
             ncol = 2
            )

dev.off()
```


<strong>pdf:</strong> 2


### Now plot with contour lines only and classified dots




```R
dat <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_metadata.rds")
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>is_637_L_800_L</th><th scope=col>pred_class</th><th scope=col>misclassified</th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640        </td><td>  875       </td><td>BTSC        </td><td>0.043428571 </td><td> 0.0957859  </td><td> 0.07303696 </td><td>S           </td><td> 0.022748941</td><td>0.10663091  </td><td>0.000000000 </td><td>⋯           </td><td>21          </td><td>21          </td><td>other       </td><td>Line        </td><td>FALSE       </td><td>-12.7289371 </td><td>-0.5538633  </td><td>-20.410265  </td><td> -7.732371  </td><td>-14.582780  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036        </td><td> 2408       </td><td>BTSC        </td><td>0.002076412 </td><td> 0.0535880  </td><td> 0.30872825 </td><td>G2M         </td><td>-0.255140251</td><td>0.19311215  </td><td>0.095520701 </td><td>⋯           </td><td>21          </td><td>21          </td><td>other       </td><td>Line        </td><td>FALSE       </td><td>  0.9425536 </td><td> 1.6079858  </td><td>-13.909896  </td><td>-13.610787  </td><td>-20.616752  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240        </td><td>10058       </td><td>BTSC        </td><td>0.078047326 </td><td> 0.2119064  </td><td>-0.22082262 </td><td>S           </td><td> 0.432729024</td><td>0.10956094  </td><td>0.007489049 </td><td>⋯           </td><td>4           </td><td>4           </td><td>other       </td><td>Line        </td><td>FALSE       </td><td>  4.2237210 </td><td>17.3566342  </td><td>  2.027612  </td><td>  1.856217  </td><td> -0.944526  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337        </td><td>10798       </td><td>BTSC        </td><td>0.061863308 </td><td>-0.1322671  </td><td>-0.20464285 </td><td>G1          </td><td> 0.072375715</td><td>0.07132768  </td><td>0.000000000 </td><td>⋯           </td><td>13          </td><td>13          </td><td>other       </td><td>Line        </td><td>FALSE       </td><td>-14.3646978 </td><td>12.1673063  </td><td> -2.144286  </td><td> -1.250715  </td><td> 10.269212  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140        </td><td>14601       </td><td>BTSC        </td><td>0.081501267 </td><td> 0.4082643  </td><td> 0.40187978 </td><td>S           </td><td> 0.006384479</td><td>0.25821914  </td><td>0.007630352 </td><td>⋯           </td><td>4           </td><td>4           </td><td>other       </td><td>Line        </td><td>FALSE       </td><td>  7.0178336 </td><td>23.3358358  </td><td> 13.651359  </td><td> -2.276415  </td><td>  2.269168  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543        </td><td>  820       </td><td>BTSC        </td><td>0.108536585 </td><td>-0.1085568  </td><td>-0.14419756 </td><td>G1          </td><td> 0.035640733</td><td>0.04803011  </td><td>0.000000000 </td><td>⋯           </td><td>21          </td><td>21          </td><td>other       </td><td>Tumour      </td><td> TRUE       </td><td>-14.7980172 </td><td>-8.7067006  </td><td>-22.498126  </td><td>-12.294728  </td><td> -9.675785  </td></tr>
</tbody>
</table>




```R
df <- dat
```


```R
###Calcualte contours

set.seed(1001)

#define 95% and 50% contour levels for GSCs
GSC <- data.frame(x=df[df$SampleType == "Line", "PC1"], 
                  y=df[df$SampleType == "Line", "PC2"]
                 )

kd <- ks::kde(GSC, compute.cont=TRUE)
GSC_contour_95 <- with(kd, contourLines(x=eval.points[[1]], 
                                        y=eval.points[[2]],
                                        z=estimate, 
                                        levels=cont["1%"])[[1]]
                      )
GSC_contour_95 <- data.frame(GSC_contour_95)

GSC_contour_50 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
GSC_contour_50 <- data.frame(GSC_contour_50)


###### Tumours

set.seed(1001)

#define 95% and 50% contour levels for Tumour
Tumour <- data.frame(x=df[df$SampleType == "Tumour", "PC1"], 
                  y=df[df$SampleType == "Tumour", "PC2"]
                 )

kd <- ks::kde(Tumour, compute.cont=TRUE)
Tumour_contour_95 <- with(kd, contourLines(x=eval.points[[1]], 
                                        y=eval.points[[2]],
                                        z=estimate, 
                                        levels=cont["1%"])[[1]]
                      )
Tumour_contour_95 <- data.frame(Tumour_contour_95)

Tumour_contour_50 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
Tumour_contour_50 <- data.frame(Tumour_contour_50)
```


```R

```


```R
    
p <- ggplot(data=dat) +
 #geom_density_2d(aes(x=PC1, y=PC2, color=SampleType), bins = 15, binwidth = 1) + 
scale_colour_manual(values=c("#225ea8", "black", "grey")) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
theme(legend.position = "none") + theme(text = element_text(size=12))+

geom_point(data = dat[dat$pred_class == "Line" & dat$misclassified == TRUE ,], aes(x=PC1, y=PC2),
          color = "black", fill = "grey",  alpha = 1, pch=21) + ggtitle("GSC-like Tumour Cells (n=2062 cells)")+

geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 1
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 1
             ) 



```


```R
    
q <- ggplot(data=dat) +
 #geom_density_2d(aes(x=PC1, y=PC2, color=SampleType), bins = 15, binwidth = 1) + 
scale_colour_manual(values=c("#225ea8", "black", "grey")) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
theme(legend.position = "none") + theme(text = element_text(size=12))+

geom_point(data = dat[dat$pred_class == "Tumour" & dat$misclassified == TRUE ,], aes(x=PC1, y=PC2),
          color = "#225ea8", fill = "white",  alpha = 1, pch=21) + ggtitle("Tumour-like GSCs (n=1153 cells)")+
geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 1
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 1
             ) 
    
```


```R
pdf("~/Desktop/PCA_classifer_noContour.pdf", width = 10, height = 5)
grid.arrange(p,
             q,
             ncol = 2
            )

dev.off()
```


<strong>pdf:</strong> 2



```R
p
```




![png](output_44_1.png)



```R

```


```R

```


```R

```


```R

```

 -----
 ## Plot proportion of misclassified cells across samples
 ----


```R
meta <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_metadata_DEgenes.rds")
head(meta)
#colnames(meta)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640        </td><td>  875       </td><td>BTSC        </td><td>0.043428571 </td><td> 0.0957859  </td><td> 0.07303696 </td><td>S           </td><td> 0.022748941</td><td>0.10663091  </td><td>0.000000000 </td><td>⋯           </td><td>1.3949986   </td><td>-1.2290755  </td><td> 0.04078296 </td><td>-0.3940755  </td><td>-4.3912124  </td><td>-0.5042649  </td><td> 3.3774399  </td><td> 7.9740340  </td><td>-0.9635511  </td><td>-3.9862054  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036        </td><td> 2408       </td><td>BTSC        </td><td>0.002076412 </td><td> 0.0535880  </td><td> 0.30872825 </td><td>G2M         </td><td>-0.255140251</td><td>0.19311215  </td><td>0.095520701 </td><td>⋯           </td><td>0.5659140   </td><td> 1.0364168  </td><td>-1.30456924 </td><td>-5.4628116  </td><td>-1.9005555  </td><td> 1.0018027  </td><td>-2.5513941  </td><td>-0.3275569  </td><td> 4.7457853  </td><td> 3.1820851  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240        </td><td>10058       </td><td>BTSC        </td><td>0.078047326 </td><td> 0.2119064  </td><td>-0.22082262 </td><td>S           </td><td> 0.432729024</td><td>0.10956094  </td><td>0.007489049 </td><td>⋯           </td><td>2.2555504   </td><td>-1.6466283  </td><td> 1.18936923 </td><td> 2.4664831  </td><td> 1.0089357  </td><td> 0.7129536  </td><td>-1.7719795  </td><td>-1.0386423  </td><td> 1.7797592  </td><td> 0.1388518  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337        </td><td>10798       </td><td>BTSC        </td><td>0.061863308 </td><td>-0.1322671  </td><td>-0.20464285 </td><td>G1          </td><td> 0.072375715</td><td>0.07132768  </td><td>0.000000000 </td><td>⋯           </td><td>2.6829950   </td><td>-3.1942123  </td><td>-1.53382334 </td><td> 1.0749278  </td><td> 3.2954332  </td><td>-2.6045133  </td><td>-0.3074987  </td><td>-0.6009474  </td><td> 1.2862561  </td><td>-2.3595514  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140        </td><td>14601       </td><td>BTSC        </td><td>0.081501267 </td><td> 0.4082643  </td><td> 0.40187978 </td><td>S           </td><td> 0.006384479</td><td>0.25821914  </td><td>0.007630352 </td><td>⋯           </td><td>2.6098085   </td><td> 0.4529079  </td><td>-1.02563061 </td><td> 0.5009398  </td><td> 2.2269572  </td><td>-1.2858385  </td><td>-2.4643894  </td><td>-0.9815629  </td><td>-0.7160265  </td><td>-0.7608637  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543        </td><td>  820       </td><td>BTSC        </td><td>0.108536585 </td><td>-0.1085568  </td><td>-0.14419756 </td><td>G1          </td><td> 0.035640733</td><td>0.04803011  </td><td>0.000000000 </td><td>⋯           </td><td>0.8944369   </td><td>-1.6530236  </td><td>-0.08019493 </td><td> 2.2799159  </td><td> 0.4626099  </td><td>-3.2290846  </td><td> 0.7074034  </td><td> 2.3155190  </td><td>-0.3342554  </td><td> 1.4165346  </td></tr>
</tbody>
</table>




```R
prop <- prop.table(table(meta$SampleID, meta$Class), 1)
head(prop)

prop_L <- prop[grep("_L", rownames(prop)), ]
prop_L <- prop_L[order(prop_L[,1 ], decreasing = T), ]
prop_L

prop_T <- prop[grep("_T", rownames(prop)), ]
prop_T <- prop_T[order(prop_T[,3 ], decreasing = T), ]
prop_T
```


             
                      GSC GSC-like_Tumour      Tumour Tumour-like_GSC
      BT127_L 0.890581717     0.000000000 0.000000000     0.109418283
      BT147_L 0.981438515     0.000000000 0.000000000     0.018561485
      BT48_L  0.994546694     0.000000000 0.000000000     0.005453306
      BT67_L  0.994176373     0.000000000 0.000000000     0.005823627
      BT73_L  0.989690722     0.000000000 0.000000000     0.010309278
      BT84_L  0.871937639     0.000000000 0.000000000     0.128062361



             
                       GSC GSC-like_Tumour       Tumour Tumour-like_GSC
      BT89_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G549_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G583_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G797_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G837_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G851_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G885_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G945_L  1.0000000000    0.0000000000 0.0000000000    0.0000000000
      G523_L  0.9994659546    0.0000000000 0.0000000000    0.0005340454
      G564_L  0.9994604802    0.0000000000 0.0000000000    0.0005395198
      G729_L  0.9993057966    0.0000000000 0.0000000000    0.0006942034
      G876_L  0.9976019185    0.0000000000 0.0000000000    0.0023980815
      BT94_L  0.9959116926    0.0000000000 0.0000000000    0.0040883074
      G637_L  0.9954622802    0.0000000000 0.0000000000    0.0045377198
      G895_L  0.9946483180    0.0000000000 0.0000000000    0.0053516820
      BT48_L  0.9945466939    0.0000000000 0.0000000000    0.0054533061
      BT67_L  0.9941763727    0.0000000000 0.0000000000    0.0058236273
      BT73_L  0.9896907216    0.0000000000 0.0000000000    0.0103092784
      BT147_L 0.9814385151    0.0000000000 0.0000000000    0.0185614849
      G799_L  0.9792843691    0.0000000000 0.0000000000    0.0207156309
      G946_L  0.9450980392    0.0000000000 0.0000000000    0.0549019608
      G566_L  0.9192343604    0.0000000000 0.0000000000    0.0807656396
      G620_L  0.9040404040    0.0000000000 0.0000000000    0.0959595960
      BT127_L 0.8905817175    0.0000000000 0.0000000000    0.1094182825
      BT84_L  0.8719376392    0.0000000000 0.0000000000    0.1280623608



             
                    GSC GSC-like_Tumour    Tumour Tumour-like_GSC
      G1003_T 0.0000000       0.0187500 0.9812500       0.0000000
      G945_T  0.0000000       0.1595300 0.8404700       0.0000000
      G946_T  0.0000000       0.2419929 0.7580071       0.0000000
      G983_T  0.0000000       0.2697228 0.7302772       0.0000000
      G910_T  0.0000000       0.3697674 0.6302326       0.0000000
      G967_T  0.0000000       0.6626506 0.3373494       0.0000000
      G620_T  0.0000000       0.7914530 0.2085470       0.0000000



```R
dff <- rbind(prop_L, prop_T)
```


```R
pdf("~/Desktop/PropCells_Classifier.pdf", width = 6, height = 3)

barplot(t(dff),
        las = 2,
        col = c("#225ea8", "grey", "black", "white"),
        ylab = "Proportion Cells (%)",
        cex.names = 0.6,
        cex.axis = 0.6,
        cex = 0.6
       )

dev.off()
```


<strong>pdf:</strong> 2

