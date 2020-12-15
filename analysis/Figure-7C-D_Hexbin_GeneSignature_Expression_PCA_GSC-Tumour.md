
---
# Hexbin PCA GSC-Tumour Paper Figures
---
L.Richards  


```R
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(RColorBrewer)
```

---
## 1.0 Expression of gene signatures on PCA plot
---



```R
pc.file <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/BTSC_TumourCell_G800Removed_PCA_allGenes.Rdata"
meta.file <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/BTSC_TumourCell_G800Removed_PCA_metdata.Rdata"
aucell <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/AUCell/BTSC_Tumour_G800Removed_AUCellScores.Rdata"

load(pc.file)  
load(meta.file)
load(aucell)

file.prefix <- "G800L_Removed"

a <- strsplit(rownames(BTSC_Tumour_AUCell), "_")
a <- matrix(unlist(a), ncol = 4, byrow = TRUE)
#head(a)

BTSC_Tumour_AUCell$SampleID <- paste(a[,2], a[,3], sep = "_")
BTSC_Tumour_AUCell$SampleType <- ifelse(grepl("_L", rownames(BTSC_Tumour_AUCell)), "GSC", "Tumour")
head(BTSC_Tumour_AUCell)

dff <- data.frame(pc@cell.embeddings)
dff <- cbind(dff, BTSC_Tumour_AUCell)
head(dff)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>Neftel_MES2_AUC</th><th scope=col>Neftel_MES1_AUC</th><th scope=col>Neftel_AC_AUC</th><th scope=col>Neftel_OPC_AUC</th><th scope=col>Neftel_NPC1_AUC</th><th scope=col>Neftel_NPC2_AUC</th><th scope=col>Neftel_G1.S_AUC</th><th scope=col>Neftel_G2.M_AUC</th><th scope=col>SampleID</th><th scope=col>SampleType</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td> 640        </td><td>  875       </td><td>BTSC        </td><td>0.043428571 </td><td> 0.0957859  </td><td> 0.07303696 </td><td>S           </td><td> 0.022748941</td><td>0.10663091  </td><td>0.000000000 </td><td>⋯           </td><td>0.05034014  </td><td>0.05806973  </td><td>0.09693621  </td><td>0.05047956  </td><td>0.03978561  </td><td>0.08501340  </td><td>0.120975227 </td><td>0.07175179  </td><td>BT127_L     </td><td>GSC         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>1036        </td><td> 2408       </td><td>BTSC        </td><td>0.002076412 </td><td> 0.0535880  </td><td> 0.30872825 </td><td>G2M         </td><td>-0.255140251</td><td>0.19311215  </td><td>0.095520701 </td><td>⋯           </td><td>0.07200577  </td><td>0.12873168  </td><td>0.11406593  </td><td>0.02158001  </td><td>0.09678417  </td><td>0.07025356  </td><td>0.109480974 </td><td>0.12056452  </td><td>BT127_L     </td><td>GSC         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>3240        </td><td>10058       </td><td>BTSC        </td><td>0.078047326 </td><td> 0.2119064  </td><td>-0.22082262 </td><td>S           </td><td> 0.432729024</td><td>0.10956094  </td><td>0.007489049 </td><td>⋯           </td><td>0.12234591  </td><td>0.16545730  </td><td>0.17809088  </td><td>0.11423103  </td><td>0.15718409  </td><td>0.11074005  </td><td>0.048511459 </td><td>0.03487903  </td><td>BT127_L     </td><td>GSC         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>3337        </td><td>10798       </td><td>BTSC        </td><td>0.061863308 </td><td>-0.1322671  </td><td>-0.20464285 </td><td>G1          </td><td> 0.072375715</td><td>0.07132768  </td><td>0.000000000 </td><td>⋯           </td><td>0.16720264  </td><td>0.19417888  </td><td>0.18160670  </td><td>0.08821302  </td><td>0.21583179  </td><td>0.19239332  </td><td>0.005675734 </td><td>0.06377688  </td><td>BT127_L     </td><td>GSC         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>4140        </td><td>14601       </td><td>BTSC        </td><td>0.081501267 </td><td> 0.4082643  </td><td> 0.40187978 </td><td>S           </td><td> 0.006384479</td><td>0.25821914  </td><td>0.007630352 </td><td>⋯           </td><td>0.14623789  </td><td>0.18199090  </td><td>0.17457506  </td><td>0.06713781  </td><td>0.09196042  </td><td>0.08276644  </td><td>0.263047048 </td><td>0.10450269  </td><td>BT127_L     </td><td>GSC         </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td> 543        </td><td>  820       </td><td>BTSC        </td><td>0.108536585 </td><td>-0.1085568  </td><td>-0.14419756 </td><td>G1          </td><td> 0.035640733</td><td>0.04803011  </td><td>0.000000000 </td><td>⋯           </td><td>0.15947227  </td><td>0.14025265  </td><td>0.20820006  </td><td>0.12483173  </td><td>0.06507937  </td><td>0.03673469  </td><td>0.000000000 </td><td>0.05219534  </td><td>BT127_L     </td><td>GSC         </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>PC9</th><th scope=col>PC10</th><th scope=col>⋯</th><th scope=col>Neftel_MES2_AUC</th><th scope=col>Neftel_MES1_AUC</th><th scope=col>Neftel_AC_AUC</th><th scope=col>Neftel_OPC_AUC</th><th scope=col>Neftel_NPC1_AUC</th><th scope=col>Neftel_NPC2_AUC</th><th scope=col>Neftel_G1.S_AUC</th><th scope=col>Neftel_G2.M_AUC</th><th scope=col>SampleID</th><th scope=col>SampleType</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>-12.7289371</td><td>-0.5538633 </td><td>-20.410265 </td><td> -7.732371 </td><td>-14.582780 </td><td>-22.689954 </td><td>-29.986860 </td><td> 7.0352326 </td><td>-2.5351561 </td><td> 1.9578116 </td><td>⋯          </td><td>0.05034014 </td><td>0.05806973 </td><td>0.09693621 </td><td>0.05047956 </td><td>0.03978561 </td><td>0.08501340 </td><td>0.120975227</td><td>0.07175179 </td><td>BT127_L    </td><td>GSC        </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>  0.9425536</td><td> 1.6079858 </td><td>-13.909896 </td><td>-13.610787 </td><td>-20.616752 </td><td> -2.893296 </td><td>  9.973656 </td><td>-0.1489149 </td><td>-0.9948348 </td><td> 9.9273695 </td><td>⋯          </td><td>0.07200577 </td><td>0.12873168 </td><td>0.11406593 </td><td>0.02158001 </td><td>0.09678417 </td><td>0.07025356 </td><td>0.109480974</td><td>0.12056452 </td><td>BT127_L    </td><td>GSC        </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>  4.2237210</td><td>17.3566342 </td><td>  2.027612 </td><td>  1.856217 </td><td> -0.944526 </td><td> -4.670833 </td><td>  3.018915 </td><td>-1.6568515 </td><td> 0.4951982 </td><td>-0.2887847 </td><td>⋯          </td><td>0.12234591 </td><td>0.16545730 </td><td>0.17809088 </td><td>0.11423103 </td><td>0.15718409 </td><td>0.11074005 </td><td>0.048511459</td><td>0.03487903 </td><td>BT127_L    </td><td>GSC        </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>-14.3646978</td><td>12.1673063 </td><td> -2.144286 </td><td> -1.250715 </td><td> 10.269212 </td><td>  5.800723 </td><td> -5.748215 </td><td> 0.8159342 </td><td> 6.4224062 </td><td>-9.6798437 </td><td>⋯          </td><td>0.16720264 </td><td>0.19417888 </td><td>0.18160670 </td><td>0.08821302 </td><td>0.21583179 </td><td>0.19239332 </td><td>0.005675734</td><td>0.06377688 </td><td>BT127_L    </td><td>GSC        </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>  7.0178336</td><td>23.3358358 </td><td> 13.651359 </td><td> -2.276415 </td><td>  2.269168 </td><td> -9.896081 </td><td> -6.497069 </td><td> 7.2201086 </td><td>-4.4535921 </td><td> 7.2550430 </td><td>⋯          </td><td>0.14623789 </td><td>0.18199090 </td><td>0.17457506 </td><td>0.06713781 </td><td>0.09196042 </td><td>0.08276644 </td><td>0.263047048</td><td>0.10450269 </td><td>BT127_L    </td><td>GSC        </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>-14.7980172</td><td>-8.7067006 </td><td>-22.498126 </td><td>-12.294728 </td><td> -9.675785 </td><td> -2.589288 </td><td>-17.117835 </td><td> 3.0031946 </td><td>-5.0322937 </td><td> 1.1042951 </td><td>⋯          </td><td>0.15947227 </td><td>0.14025265 </td><td>0.20820006 </td><td>0.12483173 </td><td>0.06507937 </td><td>0.03673469 </td><td>0.000000000</td><td>0.05219534 </td><td>BT127_L    </td><td>GSC        </td></tr>
</tbody>
</table>




```R
#### filter sigs to be only ones we care about

sigs <- colnames(dff)[grep("AUC", colnames(dff))]
sigs <- sigs[c(grep("RNA", sigs),
               grep("a1", sigs),
               grep("a2", sigs),
          grep("Neftel", sigs),
          grep("Zhong", sigs),
          grep("cahoy", sigs),
          grep("mizrak", sigs),
          grep("VERHAAK", sigs)
         )]


df <- dff[ ,c("PC1", "PC2", "SampleType", sigs)]
colnames(df)[4:7] <- c("Developmental_GSC_AUC", 
               "InjuryResponse_GSC_AUC", 
               "Liddelow_A1_ReactiveAstrocytes_AUC",
               "Liddelow_A2_ReactiveAstrocytes_AUC"
              )

head(df)

```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>SampleType</th><th scope=col>Developmental_GSC_AUC</th><th scope=col>InjuryResponse_GSC_AUC</th><th scope=col>Liddelow_A1_ReactiveAstrocytes_AUC</th><th scope=col>Liddelow_A2_ReactiveAstrocytes_AUC</th><th scope=col>Neftel_MES2_AUC</th><th scope=col>Neftel_MES1_AUC</th><th scope=col>Neftel_AC_AUC</th><th scope=col>⋯</th><th scope=col>mizrak_M_astro_lateral_upreg_AUC</th><th scope=col>mizrak_F_astro_septal_upreg_AUC</th><th scope=col>mizrak_lateral_astro_F_upreg_AUC</th><th scope=col>mizrak_septal_astrocytes_F_upreg_AUC</th><th scope=col>mizrak_M_astro_septal_upreg_AUC</th><th scope=col>mizrak_F_astro_lateral_upreg_AUC</th><th scope=col>VERHAAK_GLIOBLASTOMA_NEURAL_AUC</th><th scope=col>VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC</th><th scope=col>VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC</th><th scope=col>VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>-12.7289371</td><td>-0.5538633 </td><td>GSC        </td><td>0.1833756  </td><td>0.1691589  </td><td>0.09818548 </td><td>0.04099462 </td><td>0.05034014 </td><td>0.05806973 </td><td>0.09693621 </td><td>⋯          </td><td>0.07904693 </td><td>0.09643124 </td><td>0.07731756 </td><td>0.3355331  </td><td>0.13455172 </td><td>0.07596196 </td><td>0.007391415</td><td>0.05141122 </td><td>0.02613838 </td><td>0.04673622 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>  0.9425536</td><td> 1.6079858 </td><td>GSC        </td><td>0.1225561  </td><td>0.1561081  </td><td>0.06594982 </td><td>0.13263889 </td><td>0.07200577 </td><td>0.12873168 </td><td>0.11406593 </td><td>⋯          </td><td>0.21127084 </td><td>0.03361996 </td><td>0.01329931 </td><td>0.2280545  </td><td>0.05327586 </td><td>0.26628706 </td><td>0.033468800</td><td>0.04211048 </td><td>0.02051541 </td><td>0.01696094 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>  4.2237210</td><td>17.3566342 </td><td>GSC        </td><td>0.1403795  </td><td>0.1532729  </td><td>0.04478047 </td><td>0.10506272 </td><td>0.12234591 </td><td>0.16545730 </td><td>0.17809088 </td><td>⋯          </td><td>0.24026060 </td><td>0.05267147 </td><td>0.02247632 </td><td>0.2437539  </td><td>0.09027586 </td><td>0.30144094 </td><td>0.049005165</td><td>0.06270099 </td><td>0.03009132 </td><td>0.03191546 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>-14.3646978</td><td>12.1673063 </td><td>GSC        </td><td>0.1531699  </td><td>0.1659836  </td><td>0.17163978 </td><td>0.16758513 </td><td>0.16720264 </td><td>0.19417888 </td><td>0.18160670 </td><td>⋯          </td><td>0.23380747 </td><td>0.07009065 </td><td>0.01761790 </td><td>0.2715267  </td><td>0.09082759 </td><td>0.27062478 </td><td>0.036643807</td><td>0.09168758 </td><td>0.04939123 </td><td>0.02886570 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>  7.0178336</td><td>23.3358358 </td><td>GSC        </td><td>0.1330963  </td><td>0.1653948  </td><td>0.10385305 </td><td>0.14056900 </td><td>0.14623789 </td><td>0.18199090 </td><td>0.17457506 </td><td>⋯          </td><td>0.24030193 </td><td>0.04624324 </td><td>0.01599843 </td><td>0.2331167  </td><td>0.07724138 </td><td>0.29493037 </td><td>0.043654221</td><td>0.05192552 </td><td>0.03003003 </td><td>0.04364633 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>-14.7980172</td><td>-8.7067006 </td><td>GSC        </td><td>0.1597614  </td><td>0.2143644  </td><td>0.05864695 </td><td>0.15152330 </td><td>0.15947227 </td><td>0.14025265 </td><td>0.20820006 </td><td>⋯          </td><td>0.10231817 </td><td>0.09272743 </td><td>0.06127006 </td><td>0.2707708  </td><td>0.04217241 </td><td>0.12904384 </td><td>0.029472526</td><td>0.03041794 </td><td>0.05860963 </td><td>0.04783975 </td></tr>
</tbody>
</table>




```R
###Calcualte contours

set.seed(1001)

#define 95% and 50% contour levels for GSCs
GSC <- data.frame(x=df[df$SampleType == "GSC", "PC1"], 
                  y=df[df$SampleType == "GSC", "PC2"]
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
plot_GSC_t_hex <- function(dat, x, y, sig, name){ 


plot.title <-  paste(sig, "\n", name, "  ", nrow(dat), "cells")
    
p <- ggplot(dat, aes_string(x=x, y=y, z=sig)) +
 stat_summary_hex(bins=100, fun = "median") +
 theme_classic(base_size=8) +
 scale_fill_gradientn("Median Raw \nAUC Score", colours = c(brewer.pal(n = 8, name = "YlOrRd")), limits=c(min(df[ ,sig]),max(df[ ,sig]))) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
ggtitle(plot.title) + labs(z = "Raw AUC") +  xlim(c(min(df[ ,"PC1"])-1, max(df[ ,"PC1"])+1)) + 
    ylim(c(min(df[ ,"PC2"])-1,max(df[ ,"PC2"])+1))


 p + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    }

```


```R
dat <- df[df$SampleType == "GSC", ]
datt <- df[df$SampleType == "Tumour", ]
x <- "PC1"
y <- "PC2"
name <- "GSCs+Tumour"
sigs <- colnames(df)[grep("AUC", colnames(df))]
sigs
```


<ol class=list-inline>
	<li>'Developmental_GSC_AUC'</li>
	<li>'InjuryResponse_GSC_AUC'</li>
	<li>'Liddelow_A1_ReactiveAstrocytes_AUC'</li>
	<li>'Liddelow_A2_ReactiveAstrocytes_AUC'</li>
	<li>'Neftel_MES2_AUC'</li>
	<li>'Neftel_MES1_AUC'</li>
	<li>'Neftel_AC_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_NPC2_AUC'</li>
	<li>'Neftel_G1.S_AUC'</li>
	<li>'Neftel_G2.M_AUC'</li>
	<li>'Zhong_NPCs_upreg_AUC'</li>
	<li>'Zhong_Excitatory_neurons_upreg_AUC'</li>
	<li>'Zhong_Interneurons_upreg_AUC'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Astrocytes_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
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
	<li>'VERHAAK_GLIOBLASTOMA_NEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC'</li>
</ol>




```R
p <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
p[[i]] <- plot_GSC_t_hex(dat, x, y, sig, name)
                                        
}


q <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
q[[i]] <- plot_GSC_t_hex(datt, x, y, sig, name)
                                        
}
```


```R
pdf("~/Desktop/PCA_GSC_T_sigs.pdf", width = 10, height = 5)

for (i in 1:length(p)){
    
   dd <- ggarrange(p[[i]],
          q[[i]], 
          ncol=2,
          common.legend = T,
          legend = "right"
         )
    print(dd)
    
}

dev.off()
```


<strong>pdf:</strong> 2


---
## 2.0 Expression of top and bottom loading genes on PCA plot
---
Extract top and bottom loading genes from PC1 and PC2
Use hexbin, calculate median expression within each hex.  
Plot with contours of where tumour and GSCs would be.   
Plot each sample type separately. 
library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
       )

dat <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_seurat.rds")


gene.load <- data.frame(dat@dr$pc@gene.loadings)


PC1_top10 <- rownames(gene.load[order(gene.load$PC1, decreasing = TRUE), ])[1:10]
PC1_bottom10 <- rownames(gene.load[order(gene.load$PC1, decreasing = FALSE), ])[1:10]

PC2_top10 <- rownames(gene.load[order(gene.load$PC2, decreasing = TRUE), ])[1:10]
PC2_bottom10 <- rownames(gene.load[order(gene.load$PC2, decreasing = FALSE), ])[1:10]

genes <- list(PC1_top10,
           PC1_bottom10,
           PC2_top10,
           PC2_bottom10
           )

names(genes) <- c('PC1_top10',
                  'PC1_bottom10',
                  'PC2_top10',
                  'PC2_bottom10'
                 )

saveRDS(genes, file = "Genes_GSC_T_PCA_top10.rds")



genes <- as.character(unlist(genes))

subset <- t(data.matrix(dat@data[genes, ]))
pc <- dat@dr$pca@cell.embeddings
meta <- dat@meta.data
meta.genes <- cbind(meta, subset, pc)
head(meta.genes)
saveRDS(meta.genes, file = "logistic_regression_T_v_L_metadata_PCAgenes.rds")


```R
genes <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/Genes_GSC_T_PCA_top10.rds")
dat_gene <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/Classifier_owen/logistic_regression_T_v_L_metadata_PCAgenes.rds")

```


```R
plot_GSC_t_hex <- function(dat, x, y, sig, name){ 


plot.title <-  paste(sig, "\n", name, "  ", nrow(dat), "cells")
    
p <- ggplot(dat, aes_string(x=x, y=y, z=sig)) +
 stat_summary_hex(bins=100, fun = "median") +
 theme_classic(base_size=8) +
 scale_fill_gradientn("Normalized\nGene\nExpression", colours = c(brewer.pal(n = 8, name = "YlOrRd")), limits=c(min(dat_gene[ ,sig]),max(dat_gene[ ,sig]))) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
ggtitle(plot.title) + labs(z = "Raw AUC") +  xlim(c(min(dat_gene[ ,"PC1"])-1, max(dat_gene[ ,"PC1"])+1)) + 
    ylim(c(min(dat_gene[ ,"PC2"])-1,max(dat_gene[ ,"PC2"])+1))


 p + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    }

```


```R
line <- dat_gene[dat_gene$SampleType == "Line", ]
tumour <- dat_gene[dat_gene$SampleType == "Tumour", ]
x <- "PC1"
y <- "PC2"
name <- "GSCs+Tumour"
sigs <- as.character(unlist(genes))
sigs
```


<ol class=list-inline>
	<li>'PFN1'</li>
	<li>'HMGA1'</li>
	<li>'PPP1R14B'</li>
	<li>'PRELID1'</li>
	<li>'TXN'</li>
	<li>'H2AFZ'</li>
	<li>'C12orf75'</li>
	<li>'MRPL33'</li>
	<li>'SIVA1'</li>
	<li>'CDKN3'</li>
	<li>'CLU'</li>
	<li>'FOS'</li>
	<li>'GFAP'</li>
	<li>'MT3'</li>
	<li>'C1orf61'</li>
	<li>'PLTP'</li>
	<li>'S100B'</li>
	<li>'TTYH1'</li>
	<li>'HSPA1A'</li>
	<li>'CST3'</li>
	<li>'FGFBP3'</li>
	<li>'PAFAH1B3'</li>
	<li>'MAP2'</li>
	<li>'STMN1'</li>
	<li>'CCND2'</li>
	<li>'KCNQ2'</li>
	<li>'OLIG1'</li>
	<li>'PTN'</li>
	<li>'MARCKSL1'</li>
	<li>'OLIG2'</li>
	<li>'OCIAD2'</li>
	<li>'S100A10'</li>
	<li>'ANXA1'</li>
	<li>'HSPB1'</li>
	<li>'CAV1'</li>
	<li>'PLP2'</li>
	<li>'C9orf24'</li>
	<li>'CD44'</li>
	<li>'ZFP36L1'</li>
	<li>'IGFBP7'</li>
</ol>




```R
p <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
p[[i]] <- plot_GSC_t_hex(line, x, y, sig, name)
                                        
}


q <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
q[[i]] <- plot_GSC_t_hex(tumour, x, y, sig, name)
                                        
}
```


```R
pdf("~/Desktop/PCA_GSC_T_Top10_PCGenes.pdf", width = 10, height = 5)

for (i in 1:length(p)){
    
   dd <- ggarrange(p[[i]],
          q[[i]], 
          ncol=2,
          common.legend = T,
          legend = "right"
         )
    print(dd)
    
}

dev.off()
```


<strong>pdf:</strong> 2



```R

```
