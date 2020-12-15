
---
# Correlate cell type signatures to PC1 and PC2 on GSC-T analysis
---
L.Richards  
Same procedure as GSC only analysis in 'PCAPlots_GSCs_Remove_G800L_Nov2019' notebook


```R
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(gtable)
library(ks)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(Seurat)
```

---
## 1.0 Correlate PC loadings to cell types
---




```R
## load object with PC embeddings and AUC scores

pc.file <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/BTSC_TumourCell_G800Removed_PCA_allGenes.Rdata"
aucell <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/AUCell/BTSC_Tumour_G800Removed_AUCellScores.Rdata"

load(pc.file)  
load(aucell)

dat <- cbind(data.frame(pc@cell.embeddings), BTSC_Tumour_AUCell)
head(dat)
```

### Correlate cell type scores to PC1 cell embeddings


```R
cols <- c(grep("^PC1$", colnames(dat)),
          grep("AUC", colnames(dat))
         )

corr.input <- dat[ ,cols]
#head(corr.input)

res <- cor(corr.input, method = "spearman")
#head(res)

test <- as.matrix(res["PC1", ])
colnames(test) <- "SpearmanCorr"
#head(test)

cat(test, sep = "\n")
```

    1
    0.4983853
    0.369861
    -0.1297525
    -0.7441121
    -0.7729154
    -0.361571
    -0.3020671
    -0.3782562
    0.1610409
    -0.2404352
    0.2339584
    0.2297741
    0.3010769
    0.09588073
    0.01812613
    -0.1194779
    0.114049
    -0.004356358
    0.08645941
    0.1671193
    -0.3193611
    0.1277149
    0.4674822
    0.5392047
    -0.05358472
    -0.1611358
    0.1990237
    -0.2744688
    0.07164461
    -0.2540353
    0.1536771
    -0.4416253
    -0.04724409
    0.2024375
    -0.00166785
    -0.2269942
    0.4303661
    0.1204036
    0.4150012
    0.4122734
    -0.618486
    -0.4692893
    -0.6654496
    0.1055039
    -0.6788956
    -0.8088874
    0.1101665
    -0.6583352
    0.2920855
    0.5296432
    0.4950143
    0.4293216
    0.1968858
    -0.7007616
    -0.1035149
    0.03137179
    0.3759063
    -0.5299804
    -0.3712721
    -0.4733416
    -0.6404959
    -0.782017
    0.6310415
    -0.2711062
    0.2644512
    -0.3346885
    -0.5005175
    0.1174335
    -0.1438626
    -0.4724073
    -0.07696726
    -0.02833114
    -0.2758116
    0.4414207
    -0.5516582
    0.0793313
    0.4812211
    0.1284705
    0.3485726
    -0.710829
    -0.2511188
    -0.7324021
    0.1342745
    0.5017999
    0.4262031
    0.1243766
    -0.2068082
    -0.2393718
    -0.2199414
    -0.6240748
    0.06709088
    -0.4800341
    -0.01773238
    -0.078661
    -0.808716
    0.633446
    0.347561
    -0.12117
    -0.4956153
    0.4388963
    -0.6079586
    -0.1942492
    -0.841172
    0.7585452
    -0.02833718
    0.2417216
    -0.7537688
    -0.6360853
    -0.3257717
    -0.0898876
    0.483431
    0.2613653



```R
bb <- read.csv("~/Desktop/AllSpearmanCorr_filtered_GSC_Tumour_PC1.csv")
bb <- bb[-33,]
head(bb)
tail(bb)
table(bb$Source)
```


<table>
<thead><tr><th scope=col>GeneSignature</th><th scope=col>SpearmanCorrleation</th><th scope=col>NewName</th><th scope=col>Source</th></tr></thead>
<tbody>
	<tr><td>nowakowski_Astrocyte_upreg_AUC</td><td>-0.8088874                    </td><td>Astrocytes                    </td><td>Nowakowski                    </td></tr>
	<tr><td>RNA.GSC.c1_AUC                </td><td>-0.8087160                    </td><td>Developmental_GSC_Program     </td><td>InHouse                       </td></tr>
	<tr><td>cahoy_astro_in_vivo_AUC       </td><td>-0.7820170                    </td><td>InVivo-astrocytes             </td><td>Cahoy                         </td></tr>
	<tr><td>Zhong_Astrocytes_upreg_AUC    </td><td>-0.7729154                    </td><td>Asctrocytes                   </td><td>Zhong                         </td></tr>
	<tr><td>Neftel_AC_AUC                 </td><td>-0.7537688                    </td><td>AC-GBM                        </td><td>Neftel                        </td></tr>
	<tr><td>Zhong_OPC_upreg_AUC           </td><td>-0.7441121                    </td><td>OPCs                          </td><td>Zhong                         </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>GeneSignature</th><th scope=col>SpearmanCorrleation</th><th scope=col>NewName</th><th scope=col>Source</th></tr></thead>
<tbody>
	<tr><th scope=row>27</th><td>nowakowski_MGE.IPC1_upreg_AUC</td><td>0.4950143                    </td><td>MGE-IPC1                     </td><td>Nowakowski                   </td></tr>
	<tr><th scope=row>28</th><td>Zhong_NPCs_upreg_AUC         </td><td>0.4983853                    </td><td>NPCs                         </td><td>Zhong                        </td></tr>
	<tr><th scope=row>29</th><td>nowakowski_MGE.div_upreg_AUC </td><td>0.5296432                    </td><td>MGE-div                      </td><td>Nowakowski                   </td></tr>
	<tr><th scope=row>30</th><td>nowakowski_IPC.div2_upreg_AUC</td><td>0.5392047                    </td><td>IPC-div2                     </td><td>Nowakowski                   </td></tr>
	<tr><th scope=row>31</th><td>cahoy_cultured_astroglia_AUC </td><td>0.6310415                    </td><td>Cultures-astroglia           </td><td>Cahoy                        </td></tr>
	<tr><th scope=row>32</th><td>RNA.GSC.c2_AUC               </td><td>0.6334460                    </td><td>InjuryResponse_GSC_Program   </td><td>InHouse                      </td></tr>
</tbody>
</table>





                    Cahoy    InHouse     Neftel Nowakowski    Verhaak      Zhong
             0          6          2          3         15          3          3



```R
cols <- c()
cols[grep("Zhong", bb$Source)] <- "#1b9e77"
cols[grep("Neftel", bb$Source)] <- "#d95f02"
cols[grep("Verhaak", bb$Source)] <- "#7570b3"
cols[grep("Nowakowski", bb$Source)] <- "#666666"
cols[grep("Cahoy", bb$Source)] <- "white"
cols[grep("Liddelow", bb$Source)] <- "black"
cols[grep("InHouse", bb$Source)] <- "red"
cols
```


<ol class=list-inline>
	<li>'#666666'</li>
	<li>'red'</li>
	<li>'white'</li>
	<li>'#1b9e77'</li>
	<li>'#d95f02'</li>
	<li>'#1b9e77'</li>
	<li>'white'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'white'</li>
	<li>'#d95f02'</li>
	<li>'#666666'</li>
	<li>'#7570b3'</li>
	<li>'white'</li>
	<li>'#7570b3'</li>
	<li>'white'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#7570b3'</li>
	<li>'#666666'</li>
	<li>'#d95f02'</li>
	<li>'#666666'</li>
	<li>'#1b9e77'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'white'</li>
	<li>'red'</li>
</ol>




```R


pdf("~/Desktop/celltype.GSC_Tumour_PC1.pdf", width = 10, height = 7)

par(mar=c(4.1, 13, 4.1, 11), xpd=TRUE)

barplot(c(as.numeric(bb$SpearmanCorrleation)),
        horiz = T,
        names.arg = as.character(bb$NewName),
        las= 1,
        xlim = c(-1,1),
        cex.names = 0.8,
         col = cols,
         xlab = "Spearman Correlation Coefficient"
       )


dev.off()
```


<strong>pdf:</strong> 2


---
### Correlate cell type scores to PC2 cell embeddings


```R
cols <- c(grep("^PC2$", colnames(dat)),
          grep("AUC", colnames(dat))
         )

corr.input <- dat[ ,cols]
#head(corr.input)

res <- cor(corr.input, method = "spearman")
#head(res)

test <- as.matrix(res["PC2", ])
colnames(test) <- "SpearmanCorr"
#head(test)


```

    1
    0.2093883
    -0.1337608
    -0.01738562
    0.4447152
    0.1493202
    -0.4740875
    -0.02824179
    0.08951598
    0.2428456
    0.2690887
    0.2928444
    0.2638019
    0.4779848
    0.2797117
    0.3348782
    -0.1992915
    0.3885552
    0.2177019
    0.4700269
    0.1807017
    0.4291823
    0.3690896
    0.4123759
    0.4165157
    0.4768571
    0.3109257
    0.3468325
    -0.1484875
    0.2592433
    0.3553652
    -0.3228165
    0.3420677
    -0.197335
    -0.03680531
    -0.2039841
    -0.3068358
    0.1889728
    -0.1807188
    0.2376986
    0.03294864
    0.04836984
    0.0004406445
    0.117244
    0.1787531
    0.4550529
    0.2122165
    0.3548195
    0.0656518
    0.479981
    0.4189413
    0.3280789
    0.3980771
    0.5707849
    -0.09995479
    0.1705601
    -0.001829194
    0.455536
    -0.3330018
    0.3242444
    -0.156232
    -0.1193097
    0.02632078
    -0.1899464
    -0.3175978
    -0.3879371
    0.09881122
    0.1657623
    0.03523933
    -0.0859448
    -0.2062624
    -0.1830098
    0.02828221
    -0.1202957
    -0.1198756
    -0.1307055
    0.3005241
    0.2444179
    0.262095
    0.3557408
    0.2630621
    0.2923705
    0.03399347
    -0.4779068
    0.4146296
    0.1876264
    0.5222532
    0.4178228
    0.1367378
    0.3216925
    -0.1206136
    -0.3051117
    -0.1595951
    0.07064447
    0.2570977
    0.3850271
    -0.5559915
    0.2815494
    -0.03853899
    0.5847428
    -0.4804564
    0.2567823
    0.3213978
    -0.0745432
    0.06682802
    -0.2830367
    -0.4773132
    0.2456283
    0.4450821
    0.6253201
    0.4038164
    0.2763517
    0.09423859



```R
bb <- read.csv("~/Desktop/AllSpearmanCorr_filtered_GSC_Tumour_PC2.csv")
bb <- bb[-c(21:22),]
head(bb)
tail(bb)
table(bb$Source)
```


<table>
<thead><tr><th scope=col>GeneSignature</th><th scope=col>SpearmanCorrleation</th><th scope=col>NewName</th><th scope=col>Source</th></tr></thead>
<tbody>
	<tr><td>RNA.GSC.c2_AUC                      </td><td>-0.5559915                          </td><td>InjuryResponse_GSC_Program          </td><td>InHouse                             </td></tr>
	<tr><td>VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC</td><td>-0.4804564                          </td><td>Mesenchymal-TCGA-GBM                </td><td>Verhaak                             </td></tr>
	<tr><td>Neftel_MES1_AUC                     </td><td>-0.4773132                          </td><td>MES1-GBM                            </td><td>Neftel                              </td></tr>
	<tr><td>Zhong_Microglia_upreg_AUC           </td><td>-0.4740875                          </td><td>Microglia                           </td><td>Zhong                               </td></tr>
	<tr><td>Neftel_NPC2_AUC                     </td><td> 0.4038164                          </td><td>NPC2-GBM                            </td><td>Neftel                              </td></tr>
	<tr><td>nowakowski_IPC.div1_upreg_AUC       </td><td> 0.4123759                          </td><td>IPC-dividing-1                      </td><td>Nowakowski                          </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>GeneSignature</th><th scope=col>SpearmanCorrleation</th><th scope=col>NewName</th><th scope=col>Source</th></tr></thead>
<tbody>
	<tr><th scope=row>15</th><td>nowakowski_IPC.nEN3_upreg_AUC     </td><td>0.4768571                         </td><td>IPC-Newborn-Excitatory-Neurons-3  </td><td>Nowakowski                        </td></tr>
	<tr><th scope=row>16</th><td>nowakowski_nEN.early1_upreg_AUC   </td><td>0.4779848                         </td><td>Newborn-Excitory-Neurons-Early-1  </td><td>Nowakowski                        </td></tr>
	<tr><th scope=row>17</th><td>nowakowski_MGE.RG2_upreg_AUC      </td><td>0.4799810                         </td><td>MGE-Radial-Glia-2                 </td><td>Nowakowski                        </td></tr>
	<tr><th scope=row>18</th><td>nowakowski_MGE.IPC3_upreg_AUC     </td><td>0.5707849                         </td><td>MGE-IPC-3                         </td><td>Nowakowski                        </td></tr>
	<tr><th scope=row>19</th><td>VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC</td><td>0.5847428                         </td><td>Proneural-TCGA-GBM                </td><td>Verhaak                           </td></tr>
	<tr><th scope=row>20</th><td>Neftel_NPC1_AUC                   </td><td>0.6253201                         </td><td>NPC1-GBM                          </td><td>Neftel                            </td></tr>
</tbody>
</table>





                    Cahoy    InHouse     Neftel Nowakowski    Verhaak      Zhong
             0          1          1          4         10          2          2



```R
cols <- c()
cols[grep("Zhong", bb$Source)] <- "#1b9e77"
cols[grep("Neftel", bb$Source)] <- "#d95f02"
cols[grep("Verhaak", bb$Source)] <- "#7570b3"
cols[grep("Nowakowski", bb$Source)] <- "#666666"
cols[grep("Cahoy", bb$Source)] <- "white"
cols[grep("Liddelow", bb$Source)] <- "black"
cols[grep("InHouse", bb$Source)] <- "red"
cols

pdf("~/Desktop/celltype.GSC_Tumour_PC2.pdf", width = 10, height = 7)

par(mar=c(4.1, 13, 4.1, 11), xpd=TRUE)

barplot(c(as.numeric(bb$SpearmanCorrleation)),
        horiz = T,
        names.arg = as.character(bb$NewName),
        las= 1,
        xlim = c(-1,1),
        cex.names = 0.8,
         col = cols,
         xlab = "Spearman Correlation Coefficient"
       )


dev.off()
```


<ol class=list-inline>
	<li>'red'</li>
	<li>'#7570b3'</li>
	<li>'#d95f02'</li>
	<li>'#1b9e77'</li>
	<li>'#d95f02'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#1b9e77'</li>
	<li>'#d95f02'</li>
	<li>'#666666'</li>
	<li>'white'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#666666'</li>
	<li>'#7570b3'</li>
	<li>'#d95f02'</li>
</ol>




<strong>pdf:</strong> 2


----
# Make a correlation heatmap
---


```R
library(ComplexHeatmap)
```


```R
cols <- c(grep("^PC2$", colnames(dat)),
          grep("^PC1$", colnames(dat)),
          grep("AUC", colnames(dat))
         )

corr.input <- dat[ ,cols]
#head(corr.input)

res <- cor(corr.input, method = "spearman")
#head(res)
```


```R
aa <- data.frame(res)
head(aa)


```


<table>
<thead><tr><th></th><th scope=col>PC2</th><th scope=col>PC1</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>Zhong_Interneurons_upreg_AUC</th><th scope=col>Zhong_OPC_upreg_AUC</th><th scope=col>Zhong_Astrocytes_upreg_AUC</th><th scope=col>Zhong_Microglia_upreg_AUC</th><th scope=col>nowakowski_EN.PFC1_upreg_AUC</th><th scope=col>nowakowski_nEN.early2_upreg_AUC</th><th scope=col>⋯</th><th scope=col>scGBM_AUC</th><th scope=col>scBTSC_AUC</th><th scope=col>Neftel_MES2_AUC</th><th scope=col>Neftel_MES1_AUC</th><th scope=col>Neftel_AC_AUC</th><th scope=col>Neftel_OPC_AUC</th><th scope=col>Neftel_NPC1_AUC</th><th scope=col>Neftel_NPC2_AUC</th><th scope=col>Neftel_G1.S_AUC</th><th scope=col>Neftel_G2.M_AUC</th></tr></thead>
<tbody>
	<tr><th scope=row>PC2</th><td> 1.00000000 </td><td>-0.1085506  </td><td> 0.2093883  </td><td>-0.13376075 </td><td>-0.01738562 </td><td> 0.44471525 </td><td> 0.14932025 </td><td>-0.47408748 </td><td>-0.02824179 </td><td> 0.08951598 </td><td>⋯           </td><td>-0.0745432  </td><td> 0.06682802 </td><td>-0.28303674 </td><td>-0.47731322 </td><td> 0.24562830 </td><td> 0.44508207 </td><td> 0.625320066</td><td> 0.40381637 </td><td> 0.2763517  </td><td> 0.09423859 </td></tr>
	<tr><th scope=row>PC1</th><td>-0.10855057 </td><td> 1.0000000  </td><td> 0.4983853  </td><td> 0.36986103 </td><td>-0.12975246 </td><td>-0.74411213 </td><td>-0.77291537 </td><td>-0.36157099 </td><td>-0.30206710 </td><td>-0.37825622 </td><td>⋯           </td><td>-0.8411720  </td><td> 0.75854522 </td><td>-0.02833718 </td><td> 0.24172158 </td><td>-0.75376880 </td><td>-0.63608532 </td><td>-0.325771740</td><td>-0.08988760 </td><td> 0.4834310  </td><td> 0.26136534 </td></tr>
	<tr><th scope=row>Zhong_NPCs_upreg_AUC</th><td> 0.20938826 </td><td> 0.4983853  </td><td> 1.0000000  </td><td> 0.22048840 </td><td>-0.23803546 </td><td>-0.33537723 </td><td>-0.43299481 </td><td>-0.45438377 </td><td>-0.48795054 </td><td>-0.39888933 </td><td>⋯           </td><td>-0.5210878  </td><td> 0.51459872 </td><td>-0.11562009 </td><td> 0.05465760 </td><td>-0.38038458 </td><td>-0.27765216 </td><td>-0.025789121</td><td> 0.05057851 </td><td> 0.8036037  </td><td> 0.78894135 </td></tr>
	<tr><th scope=row>Zhong_Excitatory_neurons_upreg_AUC</th><td>-0.13376075 </td><td> 0.3698610  </td><td> 0.2204884  </td><td> 1.00000000 </td><td>-0.08540876 </td><td>-0.27070412 </td><td>-0.20192353 </td><td>-0.06543234 </td><td>-0.12905225 </td><td>-0.09155833 </td><td>⋯           </td><td>-0.2835018  </td><td> 0.29812598 </td><td> 0.05741554 </td><td> 0.22253075 </td><td>-0.18803759 </td><td>-0.22762345 </td><td>-0.146055999</td><td>-0.01717608 </td><td> 0.1684457  </td><td> 0.16130178 </td></tr>
	<tr><th scope=row>Zhong_Interneurons_upreg_AUC</th><td>-0.01738562 </td><td>-0.1297525  </td><td>-0.2380355  </td><td>-0.08540876 </td><td> 1.00000000 </td><td> 0.06418086 </td><td> 0.08649256 </td><td> 0.20638923 </td><td> 0.24300068 </td><td> 0.21377368 </td><td>⋯           </td><td> 0.1410802  </td><td>-0.10977878 </td><td> 0.11066005 </td><td> 0.04525624 </td><td> 0.03417076 </td><td> 0.03023951 </td><td> 0.009141287</td><td> 0.10175236 </td><td>-0.2055333  </td><td>-0.18295230 </td></tr>
	<tr><th scope=row>Zhong_OPC_upreg_AUC</th><td> 0.44471525 </td><td>-0.7441121  </td><td>-0.3353772  </td><td>-0.27070412 </td><td> 0.06418086 </td><td> 1.00000000 </td><td> 0.85419379 </td><td> 0.22423523 </td><td> 0.24585158 </td><td> 0.26800855 </td><td>⋯           </td><td> 0.6593317  </td><td>-0.47276092 </td><td>-0.06809916 </td><td>-0.24873545 </td><td> 0.89106440 </td><td> 0.89151011 </td><td> 0.524084399</td><td> 0.19227213 </td><td>-0.2801959  </td><td>-0.25135977 </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'PC1'</li>
	<li>'Zhong_NPCs_upreg_AUC'</li>
	<li>'Zhong_Excitatory_neurons_upreg_AUC'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Astrocytes_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
	<li>'nowakowski_EN.PFC1_upreg_AUC'</li>
	<li>'nowakowski_nEN.early2_upreg_AUC'</li>
	<li>'nowakowski_nEN.early1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div2_upreg_AUC'</li>
	<li>'nowakowski_U4_upreg_AUC'</li>
	<li>'nowakowski_U3_upreg_AUC'</li>
	<li>'nowakowski_RG.early_upreg_AUC'</li>
	<li>'nowakowski_U2_upreg_AUC'</li>
	<li>'nowakowski_oRG_upreg_AUC'</li>
	<li>'nowakowski_tRG_upreg_AUC'</li>
	<li>'nowakowski_vRG_upreg_AUC'</li>
	<li>'nowakowski_OPC_upreg_AUC'</li>
	<li>'nowakowski_Astrocyte_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG1_upreg_AUC'</li>
	<li>'nowakowski_MGE.div_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC1_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC2_upreg_AUC'</li>
	<li>'cahoy_astrocyte_AUC'</li>
	<li>'cahoy_astro_young_AUC'</li>
	<li>'cahoy_astro_mature_AUC'</li>
	<li>'cahoy_OPC_AUC'</li>
	<li>'cahoy_OL_myel_AUC'</li>
	<li>'cahoy_MOG_pos_AUC'</li>
	<li>'cahoy_astro_in_vivo_AUC'</li>
	<li>'cahoy_cultured_astroglia_AUC'</li>
	<li>'mizrak_MicrogliaA_AUC'</li>
	<li>'mizrak_MicrogliaB_1_AUC'</li>
	<li>'mizrak_Endothelial2_AUC'</li>
	<li>'mizrak_vSMC_AUC'</li>
	<li>'mizrak_Neuron1_AUC'</li>
	<li>'mizrak_Neuron3_AUC'</li>
	<li>'mizrak_Astro2_AUC'</li>
	<li>'mizrak_Astro3_AUC'</li>
	<li>'mizrak_Astro5_AUC'</li>
	<li>'mizrak_aNSC1_AUC'</li>
	<li>'mizrak_aNSC2_TAC_AUC'</li>
	<li>'mizrak_F_astro_septal_upreg_AUC'</li>
	<li>'mizrak_septal_astrocytes_F_upreg_AUC'</li>
	<li>'RNA.GSC.c1_AUC'</li>
	<li>'RNA.GSC.c2_AUC'</li>
	<li>'glioma.stem.cell_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC'</li>
	<li>'scGBM_AUC'</li>
	<li>'scBTSC_AUC'</li>
	<li>'Neftel_AC_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_G1.S_AUC'</li>
	<li>'PC2'</li>
	<li>'Zhong_OPC_upreg_AUC'</li>
	<li>'Zhong_Microglia_upreg_AUC'</li>
	<li>'nowakowski_nEN.early1_upreg_AUC'</li>
	<li>'nowakowski_EN.V1.3_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.CGE1_upreg_AUC'</li>
	<li>'nowakowski_IN.CTX.MGE1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN1_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN2_upreg_AUC'</li>
	<li>'nowakowski_IPC.div1_upreg_AUC'</li>
	<li>'nowakowski_IPC.div2_upreg_AUC'</li>
	<li>'nowakowski_IPC.nEN3_upreg_AUC'</li>
	<li>'nowakowski_nIN1_upreg_AUC'</li>
	<li>'nowakowski_nIN2_upreg_AUC'</li>
	<li>'nowakowski_nIN5_upreg_AUC'</li>
	<li>'nowakowski_Endothelial_upreg_AUC'</li>
	<li>'nowakowski_U4_upreg_AUC'</li>
	<li>'nowakowski_Microglia_upreg_AUC'</li>
	<li>'nowakowski_OPC_upreg_AUC'</li>
	<li>'nowakowski_RG.div2_upreg_AUC'</li>
	<li>'nowakowski_MGE.RG2_upreg_AUC'</li>
	<li>'nowakowski_MGE.div_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC1_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC2_upreg_AUC'</li>
	<li>'nowakowski_MGE.IPC3_upreg_AUC'</li>
	<li>'cahoy_astro_young_AUC'</li>
	<li>'cahoy_astro_mature_AUC'</li>
	<li>'cahoy_OPC_AUC'</li>
	<li>'a1.astro_AUC'</li>
	<li>'a2.astro_AUC'</li>
	<li>'mizrak_Neuron2_AUC'</li>
	<li>'mizrak_Astro2_AUC'</li>
	<li>'mizrak_Ependymal_AUC'</li>
	<li>'mizrak_aNSC1_AUC'</li>
	<li>'mizrak_NB_AUC'</li>
	<li>'mizrak_lateral_astro_M_upreg_AUC'</li>
	<li>'mizrak_M_astro_lateral_upreg_AUC'</li>
	<li>'mizrak_lateral_astro_F_upreg_AUC'</li>
	<li>'RNA.GSC.c1_AUC'</li>
	<li>'RNA.GSC.c2_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC'</li>
	<li>'VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC'</li>
	<li>'Suva_Stemness_AUC'</li>
	<li>'Neftel_MES1_AUC'</li>
	<li>'Neftel_OPC_AUC'</li>
	<li>'Neftel_NPC1_AUC'</li>
	<li>'Neftel_NPC2_AUC'</li>
</ol>




```R
pc1.sigs <- rownames(aa[abs(aa$PC1) >= 0.35, ])
#pc1.sigs

pc1 <- aa["PC1", pc1.sigs]
pc1 <- pc1[ ,-1]
head(pc1)


pc2.sigs <- rownames(aa[abs(aa$PC2) >= 0.35, ])
#pc2.sigs

pc2 <- aa["PC2", pc2.sigs]
pc2 <- pc2[ ,-19]
head(pc2)
```


<table>
<thead><tr><th></th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>Zhong_OPC_upreg_AUC</th><th scope=col>Zhong_Astrocytes_upreg_AUC</th><th scope=col>Zhong_Microglia_upreg_AUC</th><th scope=col>nowakowski_nEN.early2_upreg_AUC</th><th scope=col>nowakowski_IPC.div1_upreg_AUC</th><th scope=col>nowakowski_IPC.div2_upreg_AUC</th><th scope=col>nowakowski_U4_upreg_AUC</th><th scope=col>nowakowski_U3_upreg_AUC</th><th scope=col>⋯</th><th scope=col>nowakowski_MGE.IPC2_upreg_AUC.1</th><th scope=col>cahoy_astro_young_AUC.1</th><th scope=col>cahoy_astro_mature_AUC.1</th><th scope=col>cahoy_OPC_AUC.1</th><th scope=col>mizrak_aNSC1_AUC.1</th><th scope=col>RNA.GSC.c1_AUC.1</th><th scope=col>RNA.GSC.c2_AUC.1</th><th scope=col>VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC.1</th><th scope=col>VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC.1</th><th scope=col>Neftel_OPC_AUC.1</th></tr></thead>
<tbody>
	<tr><th scope=row>PC1</th><td>0.4983853 </td><td>0.369861  </td><td>-0.7441121</td><td>-0.7729154</td><td>-0.361571 </td><td>-0.3782562</td><td>0.4674822 </td><td>0.5392047 </td><td>-0.4416253</td><td>0.4303661 </td><td>⋯         </td><td>0.4293216 </td><td>0.3759063 </td><td>-0.5299804</td><td>-0.3712721</td><td>0.5017999 </td><td>-0.808716 </td><td>0.633446  </td><td>-0.4956153</td><td>0.4388963 </td><td>-0.6360853</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Zhong_OPC_upreg_AUC</th><th scope=col>Zhong_Microglia_upreg_AUC</th><th scope=col>nowakowski_nEN.early1_upreg_AUC</th><th scope=col>nowakowski_IPC.nEN1_upreg_AUC</th><th scope=col>nowakowski_IPC.div1_upreg_AUC</th><th scope=col>nowakowski_IPC.div2_upreg_AUC</th><th scope=col>nowakowski_OPC_upreg_AUC</th><th scope=col>nowakowski_MGE.div_upreg_AUC</th><th scope=col>nowakowski_MGE.IPC2_upreg_AUC</th><th scope=col>cahoy_astro_young_AUC</th><th scope=col>⋯</th><th scope=col>mizrak_NB_AUC</th><th scope=col>mizrak_lateral_astro_M_upreg_AUC</th><th scope=col>RNA.GSC.c1_AUC.1</th><th scope=col>RNA.GSC.c2_AUC.1</th><th scope=col>VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC.1</th><th scope=col>VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC.1</th><th scope=col>Neftel_MES1_AUC</th><th scope=col>Neftel_OPC_AUC.1</th><th scope=col>Neftel_NPC1_AUC.1</th><th scope=col>Neftel_NPC2_AUC</th></tr></thead>
<tbody>
	<tr><th scope=row>PC2</th><td>0.4447152 </td><td>-0.4740875</td><td>0.4779848 </td><td>0.4291823 </td><td>0.4123759 </td><td>0.4165157 </td><td>0.4550529 </td><td>0.4189413 </td><td>0.3980771 </td><td>0.455536  </td><td>⋯         </td><td>0.5222532 </td><td>0.4178228 </td><td>0.3850271 </td><td>-0.5559915</td><td>0.5847428 </td><td>-0.4804564</td><td>-0.4773132</td><td>0.4450821 </td><td>0.6253201 </td><td>0.4038164 </td></tr>
</tbody>
</table>




```R
grep("^PC2", colnames(pc2))
```


19



```R
#pdf("~/Desktop/CorrMatrix.pdf", width = 14, height = 14)

Heatmap(as.matrix(aa),
       border = TRUE,
       row_names_gp = gpar(fontsize = 7),
       column_names_gp = gpar(fontsize = 7),
        column_km = 2,
        #column_km=3
       )

#dev.off()
```




<strong>pdf:</strong> 2



```R
pdf("~/Desktop/CorrMatrix.pdf", width = 14, height = 4)


Heatmap(as.matrix((pc1)),
       border = TRUE,
       row_names_gp = gpar(fontsize = 7),
       column_names_gp = gpar(fontsize = 7),
        #row_km = 2,
        column_km=2
       )

Heatmap(as.matrix((pc2)),
       border = TRUE,
       row_names_gp = gpar(fontsize = 7),
       column_names_gp = gpar(fontsize = 7),
        column_km = 2,
        #column_km=3
       )

dev.off()
```






<strong>pdf:</strong> 2



```R
Heatmap(as.matrix((pc1)),
       border = TRUE,
       row_names_gp = gpar(fontsize = 7),
       column_names_gp = gpar(fontsize = 7),
        row_km = 2,
        #column_km=3
       )
```




![png](output_21_1.png)



```R
cat(c(colnames(pc1), colnames(pc2)), sep = "\n")
```

    Zhong_NPCs_upreg_AUC
    Zhong_Excitatory_neurons_upreg_AUC
    Zhong_OPC_upreg_AUC
    Zhong_Astrocytes_upreg_AUC
    Zhong_Microglia_upreg_AUC
    nowakowski_nEN.early2_upreg_AUC
    nowakowski_IPC.div1_upreg_AUC
    nowakowski_IPC.div2_upreg_AUC
    nowakowski_U4_upreg_AUC
    nowakowski_U3_upreg_AUC
    nowakowski_RG.early_upreg_AUC
    nowakowski_U2_upreg_AUC
    nowakowski_oRG_upreg_AUC
    nowakowski_tRG_upreg_AUC
    nowakowski_vRG_upreg_AUC
    nowakowski_OPC_upreg_AUC
    nowakowski_Astrocyte_upreg_AUC
    nowakowski_MGE.RG1_upreg_AUC
    nowakowski_MGE.div_upreg_AUC
    nowakowski_MGE.IPC1_upreg_AUC
    nowakowski_MGE.IPC2_upreg_AUC
    cahoy_astrocyte_AUC
    cahoy_astro_young_AUC
    cahoy_astro_mature_AUC
    cahoy_OPC_AUC
    cahoy_OL_myel_AUC
    cahoy_MOG_pos_AUC
    cahoy_astro_in_vivo_AUC
    cahoy_cultured_astroglia_AUC
    mizrak_MicrogliaB_1_AUC
    mizrak_Endothelial2_AUC
    mizrak_vSMC_AUC
    mizrak_Neuron1_AUC
    mizrak_Neuron3_AUC
    mizrak_Astro3_AUC
    mizrak_Astro5_AUC
    mizrak_aNSC1_AUC
    mizrak_aNSC2_TAC_AUC
    mizrak_F_astro_septal_upreg_AUC
    mizrak_septal_astrocytes_F_upreg_AUC
    RNA.GSC.c1_AUC
    RNA.GSC.c2_AUC
    VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC
    VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC
    VERHAAK_GLIOBLASTOMA_CLASSICAL_AUC
    scGBM_AUC
    scBTSC_AUC
    Neftel_AC_AUC
    Neftel_OPC_AUC
    Neftel_G1.S_AUC
    Zhong_OPC_upreg_AUC.1
    Zhong_Microglia_upreg_AUC.1
    nowakowski_IPC.div1_upreg_AUC.1
    nowakowski_IPC.div2_upreg_AUC.1
    nowakowski_U4_upreg_AUC.1
    nowakowski_OPC_upreg_AUC.1
    nowakowski_MGE.div_upreg_AUC.1
    nowakowski_MGE.IPC1_upreg_AUC.1
    nowakowski_MGE.IPC2_upreg_AUC.1
    cahoy_astro_young_AUC.1
    cahoy_astro_mature_AUC.1
    cahoy_OPC_AUC.1
    mizrak_aNSC1_AUC.1
    RNA.GSC.c1_AUC.1
    RNA.GSC.c2_AUC.1
    VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC.1
    VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC.1
    Neftel_OPC_AUC.1
    Zhong_OPC_upreg_AUC
    Zhong_Microglia_upreg_AUC
    nowakowski_nEN.early1_upreg_AUC
    nowakowski_IPC.nEN1_upreg_AUC
    nowakowski_IPC.div1_upreg_AUC
    nowakowski_IPC.div2_upreg_AUC
    nowakowski_OPC_upreg_AUC
    nowakowski_MGE.div_upreg_AUC
    nowakowski_MGE.IPC2_upreg_AUC
    cahoy_astro_young_AUC
    mizrak_Astro2_AUC
    mizrak_aNSC1_AUC
    RNA.GSC.c1_AUC
    RNA.GSC.c2_AUC
    VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC
    VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC
    Neftel_OPC_AUC
    Neftel_NPC1_AUC
    Zhong_OPC_upreg_AUC.1
    Zhong_Microglia_upreg_AUC.1
    nowakowski_nEN.early1_upreg_AUC.1
    nowakowski_IN.CTX.CGE1_upreg_AUC
    nowakowski_IN.CTX.MGE1_upreg_AUC
    nowakowski_IPC.nEN1_upreg_AUC.1
    nowakowski_IPC.nEN2_upreg_AUC
    nowakowski_IPC.div1_upreg_AUC.1
    nowakowski_IPC.div2_upreg_AUC.1
    nowakowski_IPC.nEN3_upreg_AUC
    nowakowski_nIN5_upreg_AUC
    nowakowski_OPC_upreg_AUC.1
    nowakowski_RG.div2_upreg_AUC
    nowakowski_MGE.RG2_upreg_AUC
    nowakowski_MGE.div_upreg_AUC.1
    nowakowski_MGE.IPC2_upreg_AUC.1
    nowakowski_MGE.IPC3_upreg_AUC
    cahoy_astro_young_AUC.1
    a2.astro_AUC
    mizrak_Astro2_AUC.1
    mizrak_Ependymal_AUC
    mizrak_aNSC1_AUC.1
    mizrak_NB_AUC
    mizrak_lateral_astro_M_upreg_AUC
    RNA.GSC.c1_AUC.1
    RNA.GSC.c2_AUC.1
    VERHAAK_GLIOBLASTOMA_PRONEURAL_AUC.1
    VERHAAK_GLIOBLASTOMA_MESENCHYMAL_AUC.1
    Neftel_MES1_AUC
    Neftel_OPC_AUC.1
    Neftel_NPC1_AUC.1
    Neftel_NPC2_AUC
