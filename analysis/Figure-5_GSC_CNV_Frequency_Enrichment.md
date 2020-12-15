
---
# Enrichment of CNVs within Developmental and Injury Response GSCs
---
L.Richards  
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed

Analysis with G800_L (outlier) removed


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed")
```

-----
## 1.0 Replot frequency
----
suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/"))

#load in CNV scores for all cells
#remove G800_L
load("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/GlobalBTSC_CNVs.RData")

BTSC.CNVs <- BTSC.CNVs[ ,-grep("G800_L", colnames(BTSC.CNVs))]
dim(BTSC.CNVs)

avg.cnv.df <- avg.cnv.df[ ,-grep("G800_L", colnames(avg.cnv.df))]
dim(avg.cnv.df)

##load in meta with AUCell scores

meta <- readRDS("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_AUC_metadata.rds")
colnames(meta)


## save all the data with G800_L removed
save(BTSC.CNVs, avg.cnv.df, meta, CNV.genes, file = "GlobalBTSC_CNVs_NoG800L.Rdata")

----
### Count number of cells with gains and losses for each genes within Dev and IR bins

25292 Developmental cells, 37630 IR cells
### classify cells as C1 or C2

dat <- meta[ ,c("orig.ident", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")]
colnames(dat) <- c("Sample", "C1_AUC", "C2_AUC")

dat$Dev <- dat$C1_AUC > 0.11
dat$IR <- dat$C2_AUC > 0.2
dat$ID <- paste0(dat$Dev, dat$IR)

dat$ID <- gsub("TRUEFALSE", "Dev", dat$ID)
dat$ID <- gsub("FALSETRUE", "IR", dat$ID)
dat$ID <- gsub("FALSEFALSE", "LowLow", dat$ID)
dat$ID <- gsub("TRUETRUE", "HiHi", dat$ID)
head(dat)
#identify cells with label of interest

ID <- "IR"

gain.thres <- 0.17
loss.thres <-  -0.15


cells <- dat[dat$ID == ID ,]
dim(cells)

# subset CNV matrix by cells of interest
cnv <- BTSC.CNVs[ ,rownames(cells)]
dim(cnv)

#classify cnvs as gains if > 0.15

CNV.genes$Index <- c(1:nrow(CNV.genes))
CNV.genes$Total.Cells <- nrow(cells)

CNV.genes$Gain <- apply(cnv, 1, function(x) sum(x > gain.thres )) 
CNV.genes$Deletion <- apply(cnv, 1, function(x) sum(x < loss.thres))
    
    
save.file <- paste0(ID, "_numCellsCNV.rds")
print(save.file)
saveRDS(CNV.genes, file = save.file)
----
### Plot gains and deletions across cell types fro Dev and IR




```R
#load in C1 cell
CNV.genes <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed/Dev_numCellsCNV.rds")
```


```R
#make this into a proportion of cell

CNV.genes$Prop.Gain <- CNV.genes$Gain / CNV.genes$Total.Cells
CNV.genes$Prop.Deletion <- (CNV.genes$Deletion / CNV.genes$Total.Cells) * -1
head(CNV.genes)

gains <- CNV.genes[ ,c(2,6)]
colnames(gains)[2] <- "Prop"
dels <- CNV.genes[ ,c(2,7)]
colnames(dels)[2] <- "Prop"
plot.dat <- rbind(gains,dels)
head(plot.dat)

# add chr arm information
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/C1_C2_enrichment/CNV_gene_positionArms.Rdata")
CNV.genes <- cbind(CNV.genes, dat)

x <- as.numeric(table(CNV.genes$Chromosome)) #chromosome breaks
y <- as.numeric(table(CNV.genes$Arm)[unique(CNV.genes$Arm)]) #arm breaks

pdf("~/Desktop/Dev_CNVs.pdf", width = 11, height = 5)

plot(plot.dat$Index, 
     plot.dat$Prop,
     col = "white",
     ylim = c(-1, 1),
     xlab = "Genomic Coordinate",
     ylab = "Proportion Cells",
     main = "Developmental GSCs \n 25,292 cells"
    )
lines(plot.dat$Index, 
      plot.dat$Prop, 
      type="h", col = ifelse(plot.dat$Prop > 0, "darkred", "darkblue")
     )
abline(h=0, col = "black")
abline(v=cumsum(x), lty = 1)
abline(v=cumsum(y), lty = 2, lwd = 1)

dev.off()
```


<table>
<thead><tr><th></th><th scope=col>Chromosome</th><th scope=col>Index</th><th scope=col>Total.Cells</th><th scope=col>Gain</th><th scope=col>Deletion</th><th scope=col>Prop.Gain</th><th scope=col>Prop.Deletion</th><th scope=col>dat$Arm</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>1         </td><td>1         </td><td>25292     </td><td>2740      </td><td>4085      </td><td>0.1083347 </td><td>-0.1615135</td><td>1p        </td></tr>
	<tr><th scope=row>TMEM201</th><td>1         </td><td>2         </td><td>25292     </td><td>3112      </td><td>3816      </td><td>0.1230429 </td><td>-0.1508777</td><td>1p        </td></tr>
	<tr><th scope=row>CLSTN1</th><td>1         </td><td>3         </td><td>25292     </td><td>3181      </td><td>3864      </td><td>0.1257710 </td><td>-0.1527756</td><td>1p        </td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>1         </td><td>4         </td><td>25292     </td><td>3170      </td><td>3826      </td><td>0.1253361 </td><td>-0.1512731</td><td>1p        </td></tr>
	<tr><th scope=row>LZIC</th><td>1         </td><td>5         </td><td>25292     </td><td>3466      </td><td>3588      </td><td>0.1370394 </td><td>-0.1418630</td><td>1p        </td></tr>
	<tr><th scope=row>NMNAT1</th><td>1         </td><td>6         </td><td>25292     </td><td>3565      </td><td>3302      </td><td>0.1409537 </td><td>-0.1305551</td><td>1p        </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Index</th><th scope=col>Prop</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>1        </td><td>0.1083347</td></tr>
	<tr><th scope=row>TMEM201</th><td>2        </td><td>0.1230429</td></tr>
	<tr><th scope=row>CLSTN1</th><td>3        </td><td>0.1257710</td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>4        </td><td>0.1253361</td></tr>
	<tr><th scope=row>LZIC</th><td>5        </td><td>0.1370394</td></tr>
	<tr><th scope=row>NMNAT1</th><td>6        </td><td>0.1409537</td></tr>
</tbody>
</table>




<strong>pdf:</strong> 2



```R
#load in IR cell
CNV.genes <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed/IR_numCellsCNV.rds")

# add chr arm information
#load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/C1_C2_enrichment/CNV_gene_positionArms.Rdata")

#CNV.genes <- cbind(CNV.genes, dat)
#head(CNV.genes)

#make this into a proportion of cells

CNV.genes$Prop.Gain <- CNV.genes$Gain / CNV.genes$Total.Cells
CNV.genes$Prop.Deletion <- (CNV.genes$Deletion / CNV.genes$Total.Cells) * -1
head(CNV.genes)

gains <- CNV.genes[ ,c(2,6)]
colnames(gains)[2] <- "Prop"
dels <- CNV.genes[ ,c(2,7)]
colnames(dels)[2] <- "Prop"
plot.dat <- rbind(gains,dels)
head(plot.dat)

# add chr arm information
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/C1_C2_enrichment/CNV_gene_positionArms.Rdata")
CNV.genes <- cbind(CNV.genes, dat)

x <- as.numeric(table(CNV.genes$Chromosome)) #chromosome breaks
y <- as.numeric(table(CNV.genes$Arm)[unique(CNV.genes$Arm)]) #arm breaks

pdf("~/Desktop/IR_CNVs.pdf", width = 11, height = 5)

plot(plot.dat$Index, 
     plot.dat$Prop,
     col = "white",
     ylim = c(-1, 1),
     xlab = "Genomic Coordinate",
     ylab = "Proportion Cells",
     main = "Injury Response GSCs \n 37,630 cells"
    )
lines(plot.dat$Index, 
      plot.dat$Prop, 
      type="h", col = ifelse(plot.dat$Prop > 0, "darkred", "darkblue")
     )
abline(h=0, col = "black")
abline(v=cumsum(x), lty = 1)
abline(v=cumsum(y), lty = 2, lwd = 1)

dev.off()
```


<table>
<thead><tr><th></th><th scope=col>Chromosome</th><th scope=col>Index</th><th scope=col>Total.Cells</th><th scope=col>Gain</th><th scope=col>Deletion</th><th scope=col>Prop.Gain</th><th scope=col>Prop.Deletion</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>1          </td><td>1          </td><td>37630      </td><td> 8528      </td><td>3474       </td><td>0.2266277  </td><td>-0.09231996</td></tr>
	<tr><th scope=row>TMEM201</th><td>1          </td><td>2          </td><td>37630      </td><td> 9163      </td><td>3151       </td><td>0.2435025  </td><td>-0.08373638</td></tr>
	<tr><th scope=row>CLSTN1</th><td>1          </td><td>3          </td><td>37630      </td><td> 9189      </td><td>3151       </td><td>0.2441935  </td><td>-0.08373638</td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>1          </td><td>4          </td><td>37630      </td><td> 9338      </td><td>3079       </td><td>0.2481531  </td><td>-0.08182301</td></tr>
	<tr><th scope=row>LZIC</th><td>1          </td><td>5          </td><td>37630      </td><td>10289      </td><td>2784       </td><td>0.2734255  </td><td>-0.07398352</td></tr>
	<tr><th scope=row>NMNAT1</th><td>1          </td><td>6          </td><td>37630      </td><td>10099      </td><td>2636       </td><td>0.2683763  </td><td>-0.07005049</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Index</th><th scope=col>Prop</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>1        </td><td>0.2266277</td></tr>
	<tr><th scope=row>TMEM201</th><td>2        </td><td>0.2435025</td></tr>
	<tr><th scope=row>CLSTN1</th><td>3        </td><td>0.2441935</td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>4        </td><td>0.2481531</td></tr>
	<tr><th scope=row>LZIC</th><td>5        </td><td>0.2734255</td></tr>
	<tr><th scope=row>NMNAT1</th><td>6        </td><td>0.2683763</td></tr>
</tbody>
</table>




<strong>pdf:</strong> 2


----
## 2.0 Redo CNV enrichment between groups
----

---
## Plot ratio of gains:losses between two subgroups
---

Anlaysis Dir: 
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed

Analysis Plan:   
> classify each gene to a chr arm     
> Average CNV profile across charm for  each cell    
> classify each chr arm/cell as either a gain or loss    
> sum number of proneural and immunomesenchymal cells with a gain or loss at each arm - turn into proportion of cells per subtype
> calculate gian to loss  ratio within each cell type  
> Plot compariso of ratios and use [cooks.distance](https://stats.stackexchange.com/questions/164099/removing-outliers-based-on-cooks-distance-in-r-language) to identify outliers  
 


```R
#load chr arm average CNVs
#remove G800_L
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/C1_C2_enrichment/BTSCs_armAveragedCNVs.RData")
arm.CNVs <- arm.CNVs[-grep("G800_L", rownames(arm.CNVs)) ,] 
dim(arm.CNVs)
head(arm.CNVs)
```


<ol class=list-inline>
	<li>65655</li>
	<li>31</li>
</ol>




<table>
<thead><tr><th></th><th scope=col>1p</th><th scope=col>1q</th><th scope=col>2p</th><th scope=col>2q</th><th scope=col>3p</th><th scope=col>3q</th><th scope=col>4p</th><th scope=col>4q</th><th scope=col>5p</th><th scope=col>5q</th><th scope=col>⋯</th><th scope=col>13p</th><th scope=col>14p</th><th scope=col>15p</th><th scope=col>16p</th><th scope=col>17p</th><th scope=col>18p</th><th scope=col>19p</th><th scope=col>20p</th><th scope=col>21p</th><th scope=col>22p</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>-0.020358796 </td><td>-0.0132178000</td><td>-0.01378511  </td><td>-0.010178565 </td><td>-0.05862751  </td><td> 0.004459093 </td><td> 0.003333794 </td><td> 0.00000000  </td><td> 0.004251096 </td><td>-0.04950182  </td><td>⋯            </td><td>-0.006278459 </td><td> 0.034453499 </td><td>-0.002637757 </td><td>-0.001086704 </td><td>-0.01594457  </td><td> 0.07429173  </td><td>-0.07974874  </td><td>0.04353808   </td><td>-0.02543742  </td><td>-0.066083338 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>-0.073001708 </td><td>-0.0127253707</td><td> 0.06763824  </td><td>-0.005244394 </td><td>-0.04459320  </td><td>-0.022521376 </td><td>-0.072142442 </td><td> 0.00000000  </td><td> 0.055161456 </td><td>-0.12344028  </td><td>⋯            </td><td>-0.114629521 </td><td>-0.096130580 </td><td> 0.145074116 </td><td>-0.008182345 </td><td>-0.01591449  </td><td>-0.11907045  </td><td>-0.01350640  </td><td>0.08002627   </td><td>-0.19986996  </td><td> 0.055798528 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>-0.012922645 </td><td>-0.1089236711</td><td> 0.05801856  </td><td> 0.145588374 </td><td>-0.14788876  </td><td>-0.264993256 </td><td>-0.114694127 </td><td>-0.17847101  </td><td>-0.010662029 </td><td>-0.09495720  </td><td>⋯            </td><td>-0.288225028 </td><td>-0.115265774 </td><td>-0.002031741 </td><td>-0.044850600 </td><td>-0.03275135  </td><td> 0.00000000  </td><td>-0.00248200  </td><td>0.18312963   </td><td>-0.12726369  </td><td> 0.076537236 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>-0.048447659 </td><td> 0.0009027557</td><td>-0.01901104  </td><td>-0.058429976 </td><td>-0.11040428  </td><td> 0.082372006 </td><td>-0.099069501 </td><td>-0.03547551  </td><td>-0.016390138 </td><td> 0.03289960  </td><td>⋯            </td><td>-0.369210928 </td><td>-0.169119966 </td><td>-0.053774621 </td><td>-0.084354206 </td><td> 0.02378564  </td><td>-0.15613218  </td><td>-0.02831435  </td><td>0.13764076   </td><td> 0.18304052  </td><td>-0.060691440 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td> 0.009780757 </td><td>-0.1083831052</td><td>-0.06960953  </td><td>-0.002279644 </td><td>-0.13749249  </td><td>-0.171111547 </td><td>-0.066089129 </td><td>-0.17825962  </td><td>-0.093245210 </td><td>-0.01540262  </td><td>⋯            </td><td>-0.254078210 </td><td> 0.011841104 </td><td>-0.085081740 </td><td>-0.015204833 </td><td> 0.02730335  </td><td> 0.18135757  </td><td> 0.07664405  </td><td>0.31305420   </td><td>-0.03039282  </td><td> 0.007116073 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 0.003512711 </td><td>-0.0605612505</td><td>-0.01212824  </td><td> 0.033361497 </td><td>-0.04341334  </td><td> 0.017562216 </td><td> 0.021160570 </td><td>-0.04960070  </td><td> 0.061341793 </td><td>-0.02360285  </td><td>⋯            </td><td>-0.001941482 </td><td> 0.005718923 </td><td> 0.005417582 </td><td>-0.006585626 </td><td> 0.05556598  </td><td>-0.11939263  </td><td>-0.12870088  </td><td>0.06139763   </td><td> 0.04286166  </td><td>-0.101062145 </td></tr>
</tbody>
</table>




```R
#load up classification of Dev or IR

meta <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_AUC_metadata.rds")


### classify cells as C1 or C2

dat <- meta[ ,c("orig.ident", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")]
colnames(dat) <- c("Sample", "C1_AUC", "C2_AUC")

dat$Dev <- dat$C1_AUC > 0.11
dat$IR <- dat$C2_AUC > 0.2
dat$ID <- paste0(dat$Dev, dat$IR)

dat$ID <- gsub("TRUEFALSE", "Dev", dat$ID)
dat$ID <- gsub("FALSETRUE", "IR", dat$ID)
dat$ID <- gsub("FALSEFALSE", "LowLow", dat$ID)
dat$ID <- gsub("TRUETRUE", "HiHi", dat$ID)
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>Sample</th><th scope=col>C1_AUC</th><th scope=col>C2_AUC</th><th scope=col>Dev</th><th scope=col>IR</th><th scope=col>ID</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>BT127_L  </td><td>0.1883243</td><td>0.1665197</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>BT127_L  </td><td>0.1258778</td><td>0.1515399</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>BT127_L  </td><td>0.1398942</td><td>0.1547603</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>BT127_L  </td><td>0.1518339</td><td>0.1674444</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>BT127_L  </td><td>0.1335071</td><td>0.1654191</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>BT127_L  </td><td>0.1652775</td><td>0.2021541</td><td>TRUE     </td><td> TRUE    </td><td>HiHi     </td></tr>
</tbody>
</table>




```R
#### count across arms which are gained and lost for DE 

gain.thres <- 0.17
loss.thres <- -0.15

table(dat$ID)
```


    
       Dev   HiHi     IR LowLow 
     25292   1106  37630   1627 



```R
###############################################
### Developmental Cells
### 4

ID <- "Dev"
cells <- dat[dat$ID == ID ,]
dim(cells)

Dev_gains <- data.frame(apply(arm.CNVs[rownames(cells), ], 2, function(x) sum(x > 0.17)))
colnames(Dev_gains) <- "Dev_Gains"

Dev_loss <- data.frame(apply(arm.CNVs[rownames(cells), ], 2, function(x) sum(x < -0.15)))
colnames(Dev_loss) <- "Dev_Deletions"

    
Dev <- cbind(Dev_gains, Dev_loss)
Dev_counts <- Dev
Dev_counts$Dev_Total <- nrow(cells)
head(Dev_counts)    

Dev <- Dev/nrow(cells)
head(Dev)
```


<ol class=list-inline>
	<li>25292</li>
	<li>6</li>
</ol>




<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th><th scope=col>Dev_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>  23 </td><td> 155 </td><td>25292</td></tr>
	<tr><th scope=row>1q</th><td>2172 </td><td> 711 </td><td>25292</td></tr>
	<tr><th scope=row>2p</th><td>  33 </td><td>  78 </td><td>25292</td></tr>
	<tr><th scope=row>2q</th><td>1775 </td><td> 396 </td><td>25292</td></tr>
	<tr><th scope=row>3p</th><td> 163 </td><td>1977 </td><td>25292</td></tr>
	<tr><th scope=row>3q</th><td>1950 </td><td>1382 </td><td>25292</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>0.0009093785</td><td>0.006128420 </td></tr>
	<tr><th scope=row>1q</th><td>0.0858769571</td><td>0.028111656 </td></tr>
	<tr><th scope=row>2p</th><td>0.0013047604</td><td>0.003083979 </td></tr>
	<tr><th scope=row>2q</th><td>0.0701802942</td><td>0.015657125 </td></tr>
	<tr><th scope=row>3p</th><td>0.0064447256</td><td>0.078167009 </td></tr>
	<tr><th scope=row>3q</th><td>0.0770994781</td><td>0.054641784 </td></tr>
</tbody>
</table>




```R
###############################################
### Injury Response Cells
### 

ID <- "IR"
cells <- dat[dat$ID == ID ,]
dim(cells)

IR_gains <- data.frame(apply(arm.CNVs[rownames(cells), ], 2, function(x) sum(x > 0.17)))
colnames(IR_gains) <- "IR_Gains"

IR_loss <- data.frame(apply(arm.CNVs[rownames(cells), ], 2, function(x) sum(x < -0.15)))
colnames(IR_loss) <- "IR_Deletions"
    head(Dev)
    
IR <- cbind(IR_gains, IR_loss)
IR_counts <- IR
IR_counts$IR_Total <- nrow(cells)
head(IR_counts)     

IR <- IR/nrow(cells)
head(IR)
```


<ol class=list-inline>
	<li>37630</li>
	<li>6</li>
</ol>




<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>0.0009093785</td><td>0.006128420 </td></tr>
	<tr><th scope=row>1q</th><td>0.0858769571</td><td>0.028111656 </td></tr>
	<tr><th scope=row>2p</th><td>0.0013047604</td><td>0.003083979 </td></tr>
	<tr><th scope=row>2q</th><td>0.0701802942</td><td>0.015657125 </td></tr>
	<tr><th scope=row>3p</th><td>0.0064447256</td><td>0.078167009 </td></tr>
	<tr><th scope=row>3q</th><td>0.0770994781</td><td>0.054641784 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td> 512 </td><td> 535 </td><td>37630</td></tr>
	<tr><th scope=row>1q</th><td> 741 </td><td>4425 </td><td>37630</td></tr>
	<tr><th scope=row>2p</th><td>1357 </td><td> 453 </td><td>37630</td></tr>
	<tr><th scope=row>2q</th><td>5206 </td><td>3391 </td><td>37630</td></tr>
	<tr><th scope=row>3p</th><td> 488 </td><td>2325 </td><td>37630</td></tr>
	<tr><th scope=row>3q</th><td>5116 </td><td>1213 </td><td>37630</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>0.01360617</td><td>0.01421738</td></tr>
	<tr><th scope=row>1q</th><td>0.01969174</td><td>0.11759235</td></tr>
	<tr><th scope=row>2p</th><td>0.03606165</td><td>0.01203827</td></tr>
	<tr><th scope=row>2q</th><td>0.13834706</td><td>0.09011427</td></tr>
	<tr><th scope=row>3p</th><td>0.01296838</td><td>0.06178581</td></tr>
	<tr><th scope=row>3q</th><td>0.13595535</td><td>0.03223492</td></tr>
</tbody>
</table>




```R
################
#bind the two together

df <- cbind(Dev, IR)
df <- df*100
df$Dev_GL <- round(df$Dev_Gains / df$Dev_Deletions, 2)
df$IR_GL <- round(df$IR_Gains / df$IR_Deletions,2 )
df
    
save(df, file = "ChrArmAvg_C1_C2_counts.Rdata")
```


<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th><th scope=col>Dev_GL</th><th scope=col>IR_GL</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td> 0.09093785</td><td> 0.61284201</td><td> 1.3606165 </td><td> 1.4217380 </td><td>   0.15    </td><td> 0.96      </td></tr>
	<tr><th scope=row>1q</th><td> 8.58769571</td><td> 2.81116559</td><td> 1.9691735 </td><td>11.7592347 </td><td>   3.05    </td><td> 0.17      </td></tr>
	<tr><th scope=row>2p</th><td> 0.13047604</td><td> 0.30839791</td><td> 3.6061653 </td><td> 1.2038267 </td><td>   0.42    </td><td> 3.00      </td></tr>
	<tr><th scope=row>2q</th><td> 7.01802942</td><td> 1.56571248</td><td>13.8347064 </td><td> 9.0114271 </td><td>   4.48    </td><td> 1.54      </td></tr>
	<tr><th scope=row>3p</th><td> 0.64447256</td><td> 7.81670093</td><td> 1.2968376 </td><td> 6.1785809 </td><td>   0.08    </td><td> 0.21      </td></tr>
	<tr><th scope=row>3q</th><td> 7.70994781</td><td> 5.46417840</td><td>13.5955355 </td><td> 3.2234919 </td><td>   1.41    </td><td> 4.22      </td></tr>
	<tr><th scope=row>4p</th><td> 0.14629132</td><td> 9.69476514</td><td> 3.5184693 </td><td> 7.1698113 </td><td>   0.02    </td><td> 0.49      </td></tr>
	<tr><th scope=row>4q</th><td> 2.38810691</td><td>15.38035743</td><td>16.3566303 </td><td>12.4608026 </td><td>   0.16    </td><td> 1.31      </td></tr>
	<tr><th scope=row>5p</th><td> 0.24513680</td><td> 1.81875692</td><td> 5.0332182 </td><td> 3.9728940 </td><td>   0.13    </td><td> 1.27      </td></tr>
	<tr><th scope=row>5q</th><td> 2.27344615</td><td> 7.65064052</td><td>18.9476482 </td><td> 9.6146691 </td><td>   0.30    </td><td> 1.97      </td></tr>
	<tr><th scope=row>6p</th><td> 0.08303021</td><td> 6.95872213</td><td> 0.8557002 </td><td> 3.6247675 </td><td>   0.01    </td><td> 0.24      </td></tr>
	<tr><th scope=row>6q</th><td> 1.22963783</td><td>27.96931836</td><td>23.3483922 </td><td> 8.4480468 </td><td>   0.04    </td><td> 2.76      </td></tr>
	<tr><th scope=row>7p</th><td> 9.87664084</td><td> 0.05930729</td><td>24.3715121 </td><td> 0.7175126 </td><td> 166.53    </td><td>33.97      </td></tr>
	<tr><th scope=row>7q</th><td>79.19500237</td><td> 0.07116875</td><td>64.1376561 </td><td> 0.6882806 </td><td>1112.78    </td><td>93.19      </td></tr>
	<tr><th scope=row>8p</th><td> 5.95445200</td><td> 2.18646212</td><td> 7.2149880 </td><td> 5.4530959 </td><td>   2.72    </td><td> 1.32      </td></tr>
	<tr><th scope=row>8q</th><td>10.98766408</td><td> 7.87205440</td><td>19.5322881 </td><td> 8.9821951 </td><td>   1.40    </td><td> 2.17      </td></tr>
	<tr><th scope=row>9p</th><td>14.68448521</td><td> 1.30476040</td><td> 0.9992028 </td><td>30.0345469 </td><td>  11.25    </td><td> 0.03      </td></tr>
	<tr><th scope=row>9q</th><td> 5.90700617</td><td>14.32468765</td><td>16.8828063 </td><td>11.2091416 </td><td>   0.41    </td><td> 1.51      </td></tr>
	<tr><th scope=row>10p</th><td> 0.47050451</td><td>42.40471295</td><td> 3.1836301 </td><td>37.6773851 </td><td>   0.01    </td><td> 0.08      </td></tr>
	<tr><th scope=row>11p</th><td> 0.43492013</td><td> 5.00553535</td><td> 6.6276907 </td><td> 2.8381610 </td><td>   0.09    </td><td> 2.34      </td></tr>
	<tr><th scope=row>12p</th><td> 2.73999684</td><td> 0.68005693</td><td> 9.0911507 </td><td> 2.0595270 </td><td>   4.03    </td><td> 4.41      </td></tr>
	<tr><th scope=row>13p</th><td> 1.25731457</td><td>39.30491855</td><td> 1.8602179 </td><td>53.3005581 </td><td>   0.03    </td><td> 0.03      </td></tr>
	<tr><th scope=row>14p</th><td> 1.32452950</td><td>14.43934841</td><td> 6.9306404 </td><td>11.3101249 </td><td>   0.09    </td><td> 0.61      </td></tr>
	<tr><th scope=row>15p</th><td> 1.53012810</td><td> 7.19199747</td><td> 4.6266277 </td><td>13.9516343 </td><td>   0.21    </td><td> 0.33      </td></tr>
	<tr><th scope=row>16p</th><td> 0.83820971</td><td> 5.29021034</td><td> 0.8238108 </td><td>27.2761095 </td><td>   0.16    </td><td> 0.03      </td></tr>
	<tr><th scope=row>17p</th><td> 4.68132216</td><td> 2.54230587</td><td>14.3422801 </td><td>10.2923200 </td><td>   1.84    </td><td> 1.39      </td></tr>
	<tr><th scope=row>18p</th><td> 5.45627076</td><td>15.03242132</td><td>23.1411108 </td><td>26.5825140 </td><td>   0.36    </td><td> 0.87      </td></tr>
	<tr><th scope=row>19p</th><td>35.90858770</td><td> 2.12715483</td><td> 3.0002657 </td><td>46.2024980 </td><td>  16.88    </td><td> 0.06      </td></tr>
	<tr><th scope=row>20p</th><td>35.10200854</td><td> 0.90542464</td><td>22.8620781 </td><td> 2.3651342 </td><td>  38.77    </td><td> 9.67      </td></tr>
	<tr><th scope=row>21p</th><td> 4.19104855</td><td>10.72671200</td><td>17.1724688 </td><td>14.4831252 </td><td>   0.39    </td><td> 1.19      </td></tr>
	<tr><th scope=row>22p</th><td> 9.24402973</td><td>20.59939902</td><td>18.3922402 </td><td>25.0119585 </td><td>   0.45    </td><td> 0.74      </td></tr>
</tbody>
</table>




```R
df <- cbind(Dev, IR)
head(df)
```


<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>0.0009093785</td><td>0.006128420 </td><td>0.01360617  </td><td>0.01421738  </td></tr>
	<tr><th scope=row>1q</th><td>0.0858769571</td><td>0.028111656 </td><td>0.01969174  </td><td>0.11759235  </td></tr>
	<tr><th scope=row>2p</th><td>0.0013047604</td><td>0.003083979 </td><td>0.03606165  </td><td>0.01203827  </td></tr>
	<tr><th scope=row>2q</th><td>0.0701802942</td><td>0.015657125 </td><td>0.13834706  </td><td>0.09011427  </td></tr>
	<tr><th scope=row>3p</th><td>0.0064447256</td><td>0.078167009 </td><td>0.01296838  </td><td>0.06178581  </td></tr>
	<tr><th scope=row>3q</th><td>0.0770994781</td><td>0.054641784 </td><td>0.13595535  </td><td>0.03223492  </td></tr>
</tbody>
</table>



----
## 3.0 Chi squared test
----

USe this to determine significnace of difference in cells with loss/gain


```R
counts <- cbind(Dev_counts, IR_counts)
saveRDS(counts, file = "GSC_noG900_CNVcounts_Dev_IR.rds")
```


```R
head(counts)

#set up a contingency table for each of the chr arms
```


<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th><th scope=col>Dev_Total</th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>  23 </td><td> 155 </td><td>25292</td><td> 512 </td><td> 535 </td><td>37630</td></tr>
	<tr><th scope=row>1q</th><td>2172 </td><td> 711 </td><td>25292</td><td> 741 </td><td>4425 </td><td>37630</td></tr>
	<tr><th scope=row>2p</th><td>  33 </td><td>  78 </td><td>25292</td><td>1357 </td><td> 453 </td><td>37630</td></tr>
	<tr><th scope=row>2q</th><td>1775 </td><td> 396 </td><td>25292</td><td>5206 </td><td>3391 </td><td>37630</td></tr>
	<tr><th scope=row>3p</th><td> 163 </td><td>1977 </td><td>25292</td><td> 488 </td><td>2325 </td><td>37630</td></tr>
	<tr><th scope=row>3q</th><td>1950 </td><td>1382 </td><td>25292</td><td>5116 </td><td>1213 </td><td>37630</td></tr>
</tbody>
</table>




```R
gains <- counts[ ,c(grep("Gains", colnames(counts)),c(grep("Total", colnames(counts))))]        
rownames(gains) <- paste0(rownames(gains), "_gain")
colnames(gains) <- c("Dev_Count", "IR_Count", "Dev_Total", "IR_Total")
head(gains)


del <- counts[ ,c(grep("Del", colnames(counts)),c(grep("Total", colnames(counts))))]        
rownames(del) <- paste0(rownames(del), "_del")
colnames(del) <- c("Dev_Count", "IR_Count", "Dev_Total", "IR_Total")
head(del)

df <- rbind(gains, del)
head(df)
tail(df)
```


<table>
<thead><tr><th></th><th scope=col>Dev_Count</th><th scope=col>IR_Count</th><th scope=col>Dev_Total</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p_gain</th><td>  23 </td><td> 512 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>1q_gain</th><td>2172 </td><td> 741 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2p_gain</th><td>  33 </td><td>1357 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2q_gain</th><td>1775 </td><td>5206 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3p_gain</th><td> 163 </td><td> 488 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3q_gain</th><td>1950 </td><td>5116 </td><td>25292</td><td>37630</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Dev_Count</th><th scope=col>IR_Count</th><th scope=col>Dev_Total</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p_del</th><td> 155 </td><td> 535 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>1q_del</th><td> 711 </td><td>4425 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2p_del</th><td>  78 </td><td> 453 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2q_del</th><td> 396 </td><td>3391 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3p_del</th><td>1977 </td><td>2325 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3q_del</th><td>1382 </td><td>1213 </td><td>25292</td><td>37630</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Dev_Count</th><th scope=col>IR_Count</th><th scope=col>Dev_Total</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p_gain</th><td>  23 </td><td> 512 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>1q_gain</th><td>2172 </td><td> 741 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2p_gain</th><td>  33 </td><td>1357 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>2q_gain</th><td>1775 </td><td>5206 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3p_gain</th><td> 163 </td><td> 488 </td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>3q_gain</th><td>1950 </td><td>5116 </td><td>25292</td><td>37630</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Dev_Count</th><th scope=col>IR_Count</th><th scope=col>Dev_Total</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>17p_del</th><td> 643 </td><td> 3873</td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>18p_del</th><td>3802 </td><td>10003</td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>19p_del</th><td> 538 </td><td>17386</td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>20p_del</th><td> 229 </td><td>  890</td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>21p_del</th><td>2713 </td><td> 5450</td><td>25292</td><td>37630</td></tr>
	<tr><th scope=row>22p_del</th><td>5210 </td><td> 9412</td><td>25292</td><td>37630</td></tr>
</tbody>
</table>




```R
##set up contnge matrix

chrs <- rownames(df)
chrs

i <- 15
```


<ol class=list-inline>
	<li>'1p_gain'</li>
	<li>'1q_gain'</li>
	<li>'2p_gain'</li>
	<li>'2q_gain'</li>
	<li>'3p_gain'</li>
	<li>'3q_gain'</li>
	<li>'4p_gain'</li>
	<li>'4q_gain'</li>
	<li>'5p_gain'</li>
	<li>'5q_gain'</li>
	<li>'6p_gain'</li>
	<li>'6q_gain'</li>
	<li>'7p_gain'</li>
	<li>'7q_gain'</li>
	<li>'8p_gain'</li>
	<li>'8q_gain'</li>
	<li>'9p_gain'</li>
	<li>'9q_gain'</li>
	<li>'10p_gain'</li>
	<li>'11p_gain'</li>
	<li>'12p_gain'</li>
	<li>'13p_gain'</li>
	<li>'14p_gain'</li>
	<li>'15p_gain'</li>
	<li>'16p_gain'</li>
	<li>'17p_gain'</li>
	<li>'18p_gain'</li>
	<li>'19p_gain'</li>
	<li>'20p_gain'</li>
	<li>'21p_gain'</li>
	<li>'22p_gain'</li>
	<li>'1p_del'</li>
	<li>'1q_del'</li>
	<li>'2p_del'</li>
	<li>'2q_del'</li>
	<li>'3p_del'</li>
	<li>'3q_del'</li>
	<li>'4p_del'</li>
	<li>'4q_del'</li>
	<li>'5p_del'</li>
	<li>'5q_del'</li>
	<li>'6p_del'</li>
	<li>'6q_del'</li>
	<li>'7p_del'</li>
	<li>'7q_del'</li>
	<li>'8p_del'</li>
	<li>'8q_del'</li>
	<li>'9p_del'</li>
	<li>'9q_del'</li>
	<li>'10p_del'</li>
	<li>'11p_del'</li>
	<li>'12p_del'</li>
	<li>'13p_del'</li>
	<li>'14p_del'</li>
	<li>'15p_del'</li>
	<li>'16p_del'</li>
	<li>'17p_del'</li>
	<li>'18p_del'</li>
	<li>'19p_del'</li>
	<li>'20p_del'</li>
	<li>'21p_del'</li>
	<li>'22p_del'</li>
</ol>




```R
mat <- matrix(c(df[i, "Dev_Count"],
df[i, "IR_Count"],
df[i, "Dev_Total"] - df[i, "Dev_Count"],
df[i, "IR_Total"] - df[i, "IR_Count"]
  ),
              ncol = 2, byrow = TRUE
)

colnames(mat) <- c("Dev", "IR")
rownames(mat) <- c(chrs[i], paste0("Not_",chrs[i]))
mat
```


<table>
<thead><tr><th></th><th scope=col>Dev</th><th scope=col>IR</th></tr></thead>
<tbody>
	<tr><th scope=row>8p_gain</th><td> 1506</td><td> 2715</td></tr>
	<tr><th scope=row>Not_8p_gain</th><td>23786</td><td>34915</td></tr>
</tbody>
</table>




```R
chisq.test(t(mat))
prop.test(t(mat), conf.level = 0.99)
```


    
    	Pearson's Chi-squared test with Yates' continuity correction
    
    data:  t(mat)
    X-squared = 38.202, df = 1, p-value = 6.378e-10



---
## 4.0 Plot chr arm signal between bins
---




```R
library(ggplot2)
library(reshape)
library(ggpubr)
library("effsize")
```


```R
#load chr arm average CNVs
#remove G800_L
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/C1_C2_enrichment/BTSCs_armAveragedCNVs.RData")
arm.CNVs <- arm.CNVs[-grep("G800_L", rownames(arm.CNVs)) ,] 
colnames(arm.CNVs) <- paste0("chr", colnames(arm.CNVs))
dim(arm.CNVs)
head(arm.CNVs)
```


<ol class=list-inline>
	<li>65655</li>
	<li>31</li>
</ol>




<table>
<thead><tr><th></th><th scope=col>chr1p</th><th scope=col>chr1q</th><th scope=col>chr2p</th><th scope=col>chr2q</th><th scope=col>chr3p</th><th scope=col>chr3q</th><th scope=col>chr4p</th><th scope=col>chr4q</th><th scope=col>chr5p</th><th scope=col>chr5q</th><th scope=col>⋯</th><th scope=col>chr13p</th><th scope=col>chr14p</th><th scope=col>chr15p</th><th scope=col>chr16p</th><th scope=col>chr17p</th><th scope=col>chr18p</th><th scope=col>chr19p</th><th scope=col>chr20p</th><th scope=col>chr21p</th><th scope=col>chr22p</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>-0.020358796 </td><td>-0.0132178000</td><td>-0.01378511  </td><td>-0.010178565 </td><td>-0.05862751  </td><td> 0.004459093 </td><td> 0.003333794 </td><td> 0.00000000  </td><td> 0.004251096 </td><td>-0.04950182  </td><td>⋯            </td><td>-0.006278459 </td><td> 0.034453499 </td><td>-0.002637757 </td><td>-0.001086704 </td><td>-0.01594457  </td><td> 0.07429173  </td><td>-0.07974874  </td><td>0.04353808   </td><td>-0.02543742  </td><td>-0.066083338 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>-0.073001708 </td><td>-0.0127253707</td><td> 0.06763824  </td><td>-0.005244394 </td><td>-0.04459320  </td><td>-0.022521376 </td><td>-0.072142442 </td><td> 0.00000000  </td><td> 0.055161456 </td><td>-0.12344028  </td><td>⋯            </td><td>-0.114629521 </td><td>-0.096130580 </td><td> 0.145074116 </td><td>-0.008182345 </td><td>-0.01591449  </td><td>-0.11907045  </td><td>-0.01350640  </td><td>0.08002627   </td><td>-0.19986996  </td><td> 0.055798528 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>-0.012922645 </td><td>-0.1089236711</td><td> 0.05801856  </td><td> 0.145588374 </td><td>-0.14788876  </td><td>-0.264993256 </td><td>-0.114694127 </td><td>-0.17847101  </td><td>-0.010662029 </td><td>-0.09495720  </td><td>⋯            </td><td>-0.288225028 </td><td>-0.115265774 </td><td>-0.002031741 </td><td>-0.044850600 </td><td>-0.03275135  </td><td> 0.00000000  </td><td>-0.00248200  </td><td>0.18312963   </td><td>-0.12726369  </td><td> 0.076537236 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>-0.048447659 </td><td> 0.0009027557</td><td>-0.01901104  </td><td>-0.058429976 </td><td>-0.11040428  </td><td> 0.082372006 </td><td>-0.099069501 </td><td>-0.03547551  </td><td>-0.016390138 </td><td> 0.03289960  </td><td>⋯            </td><td>-0.369210928 </td><td>-0.169119966 </td><td>-0.053774621 </td><td>-0.084354206 </td><td> 0.02378564  </td><td>-0.15613218  </td><td>-0.02831435  </td><td>0.13764076   </td><td> 0.18304052  </td><td>-0.060691440 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td> 0.009780757 </td><td>-0.1083831052</td><td>-0.06960953  </td><td>-0.002279644 </td><td>-0.13749249  </td><td>-0.171111547 </td><td>-0.066089129 </td><td>-0.17825962  </td><td>-0.093245210 </td><td>-0.01540262  </td><td>⋯            </td><td>-0.254078210 </td><td> 0.011841104 </td><td>-0.085081740 </td><td>-0.015204833 </td><td> 0.02730335  </td><td> 0.18135757  </td><td> 0.07664405  </td><td>0.31305420   </td><td>-0.03039282  </td><td> 0.007116073 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 0.003512711 </td><td>-0.0605612505</td><td>-0.01212824  </td><td> 0.033361497 </td><td>-0.04341334  </td><td> 0.017562216 </td><td> 0.021160570 </td><td>-0.04960070  </td><td> 0.061341793 </td><td>-0.02360285  </td><td>⋯            </td><td>-0.001941482 </td><td> 0.005718923 </td><td> 0.005417582 </td><td>-0.006585626 </td><td> 0.05556598  </td><td>-0.11939263  </td><td>-0.12870088  </td><td>0.06139763   </td><td> 0.04286166  </td><td>-0.101062145 </td></tr>
</tbody>
</table>




```R
meta <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_AUC_metadata.rds")


### classify cells as C1 or C2

dat <- meta[ ,c("orig.ident", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")]
colnames(dat) <- c("Sample", "C1_AUC", "C2_AUC")

dat$Dev <- dat$C1_AUC > 0.11
dat$IR <- dat$C2_AUC > 0.2
dat$ID <- paste0(dat$Dev, dat$IR)

dat$ID <- gsub("TRUEFALSE", "Dev", dat$ID)
dat$ID <- gsub("FALSETRUE", "IR", dat$ID)
dat$ID <- gsub("FALSEFALSE", "LowLow", dat$ID)
dat$ID <- gsub("TRUETRUE", "HiHi", dat$ID)
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>Sample</th><th scope=col>C1_AUC</th><th scope=col>C2_AUC</th><th scope=col>Dev</th><th scope=col>IR</th><th scope=col>ID</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>BT127_L  </td><td>0.1883243</td><td>0.1665197</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>BT127_L  </td><td>0.1258778</td><td>0.1515399</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>BT127_L  </td><td>0.1398942</td><td>0.1547603</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>BT127_L  </td><td>0.1518339</td><td>0.1674444</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>BT127_L  </td><td>0.1335071</td><td>0.1654191</td><td>TRUE     </td><td>FALSE    </td><td>Dev      </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>BT127_L  </td><td>0.1652775</td><td>0.2021541</td><td>TRUE     </td><td> TRUE    </td><td>HiHi     </td></tr>
</tbody>
</table>




```R

```


<table>
<thead><tr><th></th><th scope=col>chr1p</th><th scope=col>chr1q</th><th scope=col>chr2p</th><th scope=col>chr2q</th><th scope=col>chr3p</th><th scope=col>chr3q</th><th scope=col>chr4p</th><th scope=col>chr4q</th><th scope=col>chr5p</th><th scope=col>chr5q</th><th scope=col>⋯</th><th scope=col>chr19p</th><th scope=col>chr20p</th><th scope=col>chr21p</th><th scope=col>chr22p</th><th scope=col>Sample</th><th scope=col>C1_AUC</th><th scope=col>C2_AUC</th><th scope=col>Dev</th><th scope=col>IR</th><th scope=col>ID</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>-0.020358796 </td><td>-0.0132178000</td><td>-0.01378511  </td><td>-0.010178565 </td><td>-0.05862751  </td><td> 0.004459093 </td><td> 0.003333794 </td><td> 0.000000000 </td><td> 0.004251096 </td><td>-0.04950182  </td><td>⋯            </td><td>-0.07974874  </td><td>0.04353808   </td><td>-0.02543742  </td><td>-0.066083338 </td><td>BT127_L      </td><td>0.1883243    </td><td>0.1665197    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>-0.073001708 </td><td>-0.0127253707</td><td> 0.06763824  </td><td>-0.005244394 </td><td>-0.04459320  </td><td>-0.022521376 </td><td>-0.072142442 </td><td> 0.000000000 </td><td> 0.055161456 </td><td>-0.12344028  </td><td>⋯            </td><td>-0.01350640  </td><td>0.08002627   </td><td>-0.19986996  </td><td> 0.055798528 </td><td>BT127_L      </td><td>0.1258778    </td><td>0.1515399    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>-0.012922645 </td><td>-0.1089236711</td><td> 0.05801856  </td><td> 0.145588374 </td><td>-0.14788876  </td><td>-0.264993256 </td><td>-0.114694127 </td><td>-0.178471007 </td><td>-0.010662029 </td><td>-0.09495720  </td><td>⋯            </td><td>-0.00248200  </td><td>0.18312963   </td><td>-0.12726369  </td><td> 0.076537236 </td><td>BT127_L      </td><td>0.1398942    </td><td>0.1547603    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>-0.048447659 </td><td> 0.0009027557</td><td>-0.01901104  </td><td>-0.058429976 </td><td>-0.11040428  </td><td> 0.082372006 </td><td>-0.099069501 </td><td>-0.035475512 </td><td>-0.016390138 </td><td> 0.03289960  </td><td>⋯            </td><td>-0.02831435  </td><td>0.13764076   </td><td> 0.18304052  </td><td>-0.060691440 </td><td>BT127_L      </td><td>0.1518339    </td><td>0.1674444    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td> 0.009780757 </td><td>-0.1083831052</td><td>-0.06960953  </td><td>-0.002279644 </td><td>-0.13749249  </td><td>-0.171111547 </td><td>-0.066089129 </td><td>-0.178259619 </td><td>-0.093245210 </td><td>-0.01540262  </td><td>⋯            </td><td> 0.07664405  </td><td>0.31305420   </td><td>-0.03039282  </td><td> 0.007116073 </td><td>BT127_L      </td><td>0.1335071    </td><td>0.1654191    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCATTGGTAC</th><td> 0.008919805 </td><td> 0.0422617515</td><td>-0.02302637  </td><td> 0.073810476 </td><td>-0.05333558  </td><td>-0.008270146 </td><td>-0.112756271 </td><td>-0.005821363 </td><td>-0.017424031 </td><td>-0.12394996  </td><td>⋯            </td><td>-0.06901603  </td><td>0.03032176   </td><td>-0.03680262  </td><td>-0.049434151 </td><td>BT127_L      </td><td>0.1588045    </td><td>0.1918557    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
</tbody>
</table>




    
      Dev    IR 
    25292 37630 



```R
chrs <- colnames(arm.CNVs)
chrs
```


<ol class=list-inline>
	<li>'chr1p'</li>
	<li>'chr1q'</li>
	<li>'chr2p'</li>
	<li>'chr2q'</li>
	<li>'chr3p'</li>
	<li>'chr3q'</li>
	<li>'chr4p'</li>
	<li>'chr4q'</li>
	<li>'chr5p'</li>
	<li>'chr5q'</li>
	<li>'chr6p'</li>
	<li>'chr6q'</li>
	<li>'chr7p'</li>
	<li>'chr7q'</li>
	<li>'chr8p'</li>
	<li>'chr8q'</li>
	<li>'chr9p'</li>
	<li>'chr9q'</li>
	<li>'chr10p'</li>
	<li>'chr11p'</li>
	<li>'chr12p'</li>
	<li>'chr13p'</li>
	<li>'chr14p'</li>
	<li>'chr15p'</li>
	<li>'chr16p'</li>
	<li>'chr17p'</li>
	<li>'chr18p'</li>
	<li>'chr19p'</li>
	<li>'chr20p'</li>
	<li>'chr21p'</li>
	<li>'chr22p'</li>
</ol>




```R
effect.size <- c()

pdf("With_Hybrids_chrArm_comparison.pdf", width = 4, height = 4)

for(i in 1:length(chrs)){

chr <- chrs[i]
print(chr)

Dev <- boxplot.dat[grep("Dev", boxplot.dat$ID), chr]
IR <- boxplot.dat[grep("IR", boxplot.dat$ID), chr]

b <- cohen.d(Dev,
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
```

    [1] "chr1p"
    [1] "chr1q"
    [1] "chr2p"
    [1] "chr2q"
    [1] "chr3p"
    [1] "chr3q"
    [1] "chr4p"
    [1] "chr4q"
    [1] "chr5p"
    [1] "chr5q"
    [1] "chr6p"
    [1] "chr6q"
    [1] "chr7p"
    [1] "chr7q"
    [1] "chr8p"
    [1] "chr8q"
    [1] "chr9p"
    [1] "chr9q"
    [1] "chr10p"
    [1] "chr11p"
    [1] "chr12p"
    [1] "chr13p"
    [1] "chr14p"
    [1] "chr15p"
    [1] "chr16p"
    [1] "chr17p"
    [1] "chr18p"
    [1] "chr19p"
    [1] "chr20p"
    [1] "chr21p"
    [1] "chr22p"



<strong>pdf:</strong> 2



```R
##plot effect size across chrs
pdf("~/Desktop/CNV_Enrichment_EffectSize.pdf", height = 5, width = 8)

cols <-rep("grey", length(chrs))
cols[c(12,17,28)] <- "red"

names(effect.size) <- chrs
barplot(effect.size, 
        las = 2,
        ylim = c(-1,2),
        ylab = "Effect Size (Hedge's g)",
        col = cols
       )
abline(h=0.8, lty = 2, col = "black")
abline(h=-0.8,  lty = 2, col = "black")
abline(h=0)

dev.off()
```


<strong>pdf:</strong> 2



```R
###plot out ones with large effect size

pdf("~/Desktop/enrichedCNVs_largeEffectSize.pdf", width = 8, height = 3)
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
```




<strong>pdf:</strong> 2



```R

```


```R
IR <- boxplot.dat[boxplot.dat$ID == "Dev", "chr9p"]
table(IR < -0.15)
```


    
    FALSE  TRUE 
    24962   330 



```R
head(boxplot.dat)
```


<table>
<thead><tr><th></th><th scope=col>chr1p</th><th scope=col>chr1q</th><th scope=col>chr2p</th><th scope=col>chr2q</th><th scope=col>chr3p</th><th scope=col>chr3q</th><th scope=col>chr4p</th><th scope=col>chr4q</th><th scope=col>chr5p</th><th scope=col>chr5q</th><th scope=col>⋯</th><th scope=col>chr19p</th><th scope=col>chr20p</th><th scope=col>chr21p</th><th scope=col>chr22p</th><th scope=col>Sample</th><th scope=col>C1_AUC</th><th scope=col>C2_AUC</th><th scope=col>Dev</th><th scope=col>IR</th><th scope=col>ID</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>-0.020358796 </td><td>-0.0132178000</td><td>-0.01378511  </td><td>-0.010178565 </td><td>-0.05862751  </td><td> 0.004459093 </td><td> 0.003333794 </td><td> 0.000000000 </td><td> 0.004251096 </td><td>-0.04950182  </td><td>⋯            </td><td>-0.07974874  </td><td>0.04353808   </td><td>-0.02543742  </td><td>-0.066083338 </td><td>BT127_L      </td><td>0.1883243    </td><td>0.1665197    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>-0.073001708 </td><td>-0.0127253707</td><td> 0.06763824  </td><td>-0.005244394 </td><td>-0.04459320  </td><td>-0.022521376 </td><td>-0.072142442 </td><td> 0.000000000 </td><td> 0.055161456 </td><td>-0.12344028  </td><td>⋯            </td><td>-0.01350640  </td><td>0.08002627   </td><td>-0.19986996  </td><td> 0.055798528 </td><td>BT127_L      </td><td>0.1258778    </td><td>0.1515399    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>-0.012922645 </td><td>-0.1089236711</td><td> 0.05801856  </td><td> 0.145588374 </td><td>-0.14788876  </td><td>-0.264993256 </td><td>-0.114694127 </td><td>-0.178471007 </td><td>-0.010662029 </td><td>-0.09495720  </td><td>⋯            </td><td>-0.00248200  </td><td>0.18312963   </td><td>-0.12726369  </td><td> 0.076537236 </td><td>BT127_L      </td><td>0.1398942    </td><td>0.1547603    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>-0.048447659 </td><td> 0.0009027557</td><td>-0.01901104  </td><td>-0.058429976 </td><td>-0.11040428  </td><td> 0.082372006 </td><td>-0.099069501 </td><td>-0.035475512 </td><td>-0.016390138 </td><td> 0.03289960  </td><td>⋯            </td><td>-0.02831435  </td><td>0.13764076   </td><td> 0.18304052  </td><td>-0.060691440 </td><td>BT127_L      </td><td>0.1518339    </td><td>0.1674444    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td> 0.009780757 </td><td>-0.1083831052</td><td>-0.06960953  </td><td>-0.002279644 </td><td>-0.13749249  </td><td>-0.171111547 </td><td>-0.066089129 </td><td>-0.178259619 </td><td>-0.093245210 </td><td>-0.01540262  </td><td>⋯            </td><td> 0.07664405  </td><td>0.31305420   </td><td>-0.03039282  </td><td> 0.007116073 </td><td>BT127_L      </td><td>0.1335071    </td><td>0.1654191    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCATTGGTAC</th><td> 0.008919805 </td><td> 0.0422617515</td><td>-0.02302637  </td><td> 0.073810476 </td><td>-0.05333558  </td><td>-0.008270146 </td><td>-0.112756271 </td><td>-0.005821363 </td><td>-0.017424031 </td><td>-0.12394996  </td><td>⋯            </td><td>-0.06901603  </td><td>0.03032176   </td><td>-0.03680262  </td><td>-0.049434151 </td><td>BT127_L      </td><td>0.1588045    </td><td>0.1918557    </td><td>TRUE         </td><td>FALSE        </td><td>Dev          </td></tr>
</tbody>
</table>




```R
counts
```


<table>
<thead><tr><th></th><th scope=col>Dev_Gains</th><th scope=col>Dev_Deletions</th><th scope=col>Dev_Total</th><th scope=col>IR_Gains</th><th scope=col>IR_Deletions</th><th scope=col>IR_Total</th></tr></thead>
<tbody>
	<tr><th scope=row>1p</th><td>   23</td><td>  155</td><td>25292</td><td>  512</td><td>  535</td><td>37630</td></tr>
	<tr><th scope=row>1q</th><td> 2172</td><td>  711</td><td>25292</td><td>  741</td><td> 4425</td><td>37630</td></tr>
	<tr><th scope=row>2p</th><td>   33</td><td>   78</td><td>25292</td><td> 1357</td><td>  453</td><td>37630</td></tr>
	<tr><th scope=row>2q</th><td> 1775</td><td>  396</td><td>25292</td><td> 5206</td><td> 3391</td><td>37630</td></tr>
	<tr><th scope=row>3p</th><td>  163</td><td> 1977</td><td>25292</td><td>  488</td><td> 2325</td><td>37630</td></tr>
	<tr><th scope=row>3q</th><td> 1950</td><td> 1382</td><td>25292</td><td> 5116</td><td> 1213</td><td>37630</td></tr>
	<tr><th scope=row>4p</th><td>   37</td><td> 2452</td><td>25292</td><td> 1324</td><td> 2698</td><td>37630</td></tr>
	<tr><th scope=row>4q</th><td>  604</td><td> 3890</td><td>25292</td><td> 6155</td><td> 4689</td><td>37630</td></tr>
	<tr><th scope=row>5p</th><td>   62</td><td>  460</td><td>25292</td><td> 1894</td><td> 1495</td><td>37630</td></tr>
	<tr><th scope=row>5q</th><td>  575</td><td> 1935</td><td>25292</td><td> 7130</td><td> 3618</td><td>37630</td></tr>
	<tr><th scope=row>6p</th><td>   21</td><td> 1760</td><td>25292</td><td>  322</td><td> 1364</td><td>37630</td></tr>
	<tr><th scope=row>6q</th><td>  311</td><td> 7074</td><td>25292</td><td> 8786</td><td> 3179</td><td>37630</td></tr>
	<tr><th scope=row>7p</th><td> 2498</td><td>   15</td><td>25292</td><td> 9171</td><td>  270</td><td>37630</td></tr>
	<tr><th scope=row>7q</th><td>20030</td><td>   18</td><td>25292</td><td>24135</td><td>  259</td><td>37630</td></tr>
	<tr><th scope=row>8p</th><td> 1506</td><td>  553</td><td>25292</td><td> 2715</td><td> 2052</td><td>37630</td></tr>
	<tr><th scope=row>8q</th><td> 2779</td><td> 1991</td><td>25292</td><td> 7350</td><td> 3380</td><td>37630</td></tr>
	<tr><th scope=row>9p</th><td> 3714</td><td>  330</td><td>25292</td><td>  376</td><td>11302</td><td>37630</td></tr>
	<tr><th scope=row>9q</th><td> 1494</td><td> 3623</td><td>25292</td><td> 6353</td><td> 4218</td><td>37630</td></tr>
	<tr><th scope=row>10p</th><td>  119</td><td>10725</td><td>25292</td><td> 1198</td><td>14178</td><td>37630</td></tr>
	<tr><th scope=row>11p</th><td>  110</td><td> 1266</td><td>25292</td><td> 2494</td><td> 1068</td><td>37630</td></tr>
	<tr><th scope=row>12p</th><td>  693</td><td>  172</td><td>25292</td><td> 3421</td><td>  775</td><td>37630</td></tr>
	<tr><th scope=row>13p</th><td>  318</td><td> 9941</td><td>25292</td><td>  700</td><td>20057</td><td>37630</td></tr>
	<tr><th scope=row>14p</th><td>  335</td><td> 3652</td><td>25292</td><td> 2608</td><td> 4256</td><td>37630</td></tr>
	<tr><th scope=row>15p</th><td>  387</td><td> 1819</td><td>25292</td><td> 1741</td><td> 5250</td><td>37630</td></tr>
	<tr><th scope=row>16p</th><td>  212</td><td> 1338</td><td>25292</td><td>  310</td><td>10264</td><td>37630</td></tr>
	<tr><th scope=row>17p</th><td> 1184</td><td>  643</td><td>25292</td><td> 5397</td><td> 3873</td><td>37630</td></tr>
	<tr><th scope=row>18p</th><td> 1380</td><td> 3802</td><td>25292</td><td> 8708</td><td>10003</td><td>37630</td></tr>
	<tr><th scope=row>19p</th><td> 9082</td><td>  538</td><td>25292</td><td> 1129</td><td>17386</td><td>37630</td></tr>
	<tr><th scope=row>20p</th><td> 8878</td><td>  229</td><td>25292</td><td> 8603</td><td>  890</td><td>37630</td></tr>
	<tr><th scope=row>21p</th><td> 1060</td><td> 2713</td><td>25292</td><td> 6462</td><td> 5450</td><td>37630</td></tr>
	<tr><th scope=row>22p</th><td> 2338</td><td> 5210</td><td>25292</td><td> 6921</td><td> 9412</td><td>37630</td></tr>
</tbody>
</table>




```R

```
