
---
# Replot GSC PCAs for paper
---
L.Richards  
Plot a hexbin for the PCA plot with hexbin with counts


```R
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(gtable)
library(ks)
library(Seurat)
library(dplyr)
library(gridExtra)
```

---
# 1.0 Plot PCA of all GSCs
---


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/PCA/BTSC_PCATop100_noribo_CCregressed.Rdata")
```


```R
pc <- data.frame(PCA_global@cell.embeddings)
head(pc)
```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>PC9</th><th scope=col>PC10</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>11.341987   </td><td> -6.298533  </td><td>  7.696922  </td><td>-14.333232  </td><td>  3.1420173 </td><td> 5.668781   </td><td>-20.125500  </td><td>-14.9932983 </td><td>37.7201701  </td><td>  1.96138971</td><td>⋯           </td><td>-1.8393812  </td><td> 3.01262004 </td><td>-1.5421039  </td><td>-5.535492   </td><td> 0.3064857  </td><td> 3.4879667  </td><td>-1.1006038  </td><td> 5.48002890 </td><td> 2.3850766  </td><td>-0.7314307  </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td> 1.056595   </td><td>-15.348043  </td><td>  3.393379  </td><td>-16.688545  </td><td>-14.6037247 </td><td>15.168300   </td><td> -2.028318  </td><td> -2.7877237 </td><td> 8.0774740  </td><td>  7.48866938</td><td>⋯           </td><td> 2.3678906  </td><td>-0.25562920 </td><td>-1.0171882  </td><td> 1.765621   </td><td>-0.8290217  </td><td>-3.4059357  </td><td> 0.1872662  </td><td> 3.45681279 </td><td> 2.2525218  </td><td>-4.5820123  </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td> 6.251153   </td><td>-13.456937  </td><td> -1.247351  </td><td>  4.091611  </td><td>  6.8636853 </td><td> 3.218648   </td><td>  1.631251  </td><td>  0.1006166 </td><td>-1.2595184  </td><td> -0.29313224</td><td>⋯           </td><td>-2.2848058  </td><td> 0.53486980 </td><td> 1.7745577  </td><td> 2.541347   </td><td> 0.6057934  </td><td>-3.2096123  </td><td>-3.0717438  </td><td> 0.07519628 </td><td> 1.0784422  </td><td> 4.4650540  </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>22.724968   </td><td> -1.382396  </td><td> -3.078221  </td><td> -5.785007  </td><td>  0.8187603 </td><td>-9.864684   </td><td> -3.553908  </td><td>  5.0713773 </td><td>-0.5594375  </td><td>-11.52830568</td><td>⋯           </td><td> 0.4969856  </td><td>-3.62822516 </td><td> 0.5136418  </td><td>-2.503812   </td><td> 0.4321562  </td><td> 3.4160284  </td><td>-2.2778279  </td><td> 1.13719371 </td><td>-0.6722062  </td><td>-2.1051435  </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td> 5.451771   </td><td>-16.402341  </td><td>-13.218978  </td><td> 14.898474  </td><td> 12.6460271 </td><td>-1.169721   </td><td> -8.060404  </td><td> -5.0854755 </td><td> 1.8670333  </td><td>  5.19176806</td><td>⋯           </td><td> 0.5075883  </td><td> 0.36602465 </td><td> 2.3488703  </td><td>-0.699607   </td><td> 3.0447717  </td><td> 0.5759174  </td><td>-2.3062207  </td><td> 1.37261172 </td><td> 1.2012238  </td><td> 3.4705618  </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 9.851638   </td><td>  2.187515  </td><td>  5.636207  </td><td>-24.340537  </td><td>-10.3310061 </td><td> 1.636632   </td><td>-12.572445  </td><td>-14.0158842 </td><td>21.9702649  </td><td> -0.06540948</td><td>⋯           </td><td>-2.6919082  </td><td>-0.09718196 </td><td>-0.6177636  </td><td>-3.172306   </td><td>-0.1509536  </td><td> 0.1900219  </td><td>-1.9516371  </td><td>-0.40031201 </td><td>-0.4081784  </td><td> 2.3436508  </td></tr>
</tbody>
</table>




```R
G800 <- pc[grep("G800", rownames(pc)), ]
head(G800)
```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>PC9</th><th scope=col>PC10</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>G800_L_AAACCTGAGGTTCCTA</th><td>30.96168   </td><td>40.41155   </td><td>-19.792728 </td><td>16.03173   </td><td>-11.141831 </td><td> 2.010839  </td><td> 2.999193  </td><td>-1.5642479 </td><td>-0.2755307 </td><td>-2.5106731 </td><td>⋯          </td><td> 0.602663  </td><td>-3.9472956 </td><td> 0.03162587</td><td> 1.1536188 </td><td>-3.98186068</td><td> 3.7112783 </td><td> 4.7172965 </td><td> 0.7115465 </td><td> 0.26918028</td><td>-2.416392  </td></tr>
	<tr><th scope=row>G800_L_AAACCTGCACCAACCG</th><td>25.09465   </td><td>43.82663   </td><td>-11.659691 </td><td>20.33346   </td><td> -5.102687 </td><td>10.625527  </td><td> 1.253252  </td><td>-2.1351002 </td><td> 5.3141392 </td><td> 1.1450576 </td><td>⋯          </td><td> 2.834766  </td><td>-2.4597203 </td><td>-0.08705853</td><td>-3.4521194 </td><td>-1.99311417</td><td>-0.1411535 </td><td>-1.4443838 </td><td>-0.4561517 </td><td>-0.38227543</td><td> 2.771433  </td></tr>
	<tr><th scope=row>G800_L_AAACCTGGTCAGGACA</th><td>25.92711   </td><td>36.65315   </td><td> -7.284527 </td><td>13.54269   </td><td> -7.228114 </td><td> 8.351476  </td><td> 4.672099  </td><td> 1.9545941 </td><td> 6.2470246 </td><td> 3.3062097 </td><td>⋯          </td><td>-6.463735  </td><td>-0.8048444 </td><td>-0.72430296</td><td> 4.4364599 </td><td>-3.44774190</td><td>-0.4044743 </td><td> 0.7114410 </td><td> 1.7059131 </td><td> 0.03941389</td><td>-5.712731  </td></tr>
	<tr><th scope=row>G800_L_AAACCTGGTCGGCACT</th><td>30.39746   </td><td>38.34856   </td><td> -2.437741 </td><td>21.96402   </td><td> -2.293559 </td><td> 7.425811  </td><td>-3.189295  </td><td> 2.0868445 </td><td> 1.2470134 </td><td> 0.5137708 </td><td>⋯          </td><td> 2.224422  </td><td> 0.1466515 </td><td>-0.45660816</td><td>-0.1625266 </td><td>-0.09723736</td><td> 1.9700645 </td><td>-0.5183455 </td><td>-0.5981441 </td><td> 3.40434503</td><td> 1.291191  </td></tr>
	<tr><th scope=row>G800_L_AAACCTGGTCTCAACA</th><td>23.33285   </td><td>43.38783   </td><td>-17.173407 </td><td>15.64226   </td><td>-12.261644 </td><td> 3.448182  </td><td> 5.347407  </td><td>-0.9250732 </td><td> 3.8257089 </td><td> 1.2992732 </td><td>⋯          </td><td>-6.369107  </td><td> 0.9359698 </td><td> 0.36200189</td><td> 2.9348347 </td><td>-3.34111062</td><td>-2.1498103 </td><td> 5.7272239 </td><td>-0.3862056 </td><td>-1.10801449</td><td>-3.364719  </td></tr>
	<tr><th scope=row>G800_L_AAACCTGGTTCAGTAC</th><td>17.38672   </td><td>34.65480   </td><td>-10.191308 </td><td>23.94100   </td><td>  1.283793 </td><td>13.410420  </td><td>-1.346395  </td><td>-3.5943052 </td><td>10.1041006 </td><td> 7.1787956 </td><td>⋯          </td><td>-1.128722  </td><td> 2.5503190 </td><td> 0.81158361</td><td> 2.9838280 </td><td>-1.58995361</td><td>-2.4754763 </td><td>-2.1052033 </td><td> 3.7220453 </td><td>-0.54080748</td><td> 1.408016  </td></tr>
</tbody>
</table>




```R
pdf("~/Desktop/PCA_SeuppFig_G800gone.pdf", width = 10, height = 7)
par(mfrow=c(2,3))
#plot all GSCs

smoothScatter(pc$PC1, 
     pc$PC2, 
     cex = 0.3, 
     xlab = "PC1",
     ylab = "PC2",
     main = "All GSCs"
    )

#highlight G800_L

plot(pc$PC1, 
     pc$PC2, 
     cex = 0.5, 
     xlab = "PC1",
     ylab = "PC2",
     main = "Color G800",
     col= ifelse(grepl("G800_L", rownames(pc)), "red", "grey"),
     pch = 18
    )

## justify as an outlier


par(cex.axis = 0.7)
boxplot(Scaled.PC2 ~ Sample, 
        data = dat, 
        las = 2, 
        ylab = "Deviation from PC2 mean (Z-score)",
        col = c(rep("grey", 18), "red", rep("grey", 10)),
        cex= 0.5
       )
abline(h=2, lty = 2, col = "red")


dev.off()

```


<strong>pdf:</strong> 2



```R
dat <- scale(pc$PC2)
colnames(dat) <- "Scaled.PC2"
dat <- data.frame(dat)
dat$Sample <- gsub('.{17}$', '', rownames(pc))
dat$Above2 <- dat$Scaled.PC2 >= 2
head(dat)
```


<table>
<thead><tr><th scope=col>Scaled.PC2</th><th scope=col>Sample</th><th scope=col>Above2</th></tr></thead>
<tbody>
	<tr><td>-0.4753536</td><td>BT127_L   </td><td>FALSE     </td></tr>
	<tr><td>-1.1623714</td><td>BT127_L   </td><td>FALSE     </td></tr>
	<tr><td>-1.0188030</td><td>BT127_L   </td><td>FALSE     </td></tr>
	<tr><td>-0.1021319</td><td>BT127_L   </td><td>FALSE     </td></tr>
	<tr><td>-1.2424112</td><td>BT127_L   </td><td>FALSE     </td></tr>
	<tr><td> 0.1688875</td><td>BT127_L   </td><td>FALSE     </td></tr>
</tbody>
</table>




```R
pdf("~/Desktop/G800_L_Outlier_GSCs.pdf", width = 12, height = 5)
par(cex.axis = 0.9)
boxplot(Scaled.PC2 ~ Sample, 
        data = dat, 
        las = 2, 
        ylab = "Deviation from PC2 mean (Z-score)",
        col = c(rep("grey", 18), "red", rep("grey", 10)),
        cex= 0.5
       )
abline(h=2, lty = 2, col = "red")

counts <- prop.table(table(dat$Above2, dat$Sample), 2)
head(counts)

barplot(counts, 
        las = 2,
        ylab = "% cells >2 deviations from PC2 mean",
        col = c("grey", "red")
       )

dev.off()


pdf("~/Desktop/G800_L_Outlier_GSCs_barplot.pdf", width = 4, height = 4)
barplot(counts, 
        las = 2,
        ylab = "% cells >2 deviations from PC2 mean",
        col = c("grey", "red")
       )

dev.off()
```


           
                 BT127_L      BT147_L       BT48_L       BT67_L       BT73_L
      FALSE 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000
      TRUE  0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
           
                  BT84_L       BT89_L       BT94_L       G523_L       G549_L
      FALSE 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000
      TRUE  0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
           
                  G564_L       G566_L       G583_L       G620_L       G637_L
      FALSE 1.0000000000 0.9845938375 0.9994288978 1.0000000000 0.9940442428
      TRUE  0.0000000000 0.0154061625 0.0005711022 0.0000000000 0.0059557572
           
                  G729_L       G797_L       G799_L       G800_L       G837_L
      FALSE 0.9878514405 1.0000000000 0.9934086629 0.1112894596 1.0000000000
      TRUE  0.0121485595 0.0000000000 0.0065913371 0.8887105404 0.0000000000
           
                  G851_L       G876_L       G885_L       G895_L     G945-I_L
      FALSE 1.0000000000 1.0000000000 1.0000000000 1.0000000000 1.0000000000
      TRUE  0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
           
                G945-J_L     G945-K_L     G946-J_L     G946-K_L
      FALSE 1.0000000000 1.0000000000 1.0000000000 1.0000000000
      TRUE  0.0000000000 0.0000000000 0.0000000000 0.0000000000



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2


---
# 2.0 Plot PCA without G800_L
---


```R
pca <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_PCA.rds")
pca <- data.frame(pca@cell.embeddings)
#head(pca)

auc <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_AUC_metadata.rds")
#head(auc)

dat <- cbind(auc, pca)
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>percent.mito</th><th scope=col>S.Score</th><th scope=col>G2M.Score</th><th scope=col>Phase</th><th scope=col>CC.Difference</th><th scope=col>Zhong_NPCs_upreg_AUC</th><th scope=col>Zhong_Excitatory_neurons_upreg_AUC</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td> 640       </td><td>  875      </td><td>BT127_L    </td><td>0.043428571</td><td> 0.10188958</td><td> 0.05323467</td><td>S          </td><td> 0.04865491</td><td>0.10723880 </td><td>0.00000000 </td><td>⋯          </td><td> 1.2268702 </td><td>-6.6163387 </td><td>-0.7391893 </td><td>-0.4254762 </td><td>-1.6319558 </td><td>-2.3351249 </td><td> 0.17579454</td><td>-1.3976101 </td><td>-3.085101  </td><td>-0.2152068 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>1036       </td><td> 2408      </td><td>BT127_L    </td><td>0.002076412</td><td> 0.04110671</td><td> 0.32122382</td><td>G2M        </td><td>-0.28011711</td><td>0.19251951 </td><td>0.10632035 </td><td>⋯          </td><td> 0.1520107 </td><td> 0.3222401 </td><td>-4.8567436 </td><td> 2.4735756 </td><td> 1.1667153 </td><td>-2.3513450 </td><td>-0.98279624</td><td>-3.1559862 </td><td>-4.092167  </td><td> 1.0089448 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>3240       </td><td>10058      </td><td>BT127_L    </td><td>0.078047326</td><td> 0.19368572</td><td>-0.19847450</td><td>S          </td><td> 0.39216023</td><td>0.10578866 </td><td>0.01610390 </td><td>⋯          </td><td>-2.7361134 </td><td> 2.6413471 </td><td>-5.1265480 </td><td>-1.3310544 </td><td>-1.6752804 </td><td>-2.3433830 </td><td>-2.94998720</td><td> 1.0941308 </td><td> 1.013161  </td><td> 0.2452560 </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>3337       </td><td>10798      </td><td>BT127_L    </td><td>0.061863308</td><td>-0.13859617</td><td>-0.23440054</td><td>G1         </td><td> 0.09580437</td><td>0.06700741 </td><td>0.00000000 </td><td>⋯          </td><td> 0.5666615 </td><td>-1.5162332 </td><td> 6.7010819 </td><td>-0.7774230 </td><td>-1.5631920 </td><td> 0.5793188 </td><td>-0.02352712</td><td>-0.3816066 </td><td>-6.154046  </td><td> 3.0005149 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>4140       </td><td>14601      </td><td>BT127_L    </td><td>0.081501267</td><td> 0.37305820</td><td> 0.40968637</td><td>G2M        </td><td>-0.03662817</td><td>0.25095289 </td><td>0.01142857 </td><td>⋯          </td><td>-1.6911843 </td><td>-2.6541281 </td><td>-1.3874565 </td><td> 1.5717118 </td><td>-0.6124571 </td><td>-1.9403558 </td><td>-0.12033434</td><td> 2.5073746 </td><td> 1.524806  </td><td>-2.2212321 </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 543       </td><td>  820      </td><td>BT127_L    </td><td>0.108536585</td><td>-0.11929411</td><td>-0.14126218</td><td>G1         </td><td> 0.02196808</td><td>0.05839376 </td><td>0.00000000 </td><td>⋯          </td><td>-0.3203109 </td><td> 0.5971691 </td><td> 0.5093098 </td><td>-1.6936641 </td><td>-0.2763692 </td><td> 1.0759517 </td><td>-4.69101274</td><td> 1.1368124 </td><td> 1.580629  </td><td> 2.1385084 </td></tr>
</tbody>
</table>




```R



smoothScatter(dat$PC1, 
     dat$PC2, 
     cex = 0.3, 
     xlab = "PC1",
     ylab = "PC2",
     main = "All GSCs"
    )


```


![png](output_11_0.png)



```R
##plot hexbin for the AUC scores
x <- "PC1"
y <- "PC2"
sig <- "RNA.GSC.c1_AUC"

p <- ggplot(dat, aes_string(x=x, y=y, z=sig)) +
 stat_summary_hex(bins=100, fun = "median") +
 theme_classic(base_size=13) +
 scale_fill_gradientn("Median Raw \nAUC Score", colours = c("white", brewer.pal(n = 8, name = "YlOrRd"))) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(z = "Raw AUC") + #ggtitle(sig) +
theme(legend.position = "none") 
p


sig <- "RNA.GSC.c2_AUC"

q <- ggplot(dat, aes_string(x=x, y=y, z=sig)) +
 stat_summary_hex(bins=100, fun = "median") +
 theme_classic(base_size=13) +
 scale_fill_gradientn("Median Raw \nAUC Score", colours = c("white", brewer.pal(n = 8, name = "YlOrRd"))) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(z = "Raw AUC") + #ggtitle(sig) +
theme(legend.position = "none") 
q
```






![png](output_12_2.png)



![png](output_12_3.png)



```R
x <- ggplot(dat, aes(PC1, PC2)) + geom_point(size=0.4, alpha = 0.3, color = "black") +  
geom_density_2d(bins = 3, color = "grey") + theme_classic(base_size=13) +
 #scale_fill_gradient("Median Raw \nAUC Score", low = "white", high = "red") +
     #scale_fill_viridis("Median Raw \nAUC Score", begin = 0, end = 1) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) 
x
```




![png](output_13_1.png)

xdensity <- ggplot(dat, aes(PC1, fill = "grey")) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999')) + 
  theme(legend.position = "none")
xdensity


```R
pdf("~/Desktop/Figure1_GSCS_noG800_l.pdf", height = 5, width = 17)

grid.arrange(p, x, q,
             ncol=3, 
             nrow=1,
             widths=c(3, 3,3), 
             heights=c(3)
            )

dev.off()
```

---
## 3.0 Identify top correlated genes to PC1
---


```R
pca <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_PCA.rds")
pca.genes <- data.frame(pca@gene.loadings)
```


```R
nGenes <- 100

top <- rownames(pca.genes[order(pca.genes$PC1, decreasing = TRUE), ])[1:nGenes]
bottom <- rownames(pca.genes[order(pca.genes$PC1, decreasing = FALSE), ])[1:nGenes]

top
bottom

```


<ol class=list-inline>
	<li>'OCIAD2'</li>
	<li>'CAV1'</li>
	<li>'TXN'</li>
	<li>'PLP2'</li>
	<li>'PFN1'</li>
	<li>'S100A10'</li>
	<li>'PITX1'</li>
	<li>'TMEM14B'</li>
	<li>'C12orf75'</li>
	<li>'UPP1'</li>
	<li>'PRELID1'</li>
	<li>'S100A6'</li>
	<li>'CAST'</li>
	<li>'TMSB4X'</li>
	<li>'HEBP2'</li>
	<li>'SIVA1'</li>
	<li>'CDKN2A'</li>
	<li>'MDK'</li>
	<li>'CTSC'</li>
	<li>'GCHFR'</li>
	<li>'MRPL33'</li>
	<li>'PMEPA1'</li>
	<li>'ANXA2'</li>
	<li>'S100A11'</li>
	<li>'TNFRSF12A'</li>
	<li>'PHLDA2'</li>
	<li>'LAYN'</li>
	<li>'HOXB2'</li>
	<li>'GLRX5'</li>
	<li>'TAGLN2'</li>
	<li>'HMGN1'</li>
	<li>'TFDP1'</li>
	<li>'ITGA3'</li>
	<li>'POLE4'</li>
	<li>'PTRF'</li>
	<li>'RNF7'</li>
	<li>'HMGA1'</li>
	<li>'DKK1'</li>
	<li>'CDKN3'</li>
	<li>'MRPL14'</li>
	<li>'EIF1AX'</li>
	<li>'EIF5A'</li>
	<li>'DCBLD2'</li>
	<li>'LINC00707'</li>
	<li>'MIR4435-2HG'</li>
	<li>'HINT1'</li>
	<li>'CD44'</li>
	<li>'CENPW'</li>
	<li>'NRP2'</li>
	<li>'PAWR'</li>
	<li>'SYAP1'</li>
	<li>'TUBA1C'</li>
	<li>'COPRS'</li>
	<li>'CHCHD7'</li>
	<li>'DHCR24'</li>
	<li>'ENY2'</li>
	<li>'AXL'</li>
	<li>'ZDHHC12'</li>
	<li>'HTATIP2'</li>
	<li>'ARHGAP29'</li>
	<li>'CHCHD2'</li>
	<li>'KTN1'</li>
	<li>'CGNL1'</li>
	<li>'TMEM167A'</li>
	<li>'ARPC1B'</li>
	<li>'HSPB1'</li>
	<li>'CYB5R2'</li>
	<li>'H2AFZ'</li>
	<li>'MRPL41'</li>
	<li>'TMSB10'</li>
	<li>'TRMT112'</li>
	<li>'STRA13'</li>
	<li>'NHP2'</li>
	<li>'EFHD2'</li>
	<li>'CSRP2'</li>
	<li>'CTNNAL1'</li>
	<li>'CAV2'</li>
	<li>'SH3BGRL3'</li>
	<li>'LSM5'</li>
	<li>'MALT1'</li>
	<li>'FBXO32'</li>
	<li>'C9orf3'</li>
	<li>'POPDC3'</li>
	<li>'ATP5G1'</li>
	<li>'EIF4EBP1'</li>
	<li>'LSM3'</li>
	<li>'TCEB1'</li>
	<li>'TPGS2'</li>
	<li>'DEK'</li>
	<li>'NDUFS6'</li>
	<li>'PFN2'</li>
	<li>'RALBP1'</li>
	<li>'FLG'</li>
	<li>'FAM65B'</li>
	<li>'KRT7'</li>
	<li>'ATOX1'</li>
	<li>'NEDD9'</li>
	<li>'FOSL1'</li>
	<li>'PPP1R14B'</li>
	<li>'CLEC2B'</li>
</ol>




<ol class=list-inline>
	<li>'PTPRZ1'</li>
	<li>'TUBB2B'</li>
	<li>'S100B'</li>
	<li>'C1orf61'</li>
	<li>'PTN'</li>
	<li>'TTYH1'</li>
	<li>'BAALC'</li>
	<li>'TUBA1A'</li>
	<li>'FEZ1'</li>
	<li>'BCAN'</li>
	<li>'FXYD6'</li>
	<li>'MAP2'</li>
	<li>'FGFBP3'</li>
	<li>'DBI'</li>
	<li>'CCND2'</li>
	<li>'GNG7'</li>
	<li>'GPM6B'</li>
	<li>'TUBB2A'</li>
	<li>'PRDX2'</li>
	<li>'OLIG1'</li>
	<li>'ADGRG1'</li>
	<li>'KLHDC8A'</li>
	<li>'CADM4'</li>
	<li>'PMP2'</li>
	<li>'ITM2B'</li>
	<li>'KCNQ2'</li>
	<li>'FAM181B'</li>
	<li>'GLUL'</li>
	<li>'ASCL1'</li>
	<li>'METRN'</li>
	<li>'TRIM9'</li>
	<li>'NOVA1'</li>
	<li>'SCRG1'</li>
	<li>'C2orf80'</li>
	<li>'OLIG2'</li>
	<li>'PCSK1N'</li>
	<li>'BMP7'</li>
	<li>'ELAVL3'</li>
	<li>'SNTG1'</li>
	<li>'SRI'</li>
	<li>'MT3'</li>
	<li>'CIRBP'</li>
	<li>'PON2'</li>
	<li>'AGT'</li>
	<li>'NLRP1'</li>
	<li>'CLU'</li>
	<li>'PLEKHB1'</li>
	<li>'SEMA5A'</li>
	<li>'GRIA2'</li>
	<li>'TSC22D4'</li>
	<li>'MARCKSL1'</li>
	<li>'NCAM1'</li>
	<li>'SNRPN'</li>
	<li>'NME4'</li>
	<li>'TCF4'</li>
	<li>'PLTP'</li>
	<li>'WSCD1'</li>
	<li>'MALAT1'</li>
	<li>'PAFAH1B3'</li>
	<li>'GLTSCR2'</li>
	<li>'BCHE'</li>
	<li>'STMN3'</li>
	<li>'CKB'</li>
	<li>'EEF2'</li>
	<li>'APP'</li>
	<li>'DPP6'</li>
	<li>'GPR37L1'</li>
	<li>'C4orf48'</li>
	<li>'AIF1L'</li>
	<li>'PLEKHA4'</li>
	<li>'PCDH9'</li>
	<li>'ZNF428'</li>
	<li>'MAPT'</li>
	<li>'TSPAN3'</li>
	<li>'SPRY1'</li>
	<li>'NAV1'</li>
	<li>'ETV1'</li>
	<li>'ATP1B2'</li>
	<li>'RFX4'</li>
	<li>'RUFY3'</li>
	<li>'RABAC1'</li>
	<li>'CSPG5'</li>
	<li>'SOX6'</li>
	<li>'APOE'</li>
	<li>'ID4'</li>
	<li>'TSPAN7'</li>
	<li>'PODXL2'</li>
	<li>'IDH1'</li>
	<li>'NDFIP1'</li>
	<li>'CPNE2'</li>
	<li>'ECH1'</li>
	<li>'SIRT2'</li>
	<li>'SOX2'</li>
	<li>'ZEB1'</li>
	<li>'PTOV1'</li>
	<li>'RNF157'</li>
	<li>'NES'</li>
	<li>'TNFRSF21'</li>
	<li>'TMEM59L'</li>
	<li>'GPRC5B'</li>
</ol>




```R
a <- pca.genes[order(pca.genes$PC1, decreasing = TRUE), ]

```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>PC9</th><th scope=col>PC10</th><th scope=col>⋯</th><th scope=col>PC91</th><th scope=col>PC92</th><th scope=col>PC93</th><th scope=col>PC94</th><th scope=col>PC95</th><th scope=col>PC96</th><th scope=col>PC97</th><th scope=col>PC98</th><th scope=col>PC99</th><th scope=col>PC100</th></tr></thead>
<tbody>
	<tr><th scope=row>OCIAD2</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CAV1</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TXN</th><td>TRUE </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>PLP2</th><td>TRUE </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>PFN1</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>S100A10</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>PITX1</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>TMEM14B</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>C12orf75</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>UPP1</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>PRELID1</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>S100A6</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CAST</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>TMSB4X</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>HEBP2</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>SIVA1</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>CDKN2A</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>MDK</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>CTSC</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>GCHFR</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>MRPL33</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>PMEPA1</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ANXA2</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>S100A11</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TNFRSF12A</th><td>TRUE </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>PHLDA2</th><td>TRUE </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>LAYN</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>HOXB2</th><td>TRUE </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>GLRX5</th><td>TRUE </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TAGLN2</th><td>TRUE </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋱</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>METRN</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>ASCL1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>GLUL</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>FAM181B</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>KCNQ2</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>ITM2B</th><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>PMP2</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>CADM4</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>KLHDC8A</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>ADGRG1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>OLIG1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>PRDX2</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TUBB2A</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>GPM6B</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>GNG7</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>CCND2</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>DBI</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>FGFBP3</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>MAP2</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>FXYD6</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>BCAN</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><th scope=row>FEZ1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TUBA1A</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>BAALC</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TTYH1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><th scope=row>PTN</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>C1orf61</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>S100B</th><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><th scope=row>TUBB2B</th><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>⋯    </td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><th scope=row>PTPRZ1</th><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td> TRUE</td><td>⋯    </td><td> TRUE</td><td>FALSE</td><td> TRUE</td><td> TRUE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
</tbody>
</table>


