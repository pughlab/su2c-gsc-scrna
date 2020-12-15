
----
# CRISPR Manuscript Figures
---

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/


```R
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(UpSetR)
library(Hmisc)
library(car)
library("lme4")
library(VennDiagram)
library(edgeR)
library(ggrepel)
```


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/")
```


    Error in setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/"): cannot change working directory
    Traceback:


    1. setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/")


---
## 1.0 Essential Difference between groups
---


```R
dat <- readRDS("./CRISPR_Avg_Diff_Matrix.rds")
dat <- dat[rev(order(dat$Diff_Z)), ]
dat$GeneRank <- 1:nrow(dat)
dat$Gene <- rownames(dat)
labels <- dat[c(1:10,(nrow(dat)-10):nrow(dat)), ]$Gene #top and bottom 10 genes
dat$Label <- ifelse(dat$Gene %in% labels, dat$Gene, "")
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th><th scope=col>Developmental_Average</th><th scope=col>InjuryResponse_Average</th><th scope=col>Diff</th><th scope=col>Diff_Z</th><th scope=col>FoldChange</th><th scope=col>GeneRank</th><th scope=col>Gene</th><th scope=col>Label</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>40.43901   </td><td>40.93868   </td><td>-16.8872393</td><td>30.37648   </td><td>47.26371   </td><td>7.051752   </td><td> -1.798783 </td><td>1          </td><td>ITGB1      </td><td>ITGB1      </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td>14.84797   </td><td>28.55795   </td><td> -6.9819679</td><td>37.50056   </td><td>44.48253   </td><td>6.636796   </td><td> -5.371059 </td><td>2          </td><td>PPP1R12A   </td><td>PPP1R12A   </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td>33.67403   </td><td>37.50733   </td><td> -8.8287541</td><td>32.77162   </td><td>41.60037   </td><td>6.206774   </td><td> -3.711919 </td><td>3          </td><td>SCAP       </td><td>SCAP       </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td>33.63656   </td><td>47.18051   </td><td> -0.7493837</td><td>39.18887   </td><td>39.93825   </td><td>5.958783   </td><td>-52.294793 </td><td>4          </td><td>EXOSC5     </td><td>EXOSC5     </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>37.95362   </td><td>53.55119   </td><td>  4.2116815</td><td>43.98490   </td><td>39.77322   </td><td>5.934161   </td><td> 10.443548 </td><td>5          </td><td>SUPT16H    </td><td>SUPT16H    </td></tr>
	<tr><th scope=row>GTF2A2</th><td>  5.5289992</td><td> -7.656784 </td><td> -9.0392470</td><td>  4.311062 </td><td>39.48553   </td><td>39.81692   </td><td>27.19230   </td><td>30.89674   </td><td>22.41802   </td><td>36.88538   </td><td>24.18497   </td><td> -1.7139924</td><td>31.96190   </td><td>33.67589   </td><td>5.024433   </td><td>-18.647633 </td><td>6          </td><td>GTF2A2     </td><td>GTF2A2     </td></tr>
</tbody>
</table>




```R
#pdf("~/Desktop/CRISPR_differential.pdf", width = 7, height = 7)

plotDiff <- ggplot(dat, aes(x=GeneRank, y = Diff_Z, label = Label)) + 
    geom_point(shape = 21, fill=alpha("black", 0.2), col = "black", alpha = 0.4) + 
    geom_text_repel(data = dat[c(1:10), ],
                    nudge_x = 2500 - dat$GeneRank[1:10], 
                    direction = "y", 
                    hjust = 0,
                    segment.size = 0.2,
                    color = "black",
                    segment.color = "grey50",
                    size= 2,
                    box.padding = 0.15
                   ) + 
 geom_text_repel(data = dat[(nrow(dat)-10):nrow(dat), ],
                    nudge_x = 13000 - dat$GeneRank[(nrow(dat)-10):nrow(dat)], 
                    direction = "y", 
                    hjust = 0,
                    segment.size = 0.2,
                    color = '#b2182b',
                     segment.color = "grey50",
                    size= 2,
                   box.padding = 0.15
                   ) +  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
       text = element_text(size=10)) +
    ylab("Differential Essentiality (z-score)") + xlab("Gene Rank") + 
    geom_hline(yintercept=c(2, -2), lwd= 0.4, lty = "longdash") + ylim(c=-10, 10)


plotDiff
#dev.off()
```




![png](output_5_1.png)



```R
pdf("~/Desktop/CRISPR_plotDiff.pdf", width = 6, height = 4)

plotDiff

dev.off()
```




<strong>pdf:</strong> 2


----
## 2.0 Unbiased clustering of GSCs
----

Cluster lines with 'variably essential genes' as defined by Graham.

>There are 11 screens in total here. The G583 and G549 screens were published in my Cell Reports paper but the rest are new. The G440 and G532R are outside of the SU2C project so haven't officially been assigned to a subtype by Owen.


```R
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(circlize)
```


```R
#### READ IN CRISPR DATA

var.genes <- read.table("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/CRISPR_data/TKOv3_GBM_qBF_subset_for_clustering.txt",
                        header = T
                       )
rownames(var.genes) <- var.genes$Gene
var.genes$Gene <- NULL
colnames(var.genes) <- c("G440_L",
                         "G523_L",
                         "G729_L",
                         "G532_L",
                         "G620_L",
                         "BT67_L",
                         "G361_L",
                         "G691_L",
                         "G809_L",
                         "G549_L",
                         "G583_L"
                        )
var.genes <- na.omit(var.genes)

scaled.dat <- scale(t(var.genes))
dat <- t(scaled.dat)
```


```R
#### READ IN BULK GSVA DATA

load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/bulkRNA_data/GSC_gsva.rdata")
GSVA <- data.frame(GSC.gsva[c("RNA.GSC.c1", "RNA.GSC.c2"), colnames(GSC.gsva) %in% colnames(var.genes)])
GSVA$G440_L <- 0
GSVA$G532_L <- 0


#order the same as var.genes
GSVA <- GSVA[ ,colnames(var.genes)]
```


```R
#pdf("~/Desktop/Cluster_CRISPR_VarGenes.pdf", width = 10, height = 10)

col_fun = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"))

column_ha <- HeatmapAnnotation(Developmental = as.numeric((GSVA[1, ])),
                               Injury_Response = as.numeric((GSVA[2, ])),
                               col = list(Developmental = col_fun,
                                          Injury_Response = col_fun
                                         ),
                                annotation_legend_param = list(Developmental = list(title = "GSVA"),
                                                               Injury_Response = list(title = "GSVA")
                                                              ),
                               gp = gpar(col = "black")
                              
                              )
Heatmap(dat,
        col = viridis(60),
        column_km = 2,
        #row_km = 2,
        border = T,
        name = "qBF \n(z-score)",
        top_annotation = column_ha,
        row_title = "Varaible Genes (n=1484)",
        column_title = "",
        show_row_names = FALSE
       )

#dev.off()


#pdf("~/Desktop/Cluster_CRISPR_Correlation.pdf", width = 12, height = 10)

#### correlation heatmap

corr.dat <- cor(var.genes)
head(corr.dat)

col_fun = colorRamp2(c(0, 1), c("white", "red"))

corr.map <- Heatmap(corr.dat,
        col = RColorBrewer::brewer.pal(name = "Blues", n = 9),
        column_km = 2,
        row_km = 2,
        border = T,
        name = "Pearson \nCorrelation",
        top_annotation = column_ha
       )


#dev.off()


```




<table>
<thead><tr><th></th><th scope=col>G440_L</th><th scope=col>G523_L</th><th scope=col>G729_L</th><th scope=col>G532_L</th><th scope=col>G620_L</th><th scope=col>BT67_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G809_L</th><th scope=col>G549_L</th><th scope=col>G583_L</th></tr></thead>
<tbody>
	<tr><th scope=row>G440_L</th><td>1.0000000</td><td>0.3482114</td><td>0.5549921</td><td>0.6861461</td><td>0.3206071</td><td>0.3608887</td><td>0.5662562</td><td>0.4300672</td><td>0.3939418</td><td>0.4643418</td><td>0.4813248</td></tr>
	<tr><th scope=row>G523_L</th><td>0.3482114</td><td>1.0000000</td><td>0.3348128</td><td>0.3556900</td><td>0.3146358</td><td>0.3813586</td><td>0.3559680</td><td>0.2221849</td><td>0.3531725</td><td>0.3040188</td><td>0.3133939</td></tr>
	<tr><th scope=row>G729_L</th><td>0.5549921</td><td>0.3348128</td><td>1.0000000</td><td>0.6022449</td><td>0.2726996</td><td>0.3476563</td><td>0.5889082</td><td>0.3469828</td><td>0.3705865</td><td>0.4935261</td><td>0.4779381</td></tr>
	<tr><th scope=row>G532_L</th><td>0.6861461</td><td>0.3556900</td><td>0.6022449</td><td>1.0000000</td><td>0.2886741</td><td>0.3736741</td><td>0.6097997</td><td>0.4423470</td><td>0.4273913</td><td>0.5087986</td><td>0.5046214</td></tr>
	<tr><th scope=row>G620_L</th><td>0.3206071</td><td>0.3146358</td><td>0.2726996</td><td>0.2886741</td><td>1.0000000</td><td>0.4973956</td><td>0.3468562</td><td>0.2571895</td><td>0.3575850</td><td>0.3140712</td><td>0.2549570</td></tr>
	<tr><th scope=row>BT67_L</th><td>0.3608887</td><td>0.3813586</td><td>0.3476563</td><td>0.3736741</td><td>0.4973956</td><td>1.0000000</td><td>0.4123455</td><td>0.2707717</td><td>0.3897286</td><td>0.3484513</td><td>0.2812992</td></tr>
</tbody>
</table>




![png](output_11_2.png)



```R


Heatmap(dat,
        col = viridis(60),
        column_km = 2,
        border = T,
        name = "qBF \n(z-score)",
        top_annotation = column_ha,
        row_title = "Varaible Genes (n=1484)",
        column_title = "",
        show_row_names = FALSE
       )
```




![png](output_12_1.png)



```R
pdf("~/Desktop/Cluster_CRISPR_Correlation.pdf", width = 6, height =4)

corr.map

dev.off()
```




<strong>pdf:</strong> 2


-----
### 3.0 Proportion Essential Genes
-----


```R
## load CRISPR data

crispr <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/CRISPR_Avg_Diff_Matrix.rds")
crispr$G440_L <- NULL
crispr$G532_L <- NULL

head(crispr)


## define list of essential IR and Dev genes
## zscore > 2 = IR; 460 genes
## zscore < 2 = Dev; 458 genes

IR_genes <- rownames(crispr[crispr$Diff_Z > 2, ])
Dev_genes <- rownames(crispr[crispr$Diff_Z < -2, ])

###load in bulk  RNA
load("./bulkRNA_data/GSC_gsva.rdata")
GSC.gsva <- GSC.gsva[c("RNA.GSC.c1", "RNA.GSC.c2"), ]
rownames(GSC.gsva) <- c("Developmental_GSVA", "InjuryResponse_GSVA")
GSC.gsva <- data.frame(t(GSC.gsva[, colnames(crispr)[1:9]]))
GSC.gsva$Diff <- c(GSC.gsva$InjuryResponse_GSVA - GSC.gsva$Developmental_GSVA)
GSC.gsva$Diff_Z <- scale(GSC.gsva$Diff)
(GSC.gsva)

samples <- rownames(GSC.gsva)
samples

IR <- c()
Dev <- c()
num.essential <- c()

for (i in 1:length(samples)){

sample <- samples[i]
print(sample)
    
num.essential[i] <- table(crispr[,sample ] > 10)[2] ###qBF cutoff of 5
IR[i] <- table(crispr[IR_genes,sample ] > 10)[2]
Dev[i] <- table(crispr[Dev_genes,sample ] > 10)[2]
    
}

GSC.gsva$Essential.Genes <- num.essential
GSC.gsva$IR_Essential <- IR
GSC.gsva$Dev_Essential <- Dev

##sort by order of Diff_Z

GSC.gsva <- GSC.gsva[order(GSC.gsva$Diff_Z, decreasing = T), ]
head(GSC.gsva)
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>Developmental_Average</th><th scope=col>InjuryResponse_Average</th><th scope=col>Diff</th><th scope=col>Diff_Z</th><th scope=col>FoldChange</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>-16.8872393</td><td>30.37648   </td><td>47.26371   </td><td>7.051752   </td><td> -1.798783 </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td> -8.8287541</td><td>32.77162   </td><td>41.60037   </td><td>6.206774   </td><td> -3.711919 </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td> -0.7493837</td><td>39.18887   </td><td>39.93825   </td><td>5.958783   </td><td>-52.294793 </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>  4.2116815</td><td>43.98490   </td><td>39.77322   </td><td>5.934161   </td><td> 10.443548 </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td> -6.9819679</td><td>37.50056   </td><td>44.48253   </td><td>6.636796   </td><td> -5.371059 </td></tr>
	<tr><th scope=row>AARS</th><td>  0.9270344</td><td>  4.774409 </td><td>-10.0122166</td><td> 34.759483 </td><td>52.83310   </td><td>47.83145   </td><td>37.41367   </td><td>39.54539   </td><td>28.55795   </td><td>  7.6121776</td><td>41.23631   </td><td>33.62414   </td><td>5.016710   </td><td>  5.417151 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Developmental_GSVA</th><th scope=col>InjuryResponse_GSVA</th><th scope=col>Diff</th><th scope=col>Diff_Z</th></tr></thead>
<tbody>
	<tr><th scope=row>BT67_L</th><td> 0.3941965</td><td>-0.3568360</td><td>-0.7510324</td><td>-1.0586786</td></tr>
	<tr><th scope=row>G620_L</th><td> 0.3953486</td><td>-0.3129144</td><td>-0.7082630</td><td>-1.0112305</td></tr>
	<tr><th scope=row>G523_L</th><td> 0.4043064</td><td>-0.3757202</td><td>-0.7800266</td><td>-1.0908447</td></tr>
	<tr><th scope=row>G809_L</th><td> 0.2936245</td><td>-0.3407813</td><td>-0.6344058</td><td>-0.9292936</td></tr>
	<tr><th scope=row>G583_L</th><td>-0.3942071</td><td> 0.4249421</td><td> 0.8191492</td><td> 0.6832733</td></tr>
	<tr><th scope=row>G549_L</th><td>-0.4958445</td><td> 0.5514643</td><td> 1.0473088</td><td> 0.9363925</td></tr>
	<tr><th scope=row>G729_L</th><td>-0.6342635</td><td> 0.6329660</td><td> 1.2672295</td><td> 1.1803714</td></tr>
	<tr><th scope=row>G361_L</th><td>-0.3245410</td><td> 0.1517567</td><td> 0.4762977</td><td> 0.3029155</td></tr>
	<tr><th scope=row>G691_L</th><td>-0.5418078</td><td> 0.5512035</td><td> 1.0930113</td><td> 0.9870946</td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'BT67_L'</li>
	<li>'G620_L'</li>
	<li>'G523_L'</li>
	<li>'G809_L'</li>
	<li>'G583_L'</li>
	<li>'G549_L'</li>
	<li>'G729_L'</li>
	<li>'G361_L'</li>
	<li>'G691_L'</li>
</ol>



    [1] "BT67_L"
    [1] "G620_L"
    [1] "G523_L"
    [1] "G809_L"
    [1] "G583_L"
    [1] "G549_L"
    [1] "G729_L"
    [1] "G361_L"
    [1] "G691_L"



<table>
<thead><tr><th></th><th scope=col>Developmental_GSVA</th><th scope=col>InjuryResponse_GSVA</th><th scope=col>Diff</th><th scope=col>Diff_Z</th><th scope=col>Essential.Genes</th><th scope=col>IR_Essential</th><th scope=col>Dev_Essential</th></tr></thead>
<tbody>
	<tr><th scope=row>G729_L</th><td>-0.6342635</td><td> 0.6329660</td><td> 1.2672295</td><td> 1.1803714</td><td>1370      </td><td>299       </td><td> 60       </td></tr>
	<tr><th scope=row>G691_L</th><td>-0.5418078</td><td> 0.5512035</td><td> 1.0930113</td><td> 0.9870946</td><td>1370      </td><td>271       </td><td> 59       </td></tr>
	<tr><th scope=row>G549_L</th><td>-0.4958445</td><td> 0.5514643</td><td> 1.0473088</td><td> 0.9363925</td><td>1370      </td><td>276       </td><td> 73       </td></tr>
	<tr><th scope=row>G583_L</th><td>-0.3942071</td><td> 0.4249421</td><td> 0.8191492</td><td> 0.6832733</td><td>1371      </td><td>296       </td><td> 69       </td></tr>
	<tr><th scope=row>G361_L</th><td>-0.3245410</td><td> 0.1517567</td><td> 0.4762977</td><td> 0.3029155</td><td>1370      </td><td>280       </td><td> 85       </td></tr>
	<tr><th scope=row>G809_L</th><td> 0.2936245</td><td>-0.3407813</td><td>-0.6344058</td><td>-0.9292936</td><td>1370      </td><td>148       </td><td>214       </td></tr>
</tbody>
</table>




```R
pdf("~/Desktop/PropEssential_test.pdf", width =4, height =4.5)

counts <- t(GSC.gsva[,6:7])
counts

plot(counts[1,]/460, 
     type = 'o', 
     ylim = c(0,1),
     ylab = "Proportion Essential Genes",
     xlab = '',
     xaxt = 'n',
     main = "Injury Response <---- Ranked Samples ----> Developmental",,
     pch = c(20,20),
     cex =1.5,
     col = "#4d4d4d",
     cex.axis = 0.8
    )
lines(counts[2,]/458, 
      type = "o", 
      col = "#b2182b",
      pch = 20,
      cex = 1.5
     )
axis(1, at=c(1:9), labels=colnames(counts), las = 2, cex.axis = 0.8)

legend("topleft", 
       c("Injury Response (460g)", "Developmental (458g)"),
       cex=0.8,
       lwd = 1,
       pch = 20,
       col = c("#4d4d4d", "#b2182b")
      )

dev.off()
```


<table>
<thead><tr><th></th><th scope=col>G729_L</th><th scope=col>G691_L</th><th scope=col>G549_L</th><th scope=col>G583_L</th><th scope=col>G361_L</th><th scope=col>G809_L</th><th scope=col>G620_L</th><th scope=col>BT67_L</th><th scope=col>G523_L</th></tr></thead>
<tbody>
	<tr><th scope=row>IR_Essential</th><td>299</td><td>271</td><td>276</td><td>296</td><td>280</td><td>148</td><td> 85</td><td>123</td><td> 86</td></tr>
	<tr><th scope=row>Dev_Essential</th><td> 60</td><td> 59</td><td> 73</td><td> 69</td><td> 85</td><td>214</td><td>239</td><td>225</td><td>250</td></tr>
</tbody>
</table>




<strong>pdf:</strong> 2



```R
pdf("~/Desktop/Lineplot_annotation.pdf", width = 5, height = 1)

Heatmap(t(GSC.gsva[ ,c(1:2)]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = T,
        rect_gp = gpar(col = "black", lwd = 1)
        )

dev.off()
```




<strong>pdf:</strong> 2


---
## 4.0 Plot integrative heat map 
---

1. CRISPR qBF
2. Bulk RNA 
3. scRNA violin plot


```R
#### READ IN BULK GSVA DATA

sample.order <- c("G729_L",
                  "G691_L",
                  "G549_L",
                  "G583_L",
                  "G361_L",
                  "G809_L",
                  "G620_L",
                  "BT67_L",
                  "G523_L"
                 )


load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/bulkRNA_data/GSC_gsva.rdata")
GSVA <- data.frame(GSC.gsva[c("RNA.GSC.c1", "RNA.GSC.c2"), colnames(GSC.gsva) %in% sample.order])
GSVA$G440_L <- 0
GSVA$G532_L <- 0


#order the same as var.genes
GSVA <- GSVA[ ,sample.order]
head(GSVA)


col_fun = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"))

column_ha <- HeatmapAnnotation(Developmental = as.numeric((GSVA[1, ])),
                               Injury_Response = as.numeric((GSVA[2, ])),
                               col = list(Developmental = col_fun,
                                          Injury_Response = col_fun
                                         ),
                                annotation_legend_param = list(Developmental = list(title = "GSVA"),
                                                               Injury_Response = list(title = "GSVA")
                                                              ),
                               gp = gpar(col = "black")
                              
                              )
```


<table>
<thead><tr><th></th><th scope=col>G729_L</th><th scope=col>G691_L</th><th scope=col>G549_L</th><th scope=col>G583_L</th><th scope=col>G361_L</th><th scope=col>G809_L</th><th scope=col>G620_L</th><th scope=col>BT67_L</th><th scope=col>G523_L</th></tr></thead>
<tbody>
	<tr><th scope=row>RNA.GSC.c1</th><td>-0.6342635</td><td>-0.5418078</td><td>-0.4958445</td><td>-0.3942071</td><td>-0.3245410</td><td> 0.2936245</td><td> 0.3953486</td><td> 0.3941965</td><td> 0.4043064</td></tr>
	<tr><th scope=row>RNA.GSC.c2</th><td> 0.6329660</td><td> 0.5512035</td><td> 0.5514643</td><td> 0.4249421</td><td> 0.1517567</td><td>-0.3407813</td><td>-0.3129144</td><td>-0.3568360</td><td>-0.3757202</td></tr>
</tbody>
</table>




```R
min(crispr.dat)
```


-29.61845628



```R
###LOAD IN CRISPR data

dat <- readRDS("./CRISPR_Avg_Diff_Matrix.rds")
dat <- dat[rev(order(dat$Diff_Z)), ]
dat$GeneRank <- 1:nrow(dat)
dat$Gene <- rownames(dat)
labels <- dat[c(1:10,(nrow(dat)-10):nrow(dat)), ]$Gene #top and bottom 10 genes
dat$Label <- ifelse(dat$Gene %in% labels, dat$Gene, "")
#head(dat)

labels <- dat[c(1:10,(nrow(dat)-9):nrow(dat)), ]$Gene #top and bottom 10 genes

crispr.dat <- dat[labels, grep("_L", colnames(dat))]
crispr.dat$G532_L <- NULL
crispr.dat$G440_L <- NULL
crispr.dat <- as.matrix(crispr.dat[ ,sample.order])
head(crispr.dat)


##### CRISPR heatmap
#### Plot top 10 essential genes

library(circlize)
col_fun = colorRamp2(seq(from = -60, to = 60, len = 11), 
                     rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
                    )

crispr.heat <- Heatmap(crispr.dat,
                        cluster_columns = FALSE,
                        row_km = 2,
                       col = col_fun,
                        row_names_side = "right",
                        #row_split = c("Dev", "IR"),
                        border = TRUE,
                        row_gap = unit(0, "mm"),
                       name = "Fitness Score\n(qBF)",
                       heatmap_legend_param = list(direction = "horizontal"),
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       top_annotation = column_ha
                         )

crispr.heat
gene.order <- rownames(crispr.dat)

```


<table>
<thead><tr><th></th><th scope=col>G729_L</th><th scope=col>G691_L</th><th scope=col>G549_L</th><th scope=col>G583_L</th><th scope=col>G361_L</th><th scope=col>G809_L</th><th scope=col>G620_L</th><th scope=col>BT67_L</th><th scope=col>G523_L</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>29.13527   </td><td>31.33191   </td><td>35.34420   </td><td>27.04493   </td><td>29.02607   </td><td> -4.322983 </td><td>-15.353046 </td><td>-28.2717311</td><td>-19.6011966</td></tr>
	<tr><th scope=row>PPP1R12A</th><td>34.50681   </td><td>42.51892   </td><td>55.71300   </td><td>35.65962   </td><td>19.10447   </td><td> -6.449141 </td><td> -4.982018 </td><td>  0.1890148</td><td>-16.6857273</td></tr>
	<tr><th scope=row>SCAP</th><td>42.90056   </td><td>17.01073   </td><td>46.03781   </td><td>38.35182   </td><td>19.55718   </td><td>-24.391670 </td><td>-14.911093 </td><td> 32.5498409</td><td>-28.5620940</td></tr>
	<tr><th scope=row>EXOSC5</th><td>41.65691   </td><td>31.50877   </td><td>30.57221   </td><td>45.94805   </td><td>46.25838   </td><td> 22.298264 </td><td>-10.239477 </td><td> -5.7701043</td><td> -9.2862183</td></tr>
	<tr><th scope=row>SUPT16H</th><td>46.14875   </td><td>48.30201   </td><td>42.16704   </td><td>39.91896   </td><td>43.38775   </td><td> 21.954226 </td><td> -1.281364 </td><td> -2.8559548</td><td> -0.9701818</td></tr>
	<tr><th scope=row>GTF2A2</th><td>27.19230   </td><td>22.41802   </td><td>39.81692   </td><td>39.48553   </td><td>30.89674   </td><td>  4.311062 </td><td> -7.656784 </td><td>  5.5289992</td><td> -9.0392470</td></tr>
</tbody>
</table>






![png](output_21_2.png)



```R
gene.order


```


<ol class=list-inline>
	<li>'ITGB1'</li>
	<li>'PPP1R12A'</li>
	<li>'SCAP'</li>
	<li>'EXOSC5'</li>
	<li>'SUPT16H'</li>
	<li>'GTF2A2'</li>
	<li>'AARS'</li>
	<li>'FOSL1'</li>
	<li>'UBL5'</li>
	<li>'ILK'</li>
	<li>'SOX6'</li>
	<li>'OGDH'</li>
	<li>'EED'</li>
	<li>'ASCL1'</li>
	<li>'CCND2'</li>
	<li>'AHR'</li>
	<li>'POLE3'</li>
	<li>'SOX2'</li>
	<li>'IRS2'</li>
	<li>'OLIG2'</li>
</ol>


# LOAD AND FORMAT SINGLE CELL RNA DATA
#### THIS WILL BE AN ANNOTATION ON 
### Need to load in single data and make two matrices

library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" )

##get list of Dev cells
## get list of IR cells

cnv < readRDS("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/G800_removed/Cell_armlevelCNVs_Dec222019.rds")
head(cnv)

dev.cells <- rownames(cnv[cnv$ID == "Dev", ])
length(dev.cells)

IR.cells <- rownames(cnv[cnv$ID == "IR", ])
length(IR.cells)

### now load the expression matrx and subset by genes in crispr

genes <- c('ITGB1',
           'PPP1R12A',
           'SCAP',
           'EXOSC5',
           'SUPT16H',
           'GTF2A2',
           'AARS',
           'FOSL1',
           'UBL5',
           'ILK',
           'SOX6',
           'OGDH',
           'EED',
           'ASCL1',
           'CCND2',
           'AHR',
           'POLE3',
           'SOX2',
           'IRS2',
           'OLIG2'
           )

BTSC <- readRDS("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSCs_RemoveG800L/BTSC_G800L_removed_AUCell_SeuratObj.rds")

DevExp <- as.matrix(BTSC@data[genes, dev.cells])
dim(DevExp)

IRExp <- as.matrix(BTSC@data[genes, IR.cells])
dim(IRExp)

saveRDS(DevExp, file = "DevCells_CRISPR_Expression.rds")
saveRDS(IRExp, file = "IRCells_CRISPR_Expression.rds")

```R
Dev <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/DevCells_CRISPR_Expression.rds")
IR <- readRDS("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/IRCells_CRISPR_Expression.rds")

```


```R
Dev <- t(Dev)
head(Dev)

IR <- t(IR)
head(IR)
```


<table>
<thead><tr><th></th><th scope=col>ITGB1</th><th scope=col>PPP1R12A</th><th scope=col>SCAP</th><th scope=col>EXOSC5</th><th scope=col>SUPT16H</th><th scope=col>GTF2A2</th><th scope=col>AARS</th><th scope=col>FOSL1</th><th scope=col>UBL5</th><th scope=col>ILK</th><th scope=col>SOX6</th><th scope=col>OGDH</th><th scope=col>EED</th><th scope=col>ASCL1</th><th scope=col>CCND2</th><th scope=col>AHR</th><th scope=col>POLE3</th><th scope=col>SOX2</th><th scope=col>IRS2</th><th scope=col>OLIG2</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>0.0000000</td><td>0.0000000</td><td>0        </td><td>0        </td><td>0.0000000</td><td>0.000000 </td><td>0        </td><td>0        </td><td>2.519998 </td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.000000 </td><td>0.0000000</td><td>0        </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>0.0000000</td><td>0.0000000</td><td>0        </td><td>0        </td><td>0.0000000</td><td>2.230622 </td><td>0        </td><td>0        </td><td>1.639545 </td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>1.6395449</td><td>0.0000000</td><td>0.0000000</td><td>0.000000 </td><td>0.0000000</td><td>0        </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>0.0000000</td><td>0.6902597</td><td>0        </td><td>0        </td><td>0.6902597</td><td>1.786942 </td><td>0        </td><td>0        </td><td>1.604814 </td><td>0.6902597</td><td>1.3819601</td><td>0.0000000</td><td>0.0000000</td><td>1.0947605</td><td>1.6048140</td><td>0.0000000</td><td>0.0000000</td><td>1.940955 </td><td>1.3819601</td><td>0        </td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>1.7281960</td><td>1.0480888</td><td>0        </td><td>0        </td><td>1.0480888</td><td>1.329272 </td><td>0        </td><td>0        </td><td>1.880470 </td><td>0.6554959</td><td>0.6554959</td><td>0.6554959</td><td>0.0000000</td><td>1.3292721</td><td>1.0480888</td><td>1.0480888</td><td>0.6554959</td><td>1.048089 </td><td>1.0480888</td><td>0        </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>0.8627926</td><td>1.3189622</td><td>0        </td><td>0        </td><td>0.5216971</td><td>1.868578 </td><td>0        </td><td>0        </td><td>1.756856 </td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>0.5216971</td><td>0.5216971</td><td>0.8627926</td><td>0.5216971</td><td>0.5216971</td><td>1.756856 </td><td>0.8627926</td><td>0        </td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCATTGGTAC</th><td>1.7276630</td><td>2.3277629</td><td>0        </td><td>0        </td><td>0.0000000</td><td>0.000000 </td><td>0        </td><td>0        </td><td>0.000000 </td><td>1.7276630</td><td>1.7276630</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>2.3277629</td><td>0.0000000</td><td>0.0000000</td><td>3.638110 </td><td>0.0000000</td><td>0        </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>ITGB1</th><th scope=col>PPP1R12A</th><th scope=col>SCAP</th><th scope=col>EXOSC5</th><th scope=col>SUPT16H</th><th scope=col>GTF2A2</th><th scope=col>AARS</th><th scope=col>FOSL1</th><th scope=col>UBL5</th><th scope=col>ILK</th><th scope=col>SOX6</th><th scope=col>OGDH</th><th scope=col>EED</th><th scope=col>ASCL1</th><th scope=col>CCND2</th><th scope=col>AHR</th><th scope=col>POLE3</th><th scope=col>SOX2</th><th scope=col>IRS2</th><th scope=col>OLIG2</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_ACCAGTAGTAGGGACT</th><td>0.6056060</td><td>0.9801007</td><td>0.0000000</td><td>0.0000000</td><td>0.6056060</td><td>1.2519303</td><td>0.0000000</td><td>0.6056060</td><td>1.465440 </td><td>0.0000000</td><td>0        </td><td>0.0000000</td><td>0        </td><td>0.0000000</td><td>1.2519303</td><td>0        </td><td>0.9801007</td><td>0.9801007</td><td>0.9801007</td><td>0.6056060</td></tr>
	<tr><th scope=row>BT127_L_AGGTCCGTCTCGGACG</th><td>1.0909711</td><td>0.6874218</td><td>0.4016518</td><td>0.4016518</td><td>0.6874218</td><td>1.0909711</td><td>0.0000000</td><td>0.4016518</td><td>1.090971 </td><td>0.6874218</td><td>0        </td><td>0.6874218</td><td>0        </td><td>0.0000000</td><td>1.4951570</td><td>0        </td><td>0.4016518</td><td>1.9360750</td><td>1.2445737</td><td>0.0000000</td></tr>
	<tr><th scope=row>BT127_L_ATAAGAGGTCACCTAA</th><td>1.1547791</td><td>0.7355653</td><td>0.0000000</td><td>0.0000000</td><td>0.0000000</td><td>1.5693062</td><td>1.1547791</td><td>0.4339428</td><td>1.942567 </td><td>0.7355653</td><td>0        </td><td>0.0000000</td><td>0        </td><td>1.1547791</td><td>1.1547791</td><td>0        </td><td>0.9669807</td><td>0.7355653</td><td>1.6764692</td><td>0.4339428</td></tr>
	<tr><th scope=row>BT127_L_CTACGTCAGTCAATAG</th><td>0.0000000</td><td>1.1517367</td><td>0.0000000</td><td>0.7332527</td><td>0.4323801</td><td>1.3095762</td><td>0.0000000</td><td>0.0000000</td><td>1.857735 </td><td>0.4323801</td><td>0        </td><td>0.4323801</td><td>0        </td><td>0.4323801</td><td>1.1517367</td><td>0        </td><td>0.9642279</td><td>1.5657879</td><td>0.0000000</td><td>0.0000000</td></tr>
	<tr><th scope=row>BT127_L_CTACGTCGTTCCACTC</th><td>1.4903151</td><td>1.2750472</td><td>0.8283052</td><td>0.6203740</td><td>0.0000000</td><td>1.0003550</td><td>0.8283052</td><td>0.0000000</td><td>2.215547 </td><td>0.6203740</td><td>0        </td><td>0.3575428</td><td>0        </td><td>0.0000000</td><td>0.8283052</td><td>0        </td><td>1.1471036</td><td>0.0000000</td><td>0.3575428</td><td>0.0000000</td></tr>
	<tr><th scope=row>BT127_L_CTGTGCTGTGAGCGAT</th><td>0.5548787</td><td>0.5548787</td><td>0.0000000</td><td>0.0000000</td><td>0.5548787</td><td>0.9096525</td><td>0.0000000</td><td>1.1709916</td><td>1.823276 </td><td>0.0000000</td><td>0        </td><td>0.5548787</td><td>0        </td><td>0.0000000</td><td>0.5548787</td><td>0        </td><td>0.5548787</td><td>1.1709916</td><td>1.1709916</td><td>0.0000000</td></tr>
</tbody>
</table>




```R
ha <- rowAnnotation(IR = anno_boxplot(t(IR), type = "violin", width = unit(1.5, "cm"), outline = F),
                   Dev = anno_boxplot(t(Dev), type = "violin", width = unit(1.5, "cm"),outline = F)
                  )
#draw(ha)


```


```R
library(circlize)
col_fun = colorRamp2(seq(from = -2, to = 2, len = 11), 
                     rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
                    )
```


```R
min(fpkm)
max(fpkm)
```


-2.2770372037826



1.97101609127626



```R
#### LOAD AND FORMAT BULK RNA DATA

library("SummarizedExperiment")
load("~/Desktop/SU2C_scRNA_Manuscript/RNA_adult_GBM_line_fpkm_combat_nonSU2C.rdata")
fpkm <- as.matrix(RNA_adult_GBM_line_fpkm_combat_nonSU2C@assays$data$fpkm)

include <- c("A61502_Dirks_G729_Line",
             "A67980_G691_line",
             "A61517_Dirks_G549_Line",
             "A61523_Dirks_G583_Line",
             "G361_L",
             "A67956_G809_line",
             "A61516_Dirks_G620_Line",
             "A61532_Weiss_BT67_Line",
             "A61519_Dirks_G523_Line"
            )

fpkm <- fpkm[labels, include]
fpkm <- log2(fpkm+1)
fpkm <- scale(t(fpkm))
fpkm <- t(fpkm)
colnames(fpkm) <- sample.order
fpkm <- fpkm[gene.order, ]
#head(fpkm)


### bulkRNA heatmap

bulkRNA.heat <- Heatmap(fpkm,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        col = col_fun,
                        #row_km = 2,
                        row_names_side = "right",
                        column_names_side = "bottom",
                        #row_split = c("Dev", "IR"),
                        border = TRUE,
                        row_gap = unit(0, "mm"),
                        name = "Bulk RNA expression\n (log2(FPKM+1), z-score)",
                        heatmap_legend_param = list(direction = "horizontal"),
                        right_annotation = ha,
                         row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        top_annotation = column_ha
                         )
bulkRNA.heat
```




![png](output_29_1.png)



```R
## COMBINE HEATMAPS

pdf("~/Desktop/test_heatmap.pdf", height = 5, width =7)

ht_list = crispr.heat + bulkRNA.heat
draw(ht_list, ht_gap = unit(0.1, "cm"), heatmap_legend_side = "top")



dev.off()
```


<strong>pdf:</strong> 2



```R

```


```R

```
