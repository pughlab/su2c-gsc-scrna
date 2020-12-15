
---
## Define thresholds for CNV gain/loss in scRNA-seq
---

L.Richards

We want to define a threshold for a gene or chr region being gained/lost in scRNA-seq data. To do this, we will compare calls from GISTIC (WGS) and InferCNV. 

> Define cohort of scRNA-seq lines that also have WGS calls  
> Identify union of genes that is present in GISTIC and scRNA-seq  
> Plot InferCNV score against GISTIC scores  
> Fit a linear relationship  
> Use the model to determine which InferCNV values are equivalent to a single copy loss or gain  


/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs
##### Average CNV profiles across samples

options(stringsAsFactors = FALSE)
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/")
                )
library(reshape2)

#load InferCNV calls at the average cluster level
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/GlobalBTSC_CNVs.RData")
str(avg.cnv.df)


samples <- as.character(unique(intra.meta$orig.ident))
samples

sample.avg <- list()

for (i in 1:length(samples)){
    
    print(samples[i])
    dat <- BTSC.CNVs[ ,grep(samples[i], colnames(BTSC.CNVs))]
    dat <- data.frame(rowMeans(dat))
    colnames(dat) <- samples[i]
    sample.avg[[samples[i]]] <- dat

    
}

df <- do.call("cbind", sample.avg)

saveRDS(df, file = "BTSC_SampleAvgCNVs.rds")

```R
#### combine WGS and sample average

avg.sample.cnv <- readRDS("BTSC_SampleAvgCNVs.rds")

#load GISTIC results from WGS
gistic.file <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_thresholded.by_genes.txt"
gistic <- read.table(gistic.file, header = T, sep = "\t")
gistic$Locus.ID <- NULL
gistic$Cytoband <- NULL

#identify common samples between scRNA and WGS
#20 samples overlap between the WGS and single cell RNA-sequencing data
scrna.samples <- unique(colnames(avg.sample.cnv))
wgs.samples <- unique(colnames(gistic))[grep("_L", unique(colnames(gistic)))]
samples <- scrna.samples[scrna.samples %in% wgs.samples]
length(samples) 

#identfy common genes bewteen scRNA and WGS
#6351 genes in scRNa df
#6241 genes overlap between the WGS and single cell RNA-sequencing data
scrna.genes <- rownames(avg.sample.cnv)
wgs.genes <- gistic$Gene.Symbol
genes <- scrna.genes[scrna.genes %in% wgs.genes]
length(genes)

#make a dataframe combining scRNA and WGS data together

combined <- avg.sample.cnv[genes, grepl(paste(samples, collapse="|"),  colnames(avg.sample.cnv))]
combined$Gene <- rownames(combined)
combined <- melt(data = combined, id = "Gene")
colnames(combined) <- c("Gene", "scRNA_sample", "scRNA_Value")
combined$ID <- paste0(combined$Gene, "_", combined$scRNA_sample)
head(combined)

##merge the WGS data
rownames(gistic) <- gistic$Gene.Symbol
gistic$Gene.Symbol <- NULL
gistic <- gistic[genes, samples]
gistic$Gene <- rownames(gistic)
gistic <- melt(data = gistic, id = "Gene")
gistic$ID <- paste0(gistic$Gene, "_", gistic$variable)
head(gistic)

wgs <- c()
for(i in 1:nrow(combined)){    
    df <- gistic[gistic$ID == combined$ID[i], ]
    wgs[i] <- df$value 
}

combined$WGS <- wgs
head(combined)


########### Add log2 values

seg <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_data_by_genes.txt"
seg <- read.table(seg, header = T, sep = "\t")
seg$Gene.ID <- NULL
seg$Cytoband <- NULL

rownames(seg) <- seg$Gene.Symbol
seg$Gene.Symbol <- NULL
seg <- seg[genes, samples]
seg$Gene <- rownames(seg)
seg <- melt(data = seg, id = "Gene")
seg$ID <- paste0(seg$Gene, "_", seg$variable)
head(seg)

log2 <- c()
for(i in 1:nrow(combined)){    
    df <- seg[seg$ID == combined$ID[i], ]
    log2[i] <- df$value 
}

combined$log2 <- log2
head(combined)

#make the GISTIC scores into binary gain/loss
#positive = gain
#negative = loss

colnames(combined)[5] <- "GISTIC"
combined$Category <- as.numeric(combined$GISTIC)

combined$Category[combined$Category > 0] <- "Gain"
combined$Category[combined$Category == -1] <- "Deletion"
combined$Category[combined$Category == -2] <- "Deletion"
combined$Category[combined$Category == 0] <- "Neutral"

save(combined, file = "Matched_scRNA_WGS_CNV.Rdata")



```

---
## 1.0 Combine single cell and WGS CNV data
---

we excluded clusters that 


```R
options(stringsAsFactors = FALSE)
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/")
                )
library(reshape2)

#load InferCNV calls at the average cluster level
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/GlobalBTSC_CNVs.RData")
str(avg.cnv.df)



#load GISTIC results from WGS
gistic.file <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_thresholded.by_genes.txt"
gistic <- read.table(gistic.file, header = T, sep = "\t")
gistic$Locus.ID <- NULL
gistic$Cytoband <- NULL

#identify common samples between scRNA and WGS
#20 samples overlap between the WGS and single cell RNA-sequencing data
scrna.samples <- unique(gsub('.{3}$', '', colnames(avg.cnv.df)))
wgs.samples <- unique(colnames(gistic))[grep("_L", unique(colnames(gistic)))]
samples <- scrna.samples[scrna.samples %in% wgs.samples]
length(samples) 

#identfy common genes bewteen scRNA and WGS
#6351 genes in scRNa df
#6241 genes overlap between the WGS and single cell RNA-sequencing data
scrna.genes <- rownames(avg.cnv.df)
wgs.genes <- gistic$Gene.Symbol
genes <- scrna.genes[scrna.genes %in% wgs.genes]
length(genes)

#make a dataframe combining scRNA and WGS data together
```


```R
combined <- avg.cnv.df[genes, grepl(paste(samples, collapse="|"),  colnames(avg.cnv.df))]
combined$Gene <- rownames(combined)
combined <- melt(data = combined, id = "Gene")
colnames(combined) <- c("Gene", "scRNA_cluster", "scRNA_Value")
combined$ID <- paste0(combined$Gene, "_", gsub('.{3}$', '', combined$scRNA_cluster))
head(combined)
```


```R
rownames(gistic) <- gistic$Gene.Symbol
gistic$Gene.Symbol <- NULL
gistic <- gistic[genes, samples]
gistic$Gene <- rownames(gistic)
gistic <- melt(data = gistic, id = "Gene")
gistic$ID <- paste0(gistic$Gene, "_", gistic$variable)
head(gistic)

wgs <- c()
for(i in 1:nrow(combined)){    
    df <- gistic[gistic$ID == combined$ID[i], ]
    wgs[i] <- df$value 
}

combined$WGS <- wgs
head(combined)

#save(combined, file = "Matched_scRNA_WGS_CNV.Rdata")
```

#### Add log2 values for genes from GISTIC


```R
seg <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_data_by_genes.txt"
seg <- read.table(seg, header = T, sep = "\t")
seg$Gene.ID <- NULL
seg$Cytoband <- NULL

rownames(seg) <- seg$Gene.Symbol
seg$Gene.Symbol <- NULL
seg <- seg[genes, samples]
seg$Gene <- rownames(seg)
seg <- melt(data = seg, id = "Gene")
seg$ID <- paste0(seg$Gene, "_", seg$variable)
head(seg)

log2 <- c()
for(i in 1:nrow(combined)){    
    df <- seg[seg$ID == combined$ID[i], ]
    log2[i] <- df$value 
}

combined$log2 <- log2
head(combined)

#make the GISTIC scores into binary gain/loss
#positive = gain
#negative = loss

colnames(combined)[5] <- "GISTIC"
combined$Category <- as.numeric(combined$GISTIC)

combined$Category[combined$Category > 0] <- "Gain"
combined$Category[combined$Category == -1] <- "Deletion"
combined$Category[combined$Category == -2] <- "Deletion"
combined$Category[combined$Category == 0] <- "Neutral"

save(combined, file = "Matched_scRNA_WGS_CNV.Rdata")
```

---
## Correlate WGS and single-cell RNA-seq data
---


```R
library(ggpubr)
library(ggplot2)
library(gridExtra)
```

    Warning message:
    “package ‘ggpubr’ was built under R version 3.4.4”Loading required package: ggplot2
    Warning message:
    “package ‘ggplot2’ was built under R version 3.4.4”Loading required package: magrittr



```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/Matched_scRNA_WGS_CNV.Rdata")
head(combined)
```


<table>
<thead><tr><th scope=col>Gene</th><th scope=col>scRNA_cluster</th><th scope=col>scRNA_Value</th><th scope=col>ID</th><th scope=col>GISTIC</th><th scope=col>log2</th><th scope=col>Category</th></tr></thead>
<tbody>
	<tr><td>SLC25A33        </td><td>BT147_L_C1      </td><td>0.03504792      </td><td>SLC25A33_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>TMEM201         </td><td>BT147_L_C1      </td><td>0.03990202      </td><td>TMEM201_BT147_L </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CLSTN1          </td><td>BT147_L_C1      </td><td>0.03787451      </td><td>CLSTN1_BT147_L  </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CTNNBIP1        </td><td>BT147_L_C1      </td><td>0.03495229      </td><td>CTNNBIP1_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>LZIC            </td><td>BT147_L_C1      </td><td>0.04058198      </td><td>LZIC_BT147_L    </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>NMNAT1          </td><td>BT147_L_C1      </td><td>0.04250626      </td><td>NMNAT1_BT147_L  </td><td>-1              </td><td>-0.953          </td><td>Deletion        </td></tr>
</tbody>
</table>




```R
    ggscatter(combined, 
          x = "scRNA_Value", 
          y = "GISTIC",
          color = "black", 
          shape = 21, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          #title = unique(dat$scRNA_cluster)
              )
```




![png](output_12_1.png)



```R
pdf("~/Desktop/WGS_scRNA_corr.pdf", height = 4, width = 4)

smoothScatter(x = combined$scRNA_Value, 
              y = combined$log2,
              xlab = "InferCNV Score (scRNA-seq)",
              ylab = "Log2 Ratio (WGS)",
              main = "r=0.63, p < 2.2e-16",
              cex.axis = 0.8
             )
abline(lm(combined$log2 ~ combined$scRNA_Value), 
       lwd = 2,
       col = "red"
      )


smoothScatter(x = combined$scRNA_Value, 
              y = combined$GISTIC,
              xlab = "InferCNV Score (scRNA-seq)",
              ylab = "GISTIC (WGS)",
              main = "r=0.57, p < 2.2e-16",
               cex.axis = 0.8
             )
abline(lm(combined$GISTIC ~ combined$scRNA_Value), 
       lwd = 2,
       col = "red"
      )

dev.off()
```


<strong>pdf:</strong> 2


---
## 2.0 Correlate cluster CNV profiles and WGS profiles
---


```R
library(ggpubr)
library(ggplot2)
library(gridExtra)
```


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/Matched_scRNA_WGS_CNV.Rdata")
head(combined)
```


<table>
<thead><tr><th scope=col>Gene</th><th scope=col>scRNA_cluster</th><th scope=col>scRNA_Value</th><th scope=col>ID</th><th scope=col>GISTIC</th><th scope=col>log2</th><th scope=col>Category</th></tr></thead>
<tbody>
	<tr><td>SLC25A33        </td><td>BT147_L_C1      </td><td>0.03504792      </td><td>SLC25A33_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>TMEM201         </td><td>BT147_L_C1      </td><td>0.03990202      </td><td>TMEM201_BT147_L </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CLSTN1          </td><td>BT147_L_C1      </td><td>0.03787451      </td><td>CLSTN1_BT147_L  </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CTNNBIP1        </td><td>BT147_L_C1      </td><td>0.03495229      </td><td>CTNNBIP1_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>LZIC            </td><td>BT147_L_C1      </td><td>0.04058198      </td><td>LZIC_BT147_L    </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>NMNAT1          </td><td>BT147_L_C1      </td><td>0.04250626      </td><td>NMNAT1_BT147_L  </td><td>-1              </td><td>-0.953          </td><td>Deletion        </td></tr>
</tbody>
</table>




```R
clusters <- as.character(unique(combined$scRNA_cluster))
clusters

df.list <- list()

for (i in 1:length(clusters)){
    
    sample <- clusters[i]
    #print(sample)
    
    dat <- combined[combined$scRNA_cluster == clusters[i], ]
    df.list[[sample]] <- dat
    
}

#str(df.list)
```


<ol class=list-inline>
	<li>'BT147_L_C1'</li>
	<li>'BT147_L_C2'</li>
	<li>'BT67_L_C2'</li>
	<li>'BT67_L_C1'</li>
	<li>'BT73_L_C1'</li>
	<li>'BT73_L_C2'</li>
	<li>'BT84_L_C1'</li>
	<li>'BT84_L_C2'</li>
	<li>'BT84_L_C3'</li>
	<li>'BT89_L_C2'</li>
	<li>'BT89_L_C1'</li>
	<li>'BT94_L_C2'</li>
	<li>'BT94_L_C1'</li>
	<li>'G523_L_C1'</li>
	<li>'G523_L_C2'</li>
	<li>'G523_L_C3'</li>
	<li>'G549_L_C1'</li>
	<li>'G549_L_C2'</li>
	<li>'G564_L_C2'</li>
	<li>'G564_L_C1'</li>
	<li>'G566_L_C1'</li>
	<li>'G566_L_C2'</li>
	<li>'G566_L_C3'</li>
	<li>'G583_L_C1'</li>
	<li>'G583_L_C3'</li>
	<li>'G583_L_C2'</li>
	<li>'G583_L_C4'</li>
	<li>'G620_L_C1'</li>
	<li>'G620_L_C5'</li>
	<li>'G620_L_C2'</li>
	<li>'G620_L_C3'</li>
	<li>'G620_L_C6'</li>
	<li>'G620_L_C4'</li>
	<li>'G637_L_C3'</li>
	<li>'G637_L_C1'</li>
	<li>'G637_L_C2'</li>
	<li>'G637_L_C4'</li>
	<li>'G637_L_C5'</li>
	<li>'G729_L_C2'</li>
	<li>'G729_L_C1'</li>
	<li>'G799_L_C1'</li>
	<li>'G799_L_C3'</li>
	<li>'G799_L_C2'</li>
	<li>'G800_L_C1'</li>
	<li>'G800_L_C2'</li>
	<li>'G800_L_C3'</li>
	<li>'G851_L_C1'</li>
	<li>'G851_L_C2'</li>
	<li>'G876_L_C1'</li>
	<li>'G876_L_C2'</li>
	<li>'G876_L_C3'</li>
	<li>'G885_L_C2'</li>
	<li>'G885_L_C3'</li>
	<li>'G885_L_C1'</li>
	<li>'G895_L_C1'</li>
	<li>'G895_L_C2'</li>
	<li>'G895_L_C3'</li>
</ol>




```R
plot_scRNA_log2 <- function(dat){
    
    ggscatter(dat, 
          x = "scRNA_Value", 
          y = "log2",
          color = "black", 
          shape = 21, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          title = unique(dat$scRNA_cluster)
   )

    }


plot_scRNA_GISTIC <- function(dat){
    
    ggscatter(dat, 
          x = "scRNA_Value", 
          y = "GISTIC",
          color = "black", 
          shape = 21, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          title = unique(dat$scRNA_cluster)
   )

    }
```


```R
p <- list()

for(i in 1:length(df.list)){
    
    p[[i]] <- plot_scRNA_log2(df.list[[i]])
                                                             
}

x <- list()

for(i in 1:length(df.list)){
    
    x[[i]] <- plot_scRNA_GISTIC(df.list[[i]])
                                                             
}
```


```R
pdf("~/Desktop/CNV_scRNA_log2_Sept2019.pdf", height = 35, width = 35)

do.call(grid.arrange, p)
do.call(grid.arrange, x)

dev.off()
```


<strong>pdf:</strong> 2


---
## 3.0 Plot box plot of inferCNV values for each GISTIC gene state
---




```R
library(ggplot2)
library(plyr)
```


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/Matched_scRNA_WGS_CNV.Rdata")
head(combined)

combined$Category <- factor(combined$Category, levels=c("Deletion", "Neutral", "Gain"))
```


<table>
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>Gene</th><th scope=col>scRNA_sample</th><th scope=col>scRNA_Value</th><th scope=col>ID</th><th scope=col>GISTIC</th><th scope=col>log2</th><th scope=col>Category</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>SLC25A33</td><td>BT147_L</td><td>0.02513320</td><td>SLC25A33_BT147_L</td><td>-1</td><td>-0.931</td><td>Deletion</td></tr>
	<tr><th scope=row>2</th><td>TMEM201 </td><td>BT147_L</td><td>0.02908419</td><td>TMEM201_BT147_L </td><td>-1</td><td>-0.931</td><td>Deletion</td></tr>
	<tr><th scope=row>3</th><td>CLSTN1  </td><td>BT147_L</td><td>0.02562848</td><td>CLSTN1_BT147_L  </td><td>-1</td><td>-0.931</td><td>Deletion</td></tr>
	<tr><th scope=row>4</th><td>CTNNBIP1</td><td>BT147_L</td><td>0.02228264</td><td>CTNNBIP1_BT147_L</td><td>-1</td><td>-0.931</td><td>Deletion</td></tr>
	<tr><th scope=row>5</th><td>LZIC    </td><td>BT147_L</td><td>0.02864583</td><td>LZIC_BT147_L    </td><td>-1</td><td>-0.931</td><td>Deletion</td></tr>
	<tr><th scope=row>6</th><td>NMNAT1  </td><td>BT147_L</td><td>0.02988958</td><td>NMNAT1_BT147_L  </td><td>-1</td><td>-0.953</td><td>Deletion</td></tr>
</tbody>
</table>




```R
length(table(combined$scRNA_sample))

```


20



```R
a <-    ggscatter(combined, 
          x = "scRNA_Value", 
          y = "log2",
          color = "black", 
          shape = 18, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n")
   )

```


```R
plot(combined$scRNA_Value,
     combined$log2, xlab="X label", ylab="Y label", pch=19, cex=.4)

# Draw the colored contour lines
contour(combined$scRNA_Value, combined$log2, drawlabels=FALSE, add=TRUE, lwd=2)

```


    Error in contour.default(combined$scRNA_Value, combined$log2, drawlabels = FALSE, : increasing 'x' and 'y' values expected
    Traceback:


    1. contour(combined$scRNA_Value, combined$log2, drawlabels = FALSE, 
     .     add = TRUE, lwd = 2)

    2. contour.default(combined$scRNA_Value, combined$log2, drawlabels = FALSE, 
     .     add = TRUE, lwd = 2)

    3. stop("increasing 'x' and 'y' values expected")



![png](output_26_1.png)



```R
smoothScatter( combined$GISTIC, combined$scRNA_Value)
abline(h=1)

ggscatter(combined, 
          x = "GISTIC", 
          y = "scRNA_Value",
          color = "black", 
          shape = 18, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n")
   )

```




![png](output_27_1.png)



![png](output_27_2.png)



```R
combined$GISTIC <- factor(combined$GISTIC)

p_meds <- ddply(combined, .(GISTIC), summarise, med = median(scRNA_Value))
p_meds

p <- ggplot(combined, aes(x=GISTIC, y=scRNA_Value, color = GISTIC)) + geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
scale_color_manual(values=c("darkblue", "darkblue", "black", "darkred", "darkred")) + theme_classic() 



    #geom_text(data = p_meds, aes(x = WGS, y = med, label = round(med, 2)), 
     #         size = 3, vjust = -1.5)


```


<table>
<thead><tr><th scope=col>GISTIC</th><th scope=col>med</th></tr></thead>
<tbody>
	<tr><td>-2          </td><td>-0.108108039</td></tr>
	<tr><td>-1          </td><td>-0.155283495</td></tr>
	<tr><td>0           </td><td>-0.002213527</td></tr>
	<tr><td>1           </td><td> 0.174804595</td></tr>
	<tr><td>2           </td><td> 0.231754441</td></tr>
</tbody>
</table>






![png](output_28_2.png)



```R


q_meds <- ddply(combined, .(Category), summarise, med = median(scRNA_Value))
q_meds

q <- ggplot(combined, aes(x=Category, y=scRNA_Value, color = Category)) + geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
scale_color_manual(values=c("darkblue", "black", "darkred")) + theme_classic() +
 geom_hline(yintercept = q_meds$med[c(1,3)], lty = 2, col = "grey") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        legend.position = "none"
       ) + xlab("") + ylab("")



    #geom_text(data = p_meds, aes(x = WGS, y = med, label = round(med, 2)), 
     #         size = 3, vjust = -1.5)


pdf("~/Desktop/scRNA_WGS_CNV_cutoffs.pdf", height = 6, width = 6)
q
dev.off()

```


<table>
<thead><tr><th scope=col>Category</th><th scope=col>med</th></tr></thead>
<tbody>
	<tr><td>Deletion    </td><td>-0.153105381</td></tr>
	<tr><td>Neutral     </td><td>-0.002213527</td></tr>
	<tr><td>Gain        </td><td> 0.177383851</td></tr>
</tbody>
</table>






<strong>pdf:</strong> 2



```R
##### pdf("~/Desktop/scRNA_WGS_CNV.pdf")
a
p
q
dev.off()
```

---
# Correlate log2 and GISTIC scores with sample averaged scRNA-seq data
---






```R
###load data and average by sample ID
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(plyr)
```


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/Matched_scRNA_WGS_CNV.Rdata")
```


```R
head(combined)
```


<table>
<thead><tr><th scope=col>Gene</th><th scope=col>scRNA_sample</th><th scope=col>scRNA_Value</th><th scope=col>ID</th><th scope=col>GISTIC</th><th scope=col>log2</th><th scope=col>Category</th></tr></thead>
<tbody>
	<tr><td>SLC25A33        </td><td>BT147_L         </td><td>0.02513320      </td><td>SLC25A33_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>TMEM201         </td><td>BT147_L         </td><td>0.02908419      </td><td>TMEM201_BT147_L </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CLSTN1          </td><td>BT147_L         </td><td>0.02562848      </td><td>CLSTN1_BT147_L  </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>CTNNBIP1        </td><td>BT147_L         </td><td>0.02228264      </td><td>CTNNBIP1_BT147_L</td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>LZIC            </td><td>BT147_L         </td><td>0.02864583      </td><td>LZIC_BT147_L    </td><td>-1              </td><td>-0.931          </td><td>Deletion        </td></tr>
	<tr><td>NMNAT1          </td><td>BT147_L         </td><td>0.02988958      </td><td>NMNAT1_BT147_L  </td><td>-1              </td><td>-0.953          </td><td>Deletion        </td></tr>
</tbody>
</table>




```R
a <-    ggscatter(combined, 
          x = "scRNA_Value", 
          y = "GISTIC",
          color = "black", 
          shape = 18, 
          size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n")
   )

a
```




![png](output_35_1.png)



```R
pdf("~/Desktop/WGS_scRNA_corr_avgSample.pdf", height = 4, width = 4)

smoothScatter(x = combined$scRNA_Value, 
              y = combined$log2,
              xlab = "InferCNV Score (scRNA-seq)",
              ylab = "Log2 Ratio (WGS)",
              main = "r=0.68, p < 2.2e-16",
              cex.axis = 0.8
             )
abline(lm(combined$log2 ~ combined$scRNA_Value), 
       lwd = 2,
       col = "red"
      )


smoothScatter(x = combined$scRNA_Value, 
              y = combined$GISTIC,
              xlab = "InferCNV Score (scRNA-seq)",
              ylab = "GISTIC (WGS)",
              main = "r=0.62, p < 2.2e-16",
               cex.axis = 0.8
             )
abline(lm(combined$GISTIC ~ combined$scRNA_Value), 
       lwd = 2,
       col = "red"
      )

dev.off()
```


<strong>pdf:</strong> 2



```R
combined$GISTIC <- factor(combined$GISTIC)

p_meds <- ddply(combined, .(GISTIC), summarise, med = median(scRNA_Value))
p_meds

p <- ggplot(combined, aes(x=GISTIC, y=scRNA_Value, color = GISTIC)) + geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
scale_color_manual(values=c("darkblue", "darkblue", "black", "darkred", "darkred")) + theme_classic() 

p
```


<table>
<thead><tr><th scope=col>GISTIC</th><th scope=col>med</th></tr></thead>
<tbody>
	<tr><td>-2          </td><td>-0.126797012</td></tr>
	<tr><td>-1          </td><td>-0.154560723</td></tr>
	<tr><td>0           </td><td>-0.001793286</td></tr>
	<tr><td>1           </td><td> 0.179210290</td></tr>
	<tr><td>2           </td><td> 0.245607426</td></tr>
</tbody>
</table>






![png](output_37_2.png)



```R

q_meds <- ddply(combined, .(Category), summarise, med = median(scRNA_Value))
q_meds

q <- ggplot(combined, aes(x=Category, y=scRNA_Value, color = Category)) + geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
scale_color_manual(values=c("darkblue", "black", "darkred")) + theme_classic() +
 geom_hline(yintercept = q_meds$med[c(1,3)], lty = 2, col = "grey") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        legend.position = "none"
       ) + xlab("") + ylab("")
q
```


<table>
<thead><tr><th scope=col>Category</th><th scope=col>med</th></tr></thead>
<tbody>
	<tr><td>Deletion    </td><td>-0.153596420</td></tr>
	<tr><td>Gain        </td><td> 0.182539324</td></tr>
	<tr><td>Neutral     </td><td>-0.001793286</td></tr>
</tbody>
</table>






![png](output_38_2.png)

