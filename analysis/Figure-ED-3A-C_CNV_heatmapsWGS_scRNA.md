
---
# Plot scRNA and WGS CNV calls
----

Visualize CNV tracks of scRNA and WGS for samples profiled with both technologies.

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/WGS_scRNA_plotting


```R
library(pheatmap)
```

    Warning message:
    “package ‘pheatmap’ was built under R version 3.4.4”

---
## 1.0 Isolate and format CNV calls from scRNA and WGS
---

#this file needs an index plot
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/Matched_scRNA_WGS_CNV.Rdata")
head(combined)

sample.avg <- readRDS("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/defineCutoffs/BTSC_SampleAvgCNVs.rds")

###load in cluster averaged CNVs
### extract only those with WGS data

load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/GlobalBTSC_CNVs_averge.RData")
samples <- names(table(combined$scRNA_sample))
cols <- c()
for (i in 1:length(samples)){
    
    cols <- c(cols, grep(samples[i], colnames(avg.cnv.df)))
    
}

cluster.cnv <- avg.cnv.df[,cols]
head(cluster.cnv)

##add in the sample averaged CNVs
sample.cnv <- sample.avg[, samples]
head(sample.cnv)

cnv <- cbind(cluster.cnv, sample.cnv)

saveRDS(cnv, file = "Cluster_Sample_CNVs.rds")
saveRDS(CNV.genes, file = "scRNA_CNV_GenePos.rds")## now get the WGS GISTIC data

gistic.file <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_thresholded.by_genes.txt"
gistic <- read.table(gistic.file, header = T, sep = "\t")
gistic <- gistic[, c("Gene.Symbol", "Locus.ID", "Cytoband", samples)]
head(gistic)
dim(gistic)

saveRDS(gistic, file = "WGS_GISTIC.rds")

## now get the WGS log2 data

seg <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/T_L_Comparison/GISTIC/output/all_data_by_genes.txt"
seg <- read.table(seg, header = T, sep = "\t")
seg <- seg[, c("Gene.Symbol", "Cytoband", samples)]
head(seg)
dim(seg)

saveRDS(seg, file = "WGS_log2_genes.rds")
---
# 2.0 Plot
---



```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/WGS_scRNA_plotting")
```


```R
options(stringsAsFactors = FALSE)

scRNA <- readRDS("Cluster_Sample_CNVs.rds")
#head(scRNA)

gistic <- readRDS("WGS_GISTIC.rds")
#head(gistic)

log2 <- readRDS("WGS_log2_genes.rds")
#head(log2)

gene.pos <- readRDS("scRNA_CNV_GenePos.rds")
#head(gene.pos)

CNV.genes <- readRDS("scRNA_CNV_GenePos.rds")

##################
common.genes <- rownames(CNV.genes)[rownames(CNV.genes) %in% log2$Gene.Symbol]
length(common.genes) ## 6241 gens

##susbet all the data
CNV.genes$Gene <- rownames(CNV.genes)
CNV.genes <- CNV.genes[CNV.genes$Gene %in% common.genes, ]
CNV.genes$Gene <- NULL
#head(CNV.genes)
#dim(CNV.genes)



log2 <- log2[log2$Gene.Symbol %in% common.genes, ]
rownames(log2) <- log2$Gene.Symbol
log2 <- log2[common.genes, ]
log2 <- log2[ ,-c(1:2)]
colnames(log2) <- paste0(colnames(log2), "_WGS")
#head(log2)
#dim(log2)

scRNA <- scRNA[rownames(scRNA) %in% common.genes, ]
scRNA <- scRNA[common.genes, ]
#dim(scRNA)
#head(scRNA)

################
dat <- cbind(scRNA, log2)

colnames(dat)[grep("L$", colnames(dat))] <- paste0(colnames(dat)[grep("L$", colnames(dat))], "_SampleAvg")

##### reorder dat
order <- read.csv("~/Desktop/CNV_order.csv")
order <- order$Sample

dat <- dat[ ,order]


```


6241



```R
dat[dat > 1] <- 1
dat[dat < -1] <- -1
head(dat)

max(dat)
```


<table>
<thead><tr><th></th><th scope=col>BT147_L_C1</th><th scope=col>BT147_L_C2</th><th scope=col>BT147_L_SampleAvg</th><th scope=col>BT147_L_WGS</th><th scope=col>BT67_L_C1</th><th scope=col>BT67_L_C2</th><th scope=col>BT67_L_SampleAvg</th><th scope=col>BT67_L_WGS</th><th scope=col>BT73_L_C1</th><th scope=col>BT73_L_C2</th><th scope=col>⋯</th><th scope=col>G885_L_C1</th><th scope=col>G885_L_C2</th><th scope=col>G885_L_C3</th><th scope=col>G885_L_SampleAvg</th><th scope=col>G885_L_WGS</th><th scope=col>G895_L_C1</th><th scope=col>G895_L_C2</th><th scope=col>G895_L_C3</th><th scope=col>G895_L_SampleAvg</th><th scope=col>G895_L_WGS</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>0.03504792   </td><td>0.01399746   </td><td>0.02513320   </td><td>-0.931       </td><td>-0.10570094  </td><td>-0.09201973  </td><td>-0.09986196  </td><td>0.337        </td><td>-0.03080463  </td><td> 0.0008629272</td><td>⋯            </td><td>0.07449602   </td><td>0.06592014   </td><td>0.05022230   </td><td>0.06912226   </td><td>0.061        </td><td>0.07666680   </td><td>0.03863795   </td><td>0.04401388   </td><td>0.06950744   </td><td>-0.002       </td></tr>
	<tr><th scope=row>TMEM201</th><td>0.03990202   </td><td>0.01693410   </td><td>0.02908419   </td><td>-0.931       </td><td>-0.09366095  </td><td>-0.08518207  </td><td>-0.09004226  </td><td>0.337        </td><td>-0.02512797  </td><td> 0.0048454969</td><td>⋯            </td><td>0.09724205   </td><td>0.08170858   </td><td>0.06610565   </td><td>0.08883081   </td><td>0.061        </td><td>0.08619665   </td><td>0.04838392   </td><td>0.04903777   </td><td>0.07883048   </td><td>-0.002       </td></tr>
	<tr><th scope=row>CLSTN1</th><td>0.03787451   </td><td>0.01187433   </td><td>0.02562848   </td><td>-0.931       </td><td>-0.09216101  </td><td>-0.08490364  </td><td>-0.08906365  </td><td>0.291        </td><td>-0.02563131  </td><td> 0.0081523679</td><td>⋯            </td><td>0.09649358   </td><td>0.08311512   </td><td>0.06713228   </td><td>0.08898703   </td><td>0.061        </td><td>0.09086135   </td><td>0.05229157   </td><td>0.05918094   </td><td>0.08367595   </td><td>-0.002       </td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>0.03495229   </td><td>0.00805269   </td><td>0.02228264   </td><td>-0.931       </td><td>-0.08319649  </td><td>-0.07528557  </td><td>-0.07982020  </td><td>0.291        </td><td>-0.03290911  </td><td>-0.0095088953</td><td>⋯            </td><td>0.09763405   </td><td>0.08374670   </td><td>0.06881378   </td><td>0.09001282   </td><td>0.061        </td><td>0.09058943   </td><td>0.05127739   </td><td>0.06141334   </td><td>0.08343003   </td><td>-0.002       </td></tr>
	<tr><th scope=row>LZIC</th><td>0.04058198   </td><td>0.01523971   </td><td>0.02864583   </td><td>-0.931       </td><td>-0.07266495  </td><td>-0.06440919  </td><td>-0.06914149  </td><td>0.291        </td><td>-0.03097119  </td><td>-0.0223801521</td><td>⋯            </td><td>0.10958090   </td><td>0.10042074   </td><td>0.08119290   </td><td>0.10358747   </td><td>0.061        </td><td>0.10344303   </td><td>0.06084255   </td><td>0.07702164   </td><td>0.09595880   </td><td>-0.002       </td></tr>
	<tr><th scope=row>NMNAT1</th><td>0.04250626   </td><td>0.01571912   </td><td>0.02988958   </td><td>-0.953       </td><td>-0.05806666  </td><td>-0.05709092  </td><td>-0.05765023  </td><td>0.291        </td><td>-0.02500879  </td><td>-0.0133474736</td><td>⋯            </td><td>0.11480409   </td><td>0.10324650   </td><td>0.08167909   </td><td>0.10751953   </td><td>0.061        </td><td>0.09406364   </td><td>0.05282940   </td><td>0.07116908   </td><td>0.08696079   </td><td>-0.002       </td></tr>
</tbody>
</table>




1



```R
### order by number of clusters

col.order <- read.csv("~/Desktop/scRNA_WGS_order_cluster.csv")
head(col.order)

order <- col.order[order(col.order$Clusters, col.order$Sample), ]
order <- order$Sample
order


```


<table>
<thead><tr><th scope=col>Sample</th><th scope=col>Clusters</th></tr></thead>
<tbody>
	<tr><td>BT147_L_C1       </td><td>2                </td></tr>
	<tr><td>BT147_L_C2       </td><td>2                </td></tr>
	<tr><td>BT147_L_SampleAvg</td><td>2                </td></tr>
	<tr><td>BT147_L_WGS      </td><td>2                </td></tr>
	<tr><td>BT67_L_C1        </td><td>2                </td></tr>
	<tr><td>BT67_L_C2        </td><td>2                </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'BT147_L_C1'</li>
	<li>'BT147_L_C2'</li>
	<li>'BT147_L_SampleAvg'</li>
	<li>'BT147_L_WGS'</li>
	<li>'BT67_L_C1'</li>
	<li>'BT67_L_C2'</li>
	<li>'BT67_L_SampleAvg'</li>
	<li>'BT67_L_WGS'</li>
	<li>'BT73_L_C1'</li>
	<li>'BT73_L_C2'</li>
	<li>'BT73_L_SampleAvg'</li>
	<li>'BT73_L_WGS'</li>
	<li>'BT89_L_C1'</li>
	<li>'BT89_L_C2'</li>
	<li>'BT89_L_SampleAvg'</li>
	<li>'BT89_L_WGS'</li>
	<li>'BT94_L_C1'</li>
	<li>'BT94_L_C2'</li>
	<li>'BT94_L_SampleAvg'</li>
	<li>'BT94_L_WGS'</li>
	<li>'G549_L_C1'</li>
	<li>'G549_L_C2'</li>
	<li>'G549_L_SampleAvg'</li>
	<li>'G549_L_WGS'</li>
	<li>'G564_L_C1'</li>
	<li>'G564_L_C2'</li>
	<li>'G564_L_SampleAvg'</li>
	<li>'G564_L_WGS'</li>
	<li>'G729_L_C1'</li>
	<li>'G729_L_C2'</li>
	<li>'G729_L_SampleAvg'</li>
	<li>'G729_L_WGS'</li>
	<li>'G851_L_C1'</li>
	<li>'G851_L_C2'</li>
	<li>'G851_L_SampleAvg'</li>
	<li>'G851_L_WGS'</li>
	<li>'BT84_L_C1'</li>
	<li>'BT84_L_C2'</li>
	<li>'BT84_L_C3'</li>
	<li>'BT84_L_SampleAvg'</li>
	<li>'BT84_L_WGS'</li>
	<li>'G523_L_C1'</li>
	<li>'G523_L_C2'</li>
	<li>'G523_L_C3'</li>
	<li>'G523_L_SampleAvg'</li>
	<li>'G523_L_WGS'</li>
	<li>'G566_L_C1'</li>
	<li>'G566_L_C2'</li>
	<li>'G566_L_C3'</li>
	<li>'G566_L_SampleAvg'</li>
	<li>'G566_L_WGS'</li>
	<li>'G799_L_C1'</li>
	<li>'G799_L_C2'</li>
	<li>'G799_L_C3'</li>
	<li>'G799_L_SampleAvg'</li>
	<li>'G799_L_WGS'</li>
	<li>'G800_L_C1'</li>
	<li>'G800_L_C2'</li>
	<li>'G800_L_C3'</li>
	<li>'G800_L_SampleAvg'</li>
	<li>'G800_L_WGS'</li>
	<li>'G876_L_C1'</li>
	<li>'G876_L_C2'</li>
	<li>'G876_L_C3'</li>
	<li>'G876_L_SampleAvg'</li>
	<li>'G876_L_WGS'</li>
	<li>'G885_L_C1'</li>
	<li>'G885_L_C2'</li>
	<li>'G885_L_C3'</li>
	<li>'G885_L_SampleAvg'</li>
	<li>'G885_L_WGS'</li>
	<li>'G895_L_C1'</li>
	<li>'G895_L_C2'</li>
	<li>'G895_L_C3'</li>
	<li>'G895_L_SampleAvg'</li>
	<li>'G895_L_WGS'</li>
	<li>'G583_L_C1'</li>
	<li>'G583_L_C2'</li>
	<li>'G583_L_C3'</li>
	<li>'G583_L_C4'</li>
	<li>'G583_L_SampleAvg'</li>
	<li>'G583_L_WGS'</li>
	<li>'G637_L_C1'</li>
	<li>'G637_L_C2'</li>
	<li>'G637_L_C3'</li>
	<li>'G637_L_C4'</li>
	<li>'G637_L_C5'</li>
	<li>'G637_L_SampleAvg'</li>
	<li>'G637_L_WGS'</li>
	<li>'G620_L_C1'</li>
	<li>'G620_L_C2'</li>
	<li>'G620_L_C3'</li>
	<li>'G620_L_C4'</li>
	<li>'G620_L_C5'</li>
	<li>'G620_L_C6'</li>
	<li>'G620_L_SampleAvg'</li>
	<li>'G620_L_WGS'</li>
</ol>




```R
dat <- dat[ ,order]
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>BT147_L_C1</th><th scope=col>BT147_L_C2</th><th scope=col>BT147_L_SampleAvg</th><th scope=col>BT147_L_WGS</th><th scope=col>BT67_L_C1</th><th scope=col>BT67_L_C2</th><th scope=col>BT67_L_SampleAvg</th><th scope=col>BT67_L_WGS</th><th scope=col>BT73_L_C1</th><th scope=col>BT73_L_C2</th><th scope=col>⋯</th><th scope=col>G637_L_SampleAvg</th><th scope=col>G637_L_WGS</th><th scope=col>G620_L_C1</th><th scope=col>G620_L_C2</th><th scope=col>G620_L_C3</th><th scope=col>G620_L_C4</th><th scope=col>G620_L_C5</th><th scope=col>G620_L_C6</th><th scope=col>G620_L_SampleAvg</th><th scope=col>G620_L_WGS</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>0.03504792   </td><td>0.01399746   </td><td>0.02513320   </td><td>-0.931       </td><td>-0.10570094  </td><td>-0.09201973  </td><td>-0.09986196  </td><td>0.337        </td><td>-0.03080463  </td><td> 0.0008629272</td><td>⋯            </td><td>-0.1290687   </td><td>-0.021       </td><td>0.1066676    </td><td>0.09700489   </td><td>0.09891157   </td><td>-0.2660879   </td><td>0.08977328   </td><td>0.03127273   </td><td>0.07282545   </td><td>-0.145       </td></tr>
	<tr><th scope=row>TMEM201</th><td>0.03990202   </td><td>0.01693410   </td><td>0.02908419   </td><td>-0.931       </td><td>-0.09366095  </td><td>-0.08518207  </td><td>-0.09004226  </td><td>0.337        </td><td>-0.02512797  </td><td> 0.0048454969</td><td>⋯            </td><td>-0.1194177   </td><td>-0.021       </td><td>0.1213500    </td><td>0.10876507   </td><td>0.11923526   </td><td>-0.2572933   </td><td>0.10926786   </td><td>0.04352081   </td><td>0.08724775   </td><td>-0.145       </td></tr>
	<tr><th scope=row>CLSTN1</th><td>0.03787451   </td><td>0.01187433   </td><td>0.02562848   </td><td>-0.931       </td><td>-0.09216101  </td><td>-0.08490364  </td><td>-0.08906365  </td><td>0.291        </td><td>-0.02563131  </td><td> 0.0081523679</td><td>⋯            </td><td>-0.1227548   </td><td>-0.021       </td><td>0.1246295    </td><td>0.10926400   </td><td>0.12556332   </td><td>-0.2634345   </td><td>0.11102832   </td><td>0.04471645   </td><td>0.08930406   </td><td>-0.145       </td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>0.03495229   </td><td>0.00805269   </td><td>0.02228264   </td><td>-0.931       </td><td>-0.08319649  </td><td>-0.07528557  </td><td>-0.07982020  </td><td>0.291        </td><td>-0.03290911  </td><td>-0.0095088953</td><td>⋯            </td><td>-0.1213200   </td><td>-0.021       </td><td>0.1273295    </td><td>0.11228151   </td><td>0.12672267   </td><td>-0.2605796   </td><td>0.11150668   </td><td>0.04864290   </td><td>0.09178882   </td><td>-0.145       </td></tr>
	<tr><th scope=row>LZIC</th><td>0.04058198   </td><td>0.01523971   </td><td>0.02864583   </td><td>-0.931       </td><td>-0.07266495  </td><td>-0.06440919  </td><td>-0.06914149  </td><td>0.291        </td><td>-0.03097119  </td><td>-0.0223801521</td><td>⋯            </td><td>-0.1248144   </td><td>-0.021       </td><td>0.1354603    </td><td>0.11582118   </td><td>0.12833345   </td><td>-0.2605023   </td><td>0.11599754   </td><td>0.03807635   </td><td>0.09581982   </td><td>-0.145       </td></tr>
	<tr><th scope=row>NMNAT1</th><td>0.04250626   </td><td>0.01571912   </td><td>0.02988958   </td><td>-0.953       </td><td>-0.05806666  </td><td>-0.05709092  </td><td>-0.05765023  </td><td>0.291        </td><td>-0.02500879  </td><td>-0.0133474736</td><td>⋯            </td><td>-0.1216510   </td><td>-0.021       </td><td>0.1315208    </td><td>0.10699000   </td><td>0.12880729   </td><td>-0.2593459   </td><td>0.11700634   </td><td>0.03443788   </td><td>0.09162905   </td><td>-0.145       </td></tr>
</tbody>
</table>




```R



mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11))
                   )

    names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
    #mat_colors

#color scale for CNVs
cols <- c("#67001f",
          "#b2182b",
          "#d6604d",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          "#d1e5f0",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 1, by = 0.01)
cnv.cols <- colorRampPalette(cols)(length(breaksList))
#cnv.cols <- colfunc(50)

file.name <- paste0("~/Desktop/AllSampels_scRNA_WGS_ClusterOrder.png")

pheatmap(t(dat),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         cellheight = 10,
         annotation_col = CNV.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 16,
         width = 15,
         filename = file.name
        )
```

----
## plot each sample


```R
samples <- colnames(gistic)[grep("_L", colnames(gistic))]

for (i in 1:length(samples)){
    


#i <- 1

sample <- samples[i]
sample

subset <- dat[, grep(sample, colnames(dat)) ]

subset[subset > 1] <- 1
subset[subset < -1] <- -1
head(subset)

#subset <- dat[, grep(sample, colnames(dat)) ]
#subset <- subset[, -grep("WGS", colnames(subset))]
#head(subset)

#wgs.subset <- dat[, grep(sample, colnames(dat)) ]
#wgs.subset <- data.frame(wgs.subset[, grep("WGS", colnames(wgs.subset))])
#colnames(wgs.subset) <- paste0(sample, "_WGS")
#rownames(wgs.subset) <- rownames(subset)
#head(wgs.subset)
    
    ### PLOT SINGLE CELL

mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11))
                   )

    names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
    #mat_colors

#color scale for CNVs
cols <- c("#67001f",
          "#b2182b",
          "#d6604d",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          "#d1e5f0",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 1, by = 0.01)
cnv.cols <- colorRampPalette(cols)(length(breaksList))
#cnv.cols <- colfunc(50)

file.name <- paste0("~/Desktop/", sample, "_scRNA_WGS.jpeg")

pheatmap(t(subset),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         cellheight = 15,
         annotation_col = CNV.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 5,
         width = 10,
         filename = file.name
        )
}
```


```R
cat(colnames(dat), sep = "\n")
```

    BT147_L_C1
    BT147_L_C2
    BT67_L_C2
    BT67_L_C1
    BT73_L_C1
    BT73_L_C2
    BT84_L_C1
    BT84_L_C2
    BT84_L_C3
    BT89_L_C2
    BT89_L_C1
    BT94_L_C2
    BT94_L_C1
    G523_L_C1
    G523_L_C2
    G523_L_C3
    G549_L_C1
    G549_L_C2
    G564_L_C2
    G564_L_C1
    G566_L_C1
    G566_L_C2
    G566_L_C3
    G583_L_C1
    G583_L_C3
    G583_L_C2
    G583_L_C4
    G620_L_C1
    G620_L_C5
    G620_L_C2
    G620_L_C3
    G620_L_C6
    G620_L_C4
    G637_L_C3
    G637_L_C1
    G637_L_C2
    G637_L_C4
    G637_L_C5
    G729_L_C2
    G729_L_C1
    G799_L_C1
    G799_L_C3
    G799_L_C2
    G800_L_C1
    G800_L_C2
    G800_L_C3
    G851_L_C1
    G851_L_C2
    G876_L_C1
    G876_L_C2
    G876_L_C3
    G885_L_C2
    G885_L_C3
    G885_L_C1
    G895_L_C1
    G895_L_C2
    G895_L_C3
    BT147_L
    BT67_L
    BT73_L
    BT84_L
    BT89_L
    BT94_L
    G523_L
    G549_L
    G564_L
    G566_L
    G583_L
    G620_L
    G637_L
    G729_L
    G799_L
    G800_L
    G851_L
    G876_L
    G885_L
    G895_L
    BT147_L_WGS
    BT67_L_WGS
    BT73_L_WGS
    BT84_L_WGS
    BT89_L_WGS
    BT94_L_WGS
    G523_L_WGS
    G549_L_WGS
    G564_L_WGS
    G566_L_WGS
    G583_L_WGS
    G620_L_WGS
    G637_L_WGS
    G729_L_WGS
    G799_L_WGS
    G800_L_WGS
    G851_L_WGS
    G876_L_WGS
    G885_L_WGS
    G895_L_WGS



```R
colnames(dat)[grep("L$", colnames(dat))] <- paste0(colnames(dat)[grep("L$", colnames(dat))], "_SampleAvg")
```


<ol class=list-inline>
	<li>'BT147_L'</li>
	<li>'BT67_L'</li>
	<li>'BT73_L'</li>
	<li>'BT84_L'</li>
	<li>'BT89_L'</li>
	<li>'BT94_L'</li>
	<li>'G523_L'</li>
	<li>'G549_L'</li>
	<li>'G564_L'</li>
	<li>'G566_L'</li>
	<li>'G583_L'</li>
	<li>'G620_L'</li>
	<li>'G637_L'</li>
	<li>'G729_L'</li>
	<li>'G799_L'</li>
	<li>'G800_L'</li>
	<li>'G851_L'</li>
	<li>'G876_L'</li>
	<li>'G885_L'</li>
	<li>'G895_L'</li>
</ol>




```R

```


```R

```


```R

```


```R
### PLOT SINGLE CELL

mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11))
                   )

    names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
    mat_colors

#color scale for CNVs
cols <- c("#67001f",
          "#b2182b",
          "#d6604d",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 1, by = 0.01)
cnv.cols <- colorRampPalette(cols)(length(breaksList))
#cnv.cols <- colfunc(50)


pheatmap(t(subset),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = CNV.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 15,
         width = 25
        # filename = "~/Desktop/scRNA_CNV.jpeg"
        )
#### PLOT WGS

cols <- c("#67001f",
          "#67001f",
          "#67001f",
          "#67001f",
          "#b2182b",
          "#b2182b",
          "#b2182b",
          "#d6604d",
          "#d6604d",
          "#d6604d",
          "#f4a582",
          "#f4a582",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          #white",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 3, by = 0.001)
cnv.cols <- colorRampPalette(cols)(length(breaksList))

pheatmap(t(wgs.subset),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = CNV.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 15,
         width = 25
        # filename = "~/Desktop/scRNA_CNV.jpeg"
        )
```


<strong>$Chromosome</strong> = <dl class=dl-horizontal>
	<dt>1</dt>
		<dd>'lightgrey'</dd>
	<dt>2</dt>
		<dd>'darkgrey'</dd>
	<dt>3</dt>
		<dd>'lightgrey'</dd>
	<dt>4</dt>
		<dd>'darkgrey'</dd>
	<dt>5</dt>
		<dd>'lightgrey'</dd>
	<dt>6</dt>
		<dd>'darkgrey'</dd>
	<dt>7</dt>
		<dd>'lightgrey'</dd>
	<dt>8</dt>
		<dd>'darkgrey'</dd>
	<dt>9</dt>
		<dd>'lightgrey'</dd>
	<dt>10</dt>
		<dd>'darkgrey'</dd>
	<dt>11</dt>
		<dd>'lightgrey'</dd>
	<dt>12</dt>
		<dd>'darkgrey'</dd>
	<dt>13</dt>
		<dd>'lightgrey'</dd>
	<dt>14</dt>
		<dd>'darkgrey'</dd>
	<dt>15</dt>
		<dd>'lightgrey'</dd>
	<dt>16</dt>
		<dd>'darkgrey'</dd>
	<dt>17</dt>
		<dd>'lightgrey'</dd>
	<dt>18</dt>
		<dd>'darkgrey'</dd>
	<dt>19</dt>
		<dd>'lightgrey'</dd>
	<dt>20</dt>
		<dd>'darkgrey'</dd>
	<dt>21</dt>
		<dd>'lightgrey'</dd>
	<dt>22</dt>
		<dd>'darkgrey'</dd>
</dl>




![png](output_19_1.png)



![png](output_19_2.png)


-----
## WGS


```R
log2 <- readRDS("WGS_log2_genes.rds")
head(log2)

```


<table>
<thead><tr><th scope=col>Gene.Symbol</th><th scope=col>Cytoband</th><th scope=col>BT147_L</th><th scope=col>BT67_L</th><th scope=col>BT73_L</th><th scope=col>BT84_L</th><th scope=col>BT89_L</th><th scope=col>BT94_L</th><th scope=col>G523_L</th><th scope=col>G549_L</th><th scope=col>⋯</th><th scope=col>G583_L</th><th scope=col>G620_L</th><th scope=col>G637_L</th><th scope=col>G729_L</th><th scope=col>G799_L</th><th scope=col>G800_L</th><th scope=col>G851_L</th><th scope=col>G876_L</th><th scope=col>G885_L</th><th scope=col>G895_L</th></tr></thead>
<tbody>
	<tr><td>DDX11L1        </td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
	<tr><td>MIR6859-1|chr1 </td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
	<tr><td>MIR6859-2|chr1 </td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
	<tr><td>MIR6859-3|chr1 </td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
	<tr><td>WASH7P         </td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
	<tr><td>MIR1302-10|chr1</td><td>1p36.33        </td><td>0.068          </td><td>0.416          </td><td>-0.271         </td><td>0.271          </td><td>0.057          </td><td>-0.038         </td><td>0.258          </td><td>0.068          </td><td>⋯              </td><td>0.423          </td><td>0.529          </td><td>0.328          </td><td>0.193          </td><td>-0.081         </td><td>0.049          </td><td>0.674          </td><td>1.303          </td><td>0.608          </td><td>0.087          </td></tr>
</tbody>
</table>




```R
log2$Cytoband <- as.character(log2$Cytoband)

log2$Cytoband[grep("^1p", log2$Cytoband)] <- "chr1"
log2$Cytoband[grep("^1q", log2$Cytoband)] <- "chr1"

log2$Cytoband[grep("^2p", log2$Cytoband)] <- "chr2"
log2$Cytoband[grep("^2q", log2$Cytoband)] <- "chr2"

log2$Cytoband[grep("^3p", log2$Cytoband)] <- "chr3"
log2$Cytoband[grep("^3q", log2$Cytoband)] <- "chr3"

log2$Cytoband[grep("^4p", log2$Cytoband)] <- "chr4"
log2$Cytoband[grep("^4q", log2$Cytoband)] <- "chr4"

log2$Cytoband[grep("^5p", log2$Cytoband)] <- "chr5"
log2$Cytoband[grep("^5q", log2$Cytoband)] <- "chr5"

log2$Cytoband[grep("^6p", log2$Cytoband)] <- "chr6"
log2$Cytoband[grep("^6q", log2$Cytoband)] <- "chr6"

log2$Cytoband[grep("^7p", log2$Cytoband)] <- "chr7"
log2$Cytoband[grep("^7q", log2$Cytoband)] <- "chr7"

log2$Cytoband[grep("^8p", log2$Cytoband)] <- "chr8"
log2$Cytoband[grep("^8q", log2$Cytoband)] <- "chr8"

log2$Cytoband[grep("^9p", log2$Cytoband)] <- "chr9"
log2$Cytoband[grep("^9q", log2$Cytoband)] <- "chr9"

log2$Cytoband[grep("^10p", log2$Cytoband)] <- "chr10"
log2$Cytoband[grep("^10q", log2$Cytoband)] <- "chr10"

log2$Cytoband[grep("^11p", log2$Cytoband)] <- "chr11"
log2$Cytoband[grep("^11q", log2$Cytoband)] <- "chr11"

log2$Cytoband[grep("^12p", log2$Cytoband)] <- "chr12"
log2$Cytoband[grep("^12q", log2$Cytoband)] <- "chr12"

log2$Cytoband[grep("^13p", log2$Cytoband)] <- "chr13"
log2$Cytoband[grep("^13q", log2$Cytoband)] <- "chr13"

log2$Cytoband[grep("^14p", log2$Cytoband)] <- "chr14"
log2$Cytoband[grep("^14q", log2$Cytoband)] <- "chr14"

log2$Cytoband[grep("^15p", log2$Cytoband)] <- "chr15"
log2$Cytoband[grep("^15q", log2$Cytoband)] <- "chr15"

log2$Cytoband[grep("^16p", log2$Cytoband)] <- "chr16"
log2$Cytoband[grep("^16q", log2$Cytoband)] <- "chr16"

log2$Cytoband[grep("^17p", log2$Cytoband)] <- "chr17"
log2$Cytoband[grep("^17q", log2$Cytoband)] <- "chr17"

log2$Cytoband[grep("^18p", log2$Cytoband)] <- "chr18"
log2$Cytoband[grep("^18q", log2$Cytoband)] <- "chr18"

log2$Cytoband[grep("^19p", log2$Cytoband)] <- "chr19"
log2$Cytoband[grep("^19q", log2$Cytoband)] <- "chr19"

log2$Cytoband[grep("^20p", log2$Cytoband)] <- "chr20"
log2$Cytoband[grep("^20q", log2$Cytoband)] <- "chr20"

log2$Cytoband[grep("^21p", log2$Cytoband)] <- "chr21"
log2$Cytoband[grep("^21q", log2$Cytoband)] <- "chr21"

log2$Cytoband[grep("^22p", log2$Cytoband)] <- "chr22"
log2$Cytoband[grep("^22q", log2$Cytoband)] <- "chr22"

table(log2$Cytoband)

rownames(log2) <- log2(Cytoband)
```


    
     chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
     2699  1076  1635  1341   605   902   944  1095  1526   410  1750  1741   761 
    chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
      370   617  1496  1026  1250  1400  1283   976  1085 



    Error in eval(expr, envir, enclos): object 'Cytoband' not found
    Traceback:




```R
bb <- log2[,-c(1:2)]
head(bb)

rownames(bb) <- log2$Gene.Symbol

wgs.genes <- data.frame(log2$Cytoband)
colnames(wgs.genes) <- "Chromosome"
rownames(wgs.genes) <- log2$Gene.Symbol
head(wgs.genes)


cols <- c("#67001f",
          "#67001f",
          "#67001f",
          "#b2182b",
          "#b2182b",
          "#b2182b",
          "#d6604d",
          "#d6604d",
          "#d6604d",
          "#f4a582",
          "#f4a582",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          #"white",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 1, by = 0.01)
cnv.cols <- colorRampPalette(cols)(length(breaksList))

mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11))
                   )

    names(mat_colors$Chromosome) <- unique(log2$Cytoband)
    mat_colors


pheatmap(t(bb),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = wgs.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 7,
         width = 25,
         filename = "~/Desktop/WGS.jpeg"
        )
```


<table>
<thead><tr><th scope=col>BT147_L</th><th scope=col>BT67_L</th><th scope=col>BT73_L</th><th scope=col>BT84_L</th><th scope=col>BT89_L</th><th scope=col>BT94_L</th><th scope=col>G523_L</th><th scope=col>G549_L</th><th scope=col>G564_L</th><th scope=col>G566_L</th><th scope=col>G583_L</th><th scope=col>G620_L</th><th scope=col>G637_L</th><th scope=col>G729_L</th><th scope=col>G799_L</th><th scope=col>G800_L</th><th scope=col>G851_L</th><th scope=col>G876_L</th><th scope=col>G885_L</th><th scope=col>G895_L</th></tr></thead>
<tbody>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
	<tr><td>0.068 </td><td>0.416 </td><td>-0.271</td><td>0.271 </td><td>0.057 </td><td>-0.038</td><td>0.258 </td><td>0.068 </td><td>-0.062</td><td>0.278 </td><td>0.423 </td><td>0.529 </td><td>0.328 </td><td>0.193 </td><td>-0.081</td><td>0.049 </td><td>0.674 </td><td>1.303 </td><td>0.608 </td><td>0.087 </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>Chromosome</th></tr></thead>
<tbody>
	<tr><th scope=row>DDX11L1</th><td>chr1</td></tr>
	<tr><th scope=row>MIR6859-1|chr1</th><td>chr1</td></tr>
	<tr><th scope=row>MIR6859-2|chr1</th><td>chr1</td></tr>
	<tr><th scope=row>MIR6859-3|chr1</th><td>chr1</td></tr>
	<tr><th scope=row>WASH7P</th><td>chr1</td></tr>
	<tr><th scope=row>MIR1302-10|chr1</th><td>chr1</td></tr>
</tbody>
</table>




<strong>$Chromosome</strong> = <dl class=dl-horizontal>
	<dt>chr1</dt>
		<dd>'lightgrey'</dd>
	<dt>chr2</dt>
		<dd>'darkgrey'</dd>
	<dt>chr3</dt>
		<dd>'lightgrey'</dd>
	<dt>chr4</dt>
		<dd>'darkgrey'</dd>
	<dt>chr5</dt>
		<dd>'lightgrey'</dd>
	<dt>chr6</dt>
		<dd>'darkgrey'</dd>
	<dt>chr7</dt>
		<dd>'lightgrey'</dd>
	<dt>chr8</dt>
		<dd>'darkgrey'</dd>
	<dt>chr9</dt>
		<dd>'lightgrey'</dd>
	<dt>chr10</dt>
		<dd>'darkgrey'</dd>
	<dt>chr11</dt>
		<dd>'lightgrey'</dd>
	<dt>chr12</dt>
		<dd>'darkgrey'</dd>
	<dt>chr13</dt>
		<dd>'lightgrey'</dd>
	<dt>chr14</dt>
		<dd>'darkgrey'</dd>
	<dt>chr15</dt>
		<dd>'lightgrey'</dd>
	<dt>chr16</dt>
		<dd>'darkgrey'</dd>
	<dt>chr17</dt>
		<dd>'lightgrey'</dd>
	<dt>chr18</dt>
		<dd>'darkgrey'</dd>
	<dt>chr19</dt>
		<dd>'lightgrey'</dd>
	<dt>chr20</dt>
		<dd>'darkgrey'</dd>
	<dt>chr21</dt>
		<dd>'lightgrey'</dd>
	<dt>chr22</dt>
		<dd>'darkgrey'</dd>
</dl>



-----
## WGS at gene level


```R
bb <- log2[,-c(1:2)]
head(bb)

rownames(bb) <- log2$Gene.Symbol

wgs.genes <- data.frame(log2$Cytoband)
colnames(wgs.genes) <- "Chromosome"
rownames(wgs.genes) <- log2$Gene.Symbol
head(wgs.genes)


cols <- c("#67001f",
          "#67001f",
          "#67001f",
          "#b2182b",
          "#b2182b",
          "#b2182b",
          "#d6604d",
          "#d6604d",
          "#d6604d",
          "#f4a582",
          "#f4a582",
          "#f4a582",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "#fddbc7",
          "white",
          "white",
          "white",
          #"white",
          "#d1e5f0",
          "#92c5de",
          "#4393c3",
          "#2166ac",
          "#053061"
          )
breaksList = seq(-1, 1, by = 0.01)
cnv.cols <- colorRampPalette(cols)(length(breaksList))

mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11))
                   )

    names(mat_colors$Chromosome) <- unique(log2$Cytoband)
    mat_colors


pheatmap(t(bb),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = wgs.genes,
         #annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 7,
         width = 25,
         filename = "~/Desktop/WGS.jpeg"
        )
```
