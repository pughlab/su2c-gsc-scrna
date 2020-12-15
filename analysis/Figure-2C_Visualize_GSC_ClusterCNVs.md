
---
# 1.0 BTSC cohort wide CNV profiles
---
L.Richards  
Need to visualize all the samples together on one plot. To do this, collapse all the CNV profiles by transcriptional clusters. Need a similarity metric to compare cells within the same clusters. How similar are their CNV profiles. 

**Analysis Plan:**
> Average (or median) CNV profiles across cells in each transcriptional cluster.  
> Side hisotgram of proportion of cells from that sample in that clone.  
> Similarity metric for eah transcriptinal cluster to show similar or not cells within that group are.   
> Somehow need to integrate WGS of Tissue and Line onto each of these plots 


**Analysis Dir:** /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs

---
## 1.0 Average CNV profiles across transcriptional clusters
---

Using the mean for now, but could also take the median for each profile instead of the mean, because it is less sensitive to outliers. The obervations file and the expression_post_viz_transform.txt are the same values used for plotting. Just need to remove the normal cells from the post_viz matrix. 
  
Load in the data when InferCNV was run on all cells across samples at once. 

> There are 6351 genes included in this analysis across the genome to give CNV profiles for each cells.     
> 595 normal brain cells used as refernece for CNV inference in BTSCs   
> 69393 BTSCs  


```R
library(data.table)
```


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/")
```

---
## 1.1 Data prep for plotting
---

Run on Samwise
## load in data and subset BTSCs vs ref cells

obs <- data.table::fread("./output/expression_post_viz_transform.txt")

obs <- data.frame(obs)
str(obs)
rownames(obs) <- obs$V1
obs <- obs[ , -1]
obs[1:5, 1:5]

head(obs)
dim(obs)
max(obs)
min(obs)## Order the genes by genome position and order

genePos <- read.table("./input/GenePos_GRCh38.txt")
head(genePos)

#get all genes in CNV plot from genePos
CNV.genes <- genePos[genePos$V1 %in% rownames(obs), ]
dim(CNV.genes)
dim(obs)


CNV.genes <- CNV.genes[ ,1:2]
rownames(CNV.genes) <- CNV.genes[,1]
CNV.genes[,1] <- NULL
head(CNV.genes)
colnames(CNV.genes) <- c("Chromosome")

order <- c(grep("^1$", CNV.genes$Chromosome),
grep("^2$", CNV.genes$Chromosome),
grep("^3$", CNV.genes$Chromosome),
grep("^4$", CNV.genes$Chromosome),
grep("^5$", CNV.genes$Chromosome),
grep("^6$", CNV.genes$Chromosome),
grep("^7$", CNV.genes$Chromosome),
grep("^8$", CNV.genes$Chromosome),
grep("^9$", CNV.genes$Chromosome),
grep("^10$", CNV.genes$Chromosome),
grep("^11$", CNV.genes$Chromosome),
grep("^12$", CNV.genes$Chromosome),
grep("^13$", CNV.genes$Chromosome),
grep("^14$", CNV.genes$Chromosome),
grep("^15$", CNV.genes$Chromosome),
grep("^16$", CNV.genes$Chromosome),
grep("^17$", CNV.genes$Chromosome),
grep("^18$", CNV.genes$Chromosome),
grep("^19$", CNV.genes$Chromosome),
grep("^20$", CNV.genes$Chromosome),
grep("^21$", CNV.genes$Chromosome),
grep("^22$", CNV.genes$Chromosome)
           )

CNV.genes$Number <- 1:length(CNV.genes$Chromosome)
CNV.genes <- CNV.genes[match(order, CNV.genes$Number),]
CNV.genes$Number <- NULL
head(CNV.genes)


##order expression matrix by this
order <- rownames(CNV.genes)
obs <- obs[order, ]#subset out the BTSCs (69.3k cells)
BTSC.CNVs <- obs[,-grep("_T_", colnames(obs))] #normal cells are from the tumour
colnames(BTSC.CNVs) <- gsub("VIZ.", "", colnames(BTSC.CNVs))

head(BTSC.CNVs)
dim(BTSC.CNVs)
max(BTSC.CNVs)
min(BTSC.CNVs)

#subset out the normal oligo cells (595 cells)
ref.CNVs <- obs[,grep("_T_", colnames(obs))]
colnames(ref.CNVs) <- gsub("VIZ.", "", colnames(ref.CNVs))
#head(ref.CNVs)
dim(ref.CNVs)
max(ref.CNVs)
min(ref.CNVs)### Average CNV profiles across clusters from Intra-BTSC clustering results (should be 82ish clusters)

load("/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/opt_solution/BTSC_IntrCluster_metadata.RData")

intra.meta$UniqueID <- paste(intra.meta$orig.ident, intra.meta$Cluster.ID, sep = "_")
intra.meta$UniqueID[1:10]

colnames(BTSC.CNVs) <- gsub("G946.", "G946-", colnames(BTSC.CNVs))
colnames(BTSC.CNVs) <- gsub("G945.", "G945-", colnames(BTSC.CNVs))

clusters <- unique(intra.meta$UniqueID)
avg.cnv.df <- data.frame(A=NA)

for (i in 1:length(clusters)){
    
    print("-------")
    print(i)
    clust <- clusters[i]
    print(clust)
    cells <- rownames(intra.meta[grep(clust, intra.meta$UniqueID), ])
    print(length(cells))
    
    subset <- BTSC.CNVs[, cells]
    dim(subset)
    subset[1:3, 1:3]
    
    avg <- rowMeans(subset)
    ## load CNV data across cohort and 
    avg <- data.frame(avg)
    avg.cnv.df <- cbind(avg.cnv.df, avg)
    colnames(avg.cnv.df)[i+1] <- clust
     
}

avg.cnv.df <- avg.cnv.df[, -1]
dim(avg.cnv.df) #86 clusters by 6351 genes### save all the important data into one loadable R object

save(BTSC.CNVs, 
     ref.CNVs, 
     CNV.genes, 
     intra.meta,
     avg.cnv.df,
     file = "GlobalBTSC_CNVs.RData")
     
save(BTSC.CNVs, 
     CNV.genes, 
     intra.meta,
     avg.cnv.df,
     file = "GlobalBTSC_inferCNV.RData")

save(
     CNV.genes, 
     intra.meta,
     avg.cnv.df,
     file = "GlobalBTSC_CNVs_averge.RData")

save(ref.CNVs,
     file = "GlobalBTSC_CNVs_referenceCells.RData")
----
## 1.2 Plot average transcriptional clusters
---


```R
library(pheatmap)
library(RColorBrewer)
library(grid)
```

    Warning message:
    “package ‘pheatmap’ was built under R version 3.4.4”


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/")
load("GlobalBTSC_CNVs_averge.RData")
ls()
```


<ol class=list-inline>
	<li>'avg.cnv.df'</li>
	<li>'CNV.genes'</li>
	<li>'intra.meta'</li>
</ol>




```R
##annotate columns by alternating CNV colors
##annotate rows by sample ID -- use the colors from the mega tSNE.

avg.cnv.df <- avg.cnv.df[,order(colnames(avg.cnv.df))]
```


```R
head(avg.cnv.df)
```


<table>
<thead><tr><th></th><th scope=col>BT127_L_C1</th><th scope=col>BT127_L_C2</th><th scope=col>BT147_L_C1</th><th scope=col>BT147_L_C2</th><th scope=col>BT48_L_C1</th><th scope=col>BT48_L_C2</th><th scope=col>BT67_L_C1</th><th scope=col>BT67_L_C2</th><th scope=col>BT73_L_C1</th><th scope=col>BT73_L_C2</th><th scope=col>⋯</th><th scope=col>G945-J_L_C2</th><th scope=col>G945-K_L_C1</th><th scope=col>G945-K_L_C2</th><th scope=col>G945-K_L_C3</th><th scope=col>G946-J_L_C1</th><th scope=col>G946-J_L_C2</th><th scope=col>G946-J_L_C3</th><th scope=col>G946-J_L_C4</th><th scope=col>G946-K_L_C1</th><th scope=col>G946-K_L_C2</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>0.02547761   </td><td>0.020312710  </td><td>0.03504792   </td><td>0.01399746   </td><td>-0.014402775 </td><td>-0.03853899  </td><td>-0.10570094  </td><td>-0.09201973  </td><td>-0.03080463  </td><td> 0.0008629272</td><td>⋯            </td><td>-0.04217228  </td><td>0.1914731    </td><td>0.1921001    </td><td>0.07589099   </td><td>0.1030175    </td><td>-0.11005032  </td><td>0.03903750   </td><td>0.034182131  </td><td>-0.1287035   </td><td>-0.1483973   </td></tr>
	<tr><th scope=row>TMEM201</th><td>0.03122640   </td><td>0.022181641  </td><td>0.03990202   </td><td>0.01693410   </td><td>-0.009331485 </td><td>-0.03430632  </td><td>-0.09366095  </td><td>-0.08518207  </td><td>-0.02512797  </td><td> 0.0048454969</td><td>⋯            </td><td>-0.03669081  </td><td>0.2003386    </td><td>0.2001252    </td><td>0.09003453   </td><td>0.1060420    </td><td>-0.11106243  </td><td>0.04435347   </td><td>0.026516149  </td><td>-0.1311349   </td><td>-0.1447148   </td></tr>
	<tr><th scope=row>CLSTN1</th><td>0.03387397   </td><td>0.026363608  </td><td>0.03787451   </td><td>0.01187433   </td><td>-0.007656861 </td><td>-0.04092358  </td><td>-0.09216101  </td><td>-0.08490364  </td><td>-0.02563131  </td><td> 0.0081523679</td><td>⋯            </td><td>-0.03342964  </td><td>0.1990865    </td><td>0.2004908    </td><td>0.09540203   </td><td>0.1075512    </td><td>-0.11340732  </td><td>0.04536863   </td><td>0.028701350  </td><td>-0.1343353   </td><td>-0.1462613   </td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>0.03354235   </td><td>0.018246244  </td><td>0.03495229   </td><td>0.00805269   </td><td>-0.009156022 </td><td>-0.03593608  </td><td>-0.08319649  </td><td>-0.07528557  </td><td>-0.03290911  </td><td>-0.0095088953</td><td>⋯            </td><td>-0.03600168  </td><td>0.1962531    </td><td>0.1990476    </td><td>0.08756835   </td><td>0.1003282    </td><td>-0.11374022  </td><td>0.04263548   </td><td>0.016935268  </td><td>-0.1368998   </td><td>-0.1466092   </td></tr>
	<tr><th scope=row>LZIC</th><td>0.02839939   </td><td>0.007586204  </td><td>0.04058198   </td><td>0.01523971   </td><td>-0.006104221 </td><td>-0.03232243  </td><td>-0.07266495  </td><td>-0.06440919  </td><td>-0.03097119  </td><td>-0.0223801521</td><td>⋯            </td><td>-0.02649233  </td><td>0.2166301    </td><td>0.2183155    </td><td>0.10053392   </td><td>0.1239792    </td><td>-0.09231483  </td><td>0.06145141   </td><td>0.007430882  </td><td>-0.1133540   </td><td>-0.1190407   </td></tr>
	<tr><th scope=row>NMNAT1</th><td>0.03414895   </td><td>0.014567140  </td><td>0.04250626   </td><td>0.01571912   </td><td> 0.004076344 </td><td>-0.02074168  </td><td>-0.05806666  </td><td>-0.05709092  </td><td>-0.02500879  </td><td>-0.0133474736</td><td>⋯            </td><td>-0.02695886  </td><td>0.2123001    </td><td>0.2109533    </td><td>0.09369069   </td><td>0.1276709    </td><td>-0.08414843  </td><td>0.06627460   </td><td>0.014892392  </td><td>-0.1009751   </td><td>-0.1018223   </td></tr>
</tbody>
</table>




```R
annotat.row <- data.frame(colnames(avg.cnv.df))
annotat.row$Sample <- gsub('.{3}$', '', colnames(avg.cnv.df))
rownames(annotat.row) <- annotat.row$colnames.avg.cnv.df.
annotat.row$colnames.avg.cnv.df. <- NULL
head(annotat.row)
unique(annotat.row$Sample)
```


<table>
<thead><tr><th></th><th scope=col>Sample</th></tr></thead>
<tbody>
	<tr><th scope=row>BT127_L_C1</th><td>BT127_L</td></tr>
	<tr><th scope=row>BT127_L_C2</th><td>BT127_L</td></tr>
	<tr><th scope=row>BT147_L_C1</th><td>BT147_L</td></tr>
	<tr><th scope=row>BT147_L_C2</th><td>BT147_L</td></tr>
	<tr><th scope=row>BT48_L_C1</th><td>BT48_L </td></tr>
	<tr><th scope=row>BT48_L_C2</th><td>BT48_L </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'BT127_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
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
	<li>'G797_L'</li>
	<li>'G799_L'</li>
	<li>'G800_L'</li>
	<li>'G837_L'</li>
	<li>'G851_L'</li>
	<li>'G876_L'</li>
	<li>'G885_L'</li>
	<li>'G895_L'</li>
	<li>'G945-I_L'</li>
	<li>'G945-J_L'</li>
	<li>'G945-K_L'</li>
	<li>'G946-J_L'</li>
	<li>'G946-K_L'</li>
</ol>




```R
colfunc <- colorRampPalette(c("#54278f", "#bcbddc", "#084081", "#4eb3d3", "#238b45", "#ccebc5"))
dirks <- colfunc(21)
#plot(rep(1,21),col=colfunc(21),pch=19,cex=3)

colfunc <- colorRampPalette(c("#800026", "#fc4e2a", "#feb24c", "#ffeda0"))
weiss <- colfunc(8)
#plot(rep(1,8),col=colfunc(8),pch=19,cex=3)

cols <- c(weiss, dirks)
length(cols)
```


29



```R

    mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11)),
                       Sample = cols
                   )

    names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
    names(mat_colors$Sample) <- unique(annotat.row$Sample)
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
colfunc <- colorRampPalette(cols)
cnv.cols <- colfunc(50)
```


<dl>
	<dt>$Chromosome</dt>
		<dd><dl class=dl-horizontal>
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
</dd>
	<dt>$Sample</dt>
		<dd><dl class=dl-horizontal>
	<dt>BT127_L</dt>
		<dd>'#800026'</dd>
	<dt>BT147_L</dt>
		<dd>'#B52127'</dd>
	<dt>BT48_L</dt>
		<dd>'#EA4229'</dd>
	<dt>BT67_L</dt>
		<dd>'#FC6A33'</dd>
	<dt>BT73_L</dt>
		<dd>'#FD9542'</dd>
	<dt>BT84_L</dt>
		<dd>'#FEBA57'</dd>
	<dt>BT89_L</dt>
		<dd>'#FED37B'</dd>
	<dt>BT94_L</dt>
		<dd>'#FFEDA0'</dd>
	<dt>G523_L</dt>
		<dd>'#54278F'</dd>
	<dt>G549_L</dt>
		<dd>'#6E4CA2'</dd>
	<dt>G564_L</dt>
		<dd>'#8872B5'</dd>
	<dt>G566_L</dt>
		<dd>'#A297C8'</dd>
	<dt>G583_L</dt>
		<dd>'#BCBDDC'</dd>
	<dt>G620_L</dt>
		<dd>'#8F9DC5'</dd>
	<dt>G637_L</dt>
		<dd>'#617EAE'</dd>
	<dt>G729_L</dt>
		<dd>'#345F97'</dd>
	<dt>G797_L</dt>
		<dd>'#084081'</dd>
	<dt>G799_L</dt>
		<dd>'#195C95'</dd>
	<dt>G800_L</dt>
		<dd>'#2A79A9'</dd>
	<dt>G837_L</dt>
		<dd>'#3C96BE'</dd>
	<dt>G851_L</dt>
		<dd>'#4EB3D3'</dd>
	<dt>G876_L</dt>
		<dd>'#43A9AF'</dd>
	<dt>G885_L</dt>
		<dd>'#389F8B'</dd>
	<dt>G895_L</dt>
		<dd>'#2D9568'</dd>
	<dt>G945-I_L</dt>
		<dd>'#238B45'</dd>
	<dt>G945-J_L</dt>
		<dd>'#4DA365'</dd>
	<dt>G945-K_L</dt>
		<dd>'#77BB84'</dd>
	<dt>G946-J_L</dt>
		<dd>'#A1D3A5'</dd>
	<dt>G946-K_L</dt>
		<dd>'#CCEBC5'</dd>
</dl>
</dd>
</dl>




```R
#sample.order <- as.character(read.table("~/Downloads/increasingproneural.txt", sep = "\t")[,1])
#sample.order <- sample.order[!sample.order == ""]
#sample.order

#Redorder samples by cluster number

sample.order <- c('BT127_L', 
                  'BT147_L',
                  'BT48_L',
                  'BT67_L',
                  'BT73_L',
                  'BT89_L',
                  'BT94_L',
                  'G549_L',
                  'G564_L',
                  'G729_L',
                  'G851_L',
                  'G945-J_L',
                  'G946-K_L',
                  'BT84_L',
                  'G523_L',
                  'G566_L',
                  
                  'G799_L',
                  'G800_L',
                  'G876_L',
                  'G885_L',
                  'G895_L',
                  'G945-K_L',
                  
                  'G583_L',
                  'G945-I_L',
                  'G946-J_L',
                  'G637_L',
                  'G797_L',
                  'G837_L',
                  'G620_L'
                  )

test <- data.frame(colnames(avg.cnv.df))
colnames(test) <- "ClusterID"

test$SampleID <- gsub('.{3}$', '', test$ClusterID)

test$SampleID <- reorder.factor(test$SampleID, new.order=sample.order)

bb <- test %>%
  arrange(SampleID)

row.order <- as.character(bb$ClusterID)
row.order

avg.cnv.df <- avg.cnv.df[, row.order]
head(avg.cnv.df) ##reorder dataframe
```


<ol class=list-inline>
	<li>'BT127_L_C1'</li>
	<li>'BT127_L_C2'</li>
	<li>'BT147_L_C1'</li>
	<li>'BT147_L_C2'</li>
	<li>'BT48_L_C1'</li>
	<li>'BT48_L_C2'</li>
	<li>'BT67_L_C1'</li>
	<li>'BT67_L_C2'</li>
	<li>'BT73_L_C1'</li>
	<li>'BT73_L_C2'</li>
	<li>'BT89_L_C1'</li>
	<li>'BT89_L_C2'</li>
	<li>'BT94_L_C1'</li>
	<li>'BT94_L_C2'</li>
	<li>'G549_L_C1'</li>
	<li>'G549_L_C2'</li>
	<li>'G564_L_C1'</li>
	<li>'G564_L_C2'</li>
	<li>'G729_L_C1'</li>
	<li>'G729_L_C2'</li>
	<li>'G851_L_C1'</li>
	<li>'G851_L_C2'</li>
	<li>'G945-J_L_C1'</li>
	<li>'G945-J_L_C2'</li>
	<li>'G946-K_L_C1'</li>
	<li>'G946-K_L_C2'</li>
	<li>'BT84_L_C1'</li>
	<li>'BT84_L_C2'</li>
	<li>'BT84_L_C3'</li>
	<li>'G523_L_C1'</li>
	<li>'G523_L_C2'</li>
	<li>'G523_L_C3'</li>
	<li>'G566_L_C1'</li>
	<li>'G566_L_C2'</li>
	<li>'G566_L_C3'</li>
	<li>'G799_L_C1'</li>
	<li>'G799_L_C2'</li>
	<li>'G799_L_C3'</li>
	<li>'G800_L_C1'</li>
	<li>'G800_L_C2'</li>
	<li>'G800_L_C3'</li>
	<li>'G876_L_C1'</li>
	<li>'G876_L_C2'</li>
	<li>'G876_L_C3'</li>
	<li>'G885_L_C1'</li>
	<li>'G885_L_C2'</li>
	<li>'G885_L_C3'</li>
	<li>'G895_L_C1'</li>
	<li>'G895_L_C2'</li>
	<li>'G895_L_C3'</li>
	<li>'G945-K_L_C1'</li>
	<li>'G945-K_L_C2'</li>
	<li>'G945-K_L_C3'</li>
	<li>'G583_L_C1'</li>
	<li>'G583_L_C2'</li>
	<li>'G583_L_C3'</li>
	<li>'G583_L_C4'</li>
	<li>'G945-I_L_C1'</li>
	<li>'G945-I_L_C2'</li>
	<li>'G945-I_L_C3'</li>
	<li>'G945-I_L_C4'</li>
	<li>'G946-J_L_C1'</li>
	<li>'G946-J_L_C2'</li>
	<li>'G946-J_L_C3'</li>
	<li>'G946-J_L_C4'</li>
	<li>'G637_L_C1'</li>
	<li>'G637_L_C2'</li>
	<li>'G637_L_C3'</li>
	<li>'G637_L_C4'</li>
	<li>'G637_L_C5'</li>
	<li>'G797_L_C1'</li>
	<li>'G797_L_C2'</li>
	<li>'G797_L_C3'</li>
	<li>'G797_L_C4'</li>
	<li>'G797_L_C5'</li>
	<li>'G837_L_C1'</li>
	<li>'G837_L_C2'</li>
	<li>'G837_L_C3'</li>
	<li>'G837_L_C4'</li>
	<li>'G837_L_C5'</li>
	<li>'G620_L_C1'</li>
	<li>'G620_L_C2'</li>
	<li>'G620_L_C3'</li>
	<li>'G620_L_C4'</li>
	<li>'G620_L_C5'</li>
	<li>'G620_L_C6'</li>
</ol>




```R
file.name <- "~/Desktop/ClusterCNVs_noLEGEND_newCohort_clusterOrder.jpeg"

pheatmap(t(avg.cnv.df),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = CNV.genes,
         annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = TRUE,
         height = 15,
         width = 25,
         filename = file.name
        )
```

----
## 1.3 Plot CNV heatmap for reference brain cells
----


```R
#load in single cell data for normal brain cells used as the reference control

load("GlobalBTSC_CNVs_referenceCells.RData")

ls()

```


<ol class=list-inline>
	<li>'annotat.row'</li>
	<li>'avg.cnv.df'</li>
	<li>'cnv.cols'</li>
	<li>'CNV.genes'</li>
	<li>'colfunc'</li>
	<li>'cols'</li>
	<li>'dirks'</li>
	<li>'file.name'</li>
	<li>'intra.meta'</li>
	<li>'mat_colors'</li>
	<li>'ref.CNVs'</li>
	<li>'weiss'</li>
</ol>




```R
#order by sample

ref.CNVs <- ref.CNVs[,order(colnames(ref.CNVs))]
head(ref.CNVs)

samples
```


<table>
<thead><tr><th></th><th scope=col>G1003.A_T_AGGCCACCATTTGCTT</th><th scope=col>G1003.A_T_CAGAGAGCAAGAAAGG</th><th scope=col>G1003.A_T_CCACCTACAGTTTACG</th><th scope=col>G1003.A_T_CCACGGAGTGTTTGTG</th><th scope=col>G1003.A_T_CTACACCGTTGTGGCC</th><th scope=col>G1003.A_T_CTGCCTACAGAAGCAC</th><th scope=col>G1003.A_T_GCAATCATCAACCATG</th><th scope=col>G1003.A_T_GCAGTTACAGATGGGT</th><th scope=col>G1003.A_T_GCTGGGTCAAGCGAGT</th><th scope=col>G1003.A_T_GGCCGATCACACGCTG</th><th scope=col>⋯</th><th scope=col>G983.C_T_GTCCTCAAGCCAACAG</th><th scope=col>G983.C_T_GTCGTAAAGTCGATAA</th><th scope=col>G983.C_T_GTTAAGCCACGCCAGT</th><th scope=col>G983.C_T_TACTTGTTCGCGCCAA</th><th scope=col>G983.C_T_TCGGTAAGTCGGCTCA</th><th scope=col>G983.C_T_TGAGAGGCATGAACCT</th><th scope=col>G983.C_T_TGCTACCGTCGGATCC</th><th scope=col>G983.C_T_TGGACGCCAAACAACA</th><th scope=col>G983.C_T_TGTGGTAAGCACACAG</th><th scope=col>G983.C_T_TTATGCTCAATCAGAA</th></tr></thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>0         </td><td>-0.2087647</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3162926</td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
	<tr><th scope=row>TMEM201</th><td>0         </td><td>-0.2098387</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3173665</td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
	<tr><th scope=row>CLSTN1</th><td>0         </td><td>-0.2053132</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3128411</td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>0         </td><td>-0.2169603</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3244882</td><td>0         </td><td>-0.1010671</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
	<tr><th scope=row>LZIC</th><td>0         </td><td>-0.2310089</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td> 0.0000000</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3385368</td><td>0         </td><td>-0.1151157</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
	<tr><th scope=row>NMNAT1</th><td>0         </td><td>-0.2204849</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>-0.1033578</td><td>0         </td><td>0         </td><td>0         </td><td>⋯         </td><td>0         </td><td>-0.3280127</td><td>0         </td><td>-0.1045916</td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td><td>0         </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.A_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.B_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.C_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G1003.D_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.A_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.B_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.C_T'</li>
	<li>'G910.D_T'</li>
	<li>'G910.D_T'</li>
	<li>'G910.D_T'</li>
	<li>'G910.D_T'</li>
	<li>'G910.D_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G910.E_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.I_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.J_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G945.K_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.I_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.J_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G946.K_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.A_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.B_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.C_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G967.D_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.A_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.B_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
	<li>'G983.C_T'</li>
</ol>




```R
annotat.row <- data.frame(colnames(ref.CNVs))
annotat.row$Sample <- gsub('.{17}$', '', colnames(ref.CNVs))
rownames(annotat.row) <- annotat.row$colnames.ref.CNVs.
annotat.row$colnames.ref.CNVs.<- NULL
head(annotat.row)

```


<table>
<thead><tr><th></th><th scope=col>Sample</th></tr></thead>
<tbody>
	<tr><th scope=row>G1003.A_T_AGGCCACCATTTGCTT</th><td>G1003.A_T</td></tr>
	<tr><th scope=row>G1003.A_T_CAGAGAGCAAGAAAGG</th><td>G1003.A_T</td></tr>
	<tr><th scope=row>G1003.A_T_CCACCTACAGTTTACG</th><td>G1003.A_T</td></tr>
	<tr><th scope=row>G1003.A_T_CCACGGAGTGTTTGTG</th><td>G1003.A_T</td></tr>
	<tr><th scope=row>G1003.A_T_CTACACCGTTGTGGCC</th><td>G1003.A_T</td></tr>
	<tr><th scope=row>G1003.A_T_CTGCCTACAGAAGCAC</th><td>G1003.A_T</td></tr>
</tbody>
</table>




```R
mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11)),
                   Sample = c(rep(c("black", "grey"), 11))
                   )
names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
names(mat_colors$Sample) <- unique(annotat.row$Sample)
mat_colors
```


<dl>
	<dt>$Chromosome</dt>
		<dd><dl class=dl-horizontal>
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
</dd>
	<dt>$Sample</dt>
		<dd><dl class=dl-horizontal>
	<dt>G1003.A_T</dt>
		<dd>'black'</dd>
	<dt>G1003.B_T</dt>
		<dd>'grey'</dd>
	<dt>G1003.C_T</dt>
		<dd>'black'</dd>
	<dt>G1003.D_T</dt>
		<dd>'grey'</dd>
	<dt>G910.A_T</dt>
		<dd>'black'</dd>
	<dt>G910.B_T</dt>
		<dd>'grey'</dd>
	<dt>G910.C_T</dt>
		<dd>'black'</dd>
	<dt>G910.D_T</dt>
		<dd>'grey'</dd>
	<dt>G910.E_T</dt>
		<dd>'black'</dd>
	<dt>G945.I_T</dt>
		<dd>'grey'</dd>
	<dt>G945.J_T</dt>
		<dd>'black'</dd>
	<dt>G945.K_T</dt>
		<dd>'grey'</dd>
	<dt>G946.I_T</dt>
		<dd>'black'</dd>
	<dt>G946.J_T</dt>
		<dd>'grey'</dd>
	<dt>G946.K_T</dt>
		<dd>'black'</dd>
	<dt>G967.A_T</dt>
		<dd>'grey'</dd>
	<dt>G967.B_T</dt>
		<dd>'black'</dd>
	<dt>G967.C_T</dt>
		<dd>'grey'</dd>
	<dt>G967.D_T</dt>
		<dd>'black'</dd>
	<dt>G983.A_T</dt>
		<dd>'grey'</dd>
	<dt>G983.B_T</dt>
		<dd>'black'</dd>
	<dt>G983.C_T</dt>
		<dd>'grey'</dd>
</dl>
</dd>
</dl>




```R
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
colfunc <- colorRampPalette(cols)
cnv.cols <- colfunc(50)
```


```R
file.name <- "~/Desktop/ClusterCNVs__BTSC_reference_nolegened.jpeg"

pheatmap(t(ref.CNVs),
         color = rev(cnv.cols),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = CNV.genes,
         annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = FALSE,
         height = 15,
         width = 25,
         filename = file.name
        )
```

----
## 1.4 Output CNV heatmap at single cell level for each sample
---

Hierachal clustering of cells to identify common CNV clones within the samples

## on samwise load in the big CNV matrix and divide up by sample ID
## BTSC.CNVs

samples <- unique(gsub('.{17}$', '', colnames(BTSC.CNVs)))
samples #29

for (i in 1:length(samples)){
    
    print("")
    print("----------")
    print(samples[i])

    sample.CNV <- BTSC.CNVs[ ,grep(samples[i], colnames(BTSC.CNVs))]
    print(dim(sample.CNV))
    print(colnames(sample.CNV[1:5]))
    
    sample.meta <- intra.meta[grep(samples[i], rownames(intra.meta)), ]
    print(dim(sample.meta))
    print(rownames(sample.meta)[1:5])
    
    save.file.name <- paste0(samples[i], "_cellCNVs.RData")
    print(save.file.name)
    save(sample.CNV, sample.meta, CNV.genes, file = save.file.name)
    
}

----
  
module load R  
Rscript plotSampleCNVs.R  > plotSampleCNVs.log
library(pheatmap)
library(RColorBrewer)
library(grid)


test_match_order <- function(x,y) {

if (all(x==y)) print('Perfect match in same order')

if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')

if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}

#test_match_order(x,y)


files <- list.files("./", pattern = "_cellCNVs.RData")
samples <- gsub("_cellCNVs.RData", "", files)
samples
    
    for (i in 1:length(samples)){
    
    print("")
    print("----------")
    print(samples[i])

    load.file <- paste0(samples[i], "_cellCNVs.RData")
    print(load.file)
    load(load.file)

    print("Check cell BC vectors are in the same order:")
    test_match_order(colnames(sample.CNV), rownames(sample.meta))
    
    #annotate with cluster IDs
    print("Annotating columns and clusters:")

    annotat.row <- data.frame(colnames(sample.CNV))
    annotat.row$Cluster <- sample.meta$UniqueID
    rownames(annotat.row) <- annotat.row$colnames.sample.CNV.
    annotat.row$colnames.sample.CNV.<- NULL
    head(annotat.row)

    #get colors for plot

    num <- length((unique(annotat.row$Cluster)))

    mat_colors <- list(Chromosome = c(rep(c("lightgrey", "darkgrey"), 11)),
                   Cluster = c(brewer.pal(num, "Dark2")[1:num])
                   )
    names(mat_colors$Chromosome) <- unique(CNV.genes$Chromosome)
    names(mat_colors$Cluster) <- sort(unique(annotat.row$Cluster))
    print(mat_colors)

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
    colfunc <- colorRampPalette(cols)
    cnv.cols <- colfunc(50)


    file.name <- paste0(samples[i], "_cellCNV.jpeg")
    print(file.name)

    pheatmap(t(sample.CNV),
         color = rev(cnv.cols),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_col = CNV.genes,
         annotation_row = annotat.row,
         annotation_colors = mat_colors,
         fontsize_row = 9,
         border_color = "black",
         annotation_legend = FALSE,
         legend = FALSE,
         height = 15,
         width = 27,
         filename = file.name
        )
    
    
    }
----
## 1.4 Correlate the CNV profiles across cells within clusters
---

Do cells in the same transcriptional cluster have similar CNVs? More similar compared to cells in other clusters?

This can be a histogram. With each bar is a cluster. Bar goes up to average, and error bar is SE. Similar to figure 4A (bottom) in Venteicher et al., Science, 2017
files <- list.files("./", pattern = "_cellCNVs.RData")
samples <- gsub("_cellCNVs.RData", "", files)
samples

for (i in 1:length(samples)){


    print("")
    print("-------")
    print(samples[i])


    print("Load in CNV data")
    load.file <- paste0(samples[i], "_cellCNVs.RData")
    print(load.file)
    load(load.file)

    print("Subset CNV matrix by cluster ID")

    clusters <- names(table(sample.meta$UniqueID))
    print(clusters)

    cluster.cnv <- list()
    pearson.cor.cnv <- list()
    cor.means <- c()
    names <- c()
    
    for (j in 1:length(clusters)){
        
        print("")
        print("-------")
        print(clusters[j])
        cells <- rownames(sample.meta[sample.meta$UniqueID == clusters[j], ])
        print(length(cells))  
        
        #subet CN
        dat <- data.frame(sample.CNV[ ,cells])
        print(dim(dat))
        
        cluster.cnv[[j]] <- dat
        
        print("")
        print("Calculate pearson correlation between cells within each cluster")
        cor <- cor(cluster.cnv[[j]])
        pearson.cor.cnv[[j]] <- cor
        
        print("")
        print("Mean cluster correlation")
        
        avg <- mean(cor)
        print(avg)
        cor.means <- c(cor.means, avg)
        names <- c(names, clusters[j])

    }

    names(cluster.cnv) <- clusters
    names(pearson.cor.cnv) <- clusters

    print("")
    print("-----")
    print("Calculate pearson correlation between ALL cells in sample") 
    
    cor <- cor(sample.CNV)
    pearson.cor.cnv[[length(clusters)+1]] <- cor
    names(pearson.cor.cnv) <- c(clusters, samples[i])

    print("")
    print("Mean correlation between all cells")
    avg <- mean(cor)
    print(avg)
    cor.means <- c(cor.means, avg)
    names <- c(names, samples[i])


    print("")
    print("Save Data for sample:")

    file.name <- paste0(samples[i], "_PearsonCorrelation_CNV.Rdata")
    print(file.name)

    save(pearson.cor.cnv,
            cluster.cnv,
         file = file.name
        
        )

}

names(cor.means) <- names
print(cor.means)

save(cor.means, file = "Cluster_Sample_MeanCorrelations.Rdata")

print("END")

```R
cor.means
```


<dl class=dl-horizontal>
	<dt>BT127_L_C1</dt>
		<dd>0.422828000889101</dd>
	<dt>BT127_L_C2</dt>
		<dd>0.348598356936055</dd>
	<dt>BT127_L</dt>
		<dd>0.392988816769633</dd>
</dl>



---
## 1.5 Plot line graph of clone CNV profiles
---

### G637_L


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/")
load("GlobalBTSC_CNVs_averge.RData")
ls()

```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'avg.cnv.df'</li><li>'CNV.genes'</li><li>'intra.meta'</li></ol>




```R
#load in averaged CNV clones from global dataset
#setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/inferCNV/BTSCs/")
#load("GlobalBTSC_CNVs_averge.RData")
#ls()

chr <- table(CNV.genes$Chromosome)

breaks <- c()
for (i in 1:length(chr)){
    breaks[i] <- sum(chr[1:i])
}

```


```R
id <- "BT67_L"
```


```R
sample <- avg.cnv.df[ ,grep(id, colnames(avg.cnv.df))]
sample$Index <- seq(from = 1, to = nrow(sample), by = 1)
head(sample)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>BT67_L_C2</th><th scope=col>BT67_L_C1</th><th scope=col>Index</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>SLC25A33</th><td>-0.09201973</td><td>-0.10570094</td><td>1</td></tr>
	<tr><th scope=row>TMEM201</th><td>-0.08518207</td><td>-0.09366095</td><td>2</td></tr>
	<tr><th scope=row>CLSTN1</th><td>-0.08490364</td><td>-0.09216101</td><td>3</td></tr>
	<tr><th scope=row>CTNNBIP1</th><td>-0.07528557</td><td>-0.08319649</td><td>4</td></tr>
	<tr><th scope=row>LZIC</th><td>-0.06440919</td><td>-0.07266495</td><td>5</td></tr>
	<tr><th scope=row>NMNAT1</th><td>-0.05709092</td><td>-0.05806666</td><td>6</td></tr>
</tbody>
</table>




```R
#larger smoothing span value means more smooth of the fitted curve

loessMod25.1 <- loess(sample[,1] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.1 <- predict(loessMod25.1)


loessMod25.2 <- loess(sample[,2] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.2 <- predict(loessMod25.2)

loessMod25.3 <- loess(sample[,3] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.3 <- predict(loessMod25.3)


```


```R
pdf("~/Desktop/BT67_L_CNVclones_Nov2020.pdf", height = 4, width = 8)
#pdf("~/Desktop/G876_L_CNVclones.pdf", height = 4, width = 8)

plot(sample[,1],
     ylim = c(-0.5, 0.7),
     cex = 0.3,
     col = "white",
     pch = 19,
     ylab = "CNV Score"
    )
lines(smoothed25.1, lwd = 2, col="#1b9e77")
lines(smoothed25.2, lwd = 2, col= "#d95f02")
lines(smoothed25.3, lwd = 2, col= "#7570B3")
abline(v=c(0,breaks), lty =1, col = "black")
abline(h=c(-0.05, 0.05), lty = 2, col = "darkgrey")

par(mar=rep(0,4))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center",
       legend = c("G876_L.C1",
                  "G876_L.C2",
                  "G876_L.C3"
                 ),
       lty = 1,
       lwd = 5,
       col = c("#1b9e77", "#d95f02", "#7570B3"),
       horiz = TRUE,
       bty = 'n'
      
      )


dev.off()
```


<strong>pdf:</strong> 2

loessMod25.3 <- loess(sample[,3] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.3 <- predict(loessMod25.3)


loessMod25.4 <- loess(sample[,4] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.4 <- predict(loessMod25.4)


loessMod25.5 <- loess(sample[,5] ~ Index, 
                    dat = sample, 
                    span=0.05
                   ) # 25% smoothing span
smoothed25.5 <- predict(loessMod25.5)pdf("~/Desktop/G797_L_CNVclones.pdf", height = 4, width = 8)

plot(sample[,1],
     ylim = c(-0.5, 0.7),
     cex = 0.3,
     col = "white",
     pch = 19,
     ylab = "CNV Score"
    )
lines(smoothed25.1, lwd = 2, col="#1b9e77")
lines(smoothed25.2, lwd = 2, col= "#d95f02")
lines(smoothed25.3, lwd = 2, col= "#d95f02")
lines(smoothed25.4, lwd = 2, col= "#d95f02")
lines(smoothed25.5, lwd = 2, col= "#d95f02")
abline(v=c(0,breaks), lty =1, col = "black")
abline(h=c(-0.05, 0.05), lty = 2, col = "darkgrey")

par(mar=rep(0,4))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center",
       legend = c("BT67_L.C1",
                  "BT67_L.C2"
                 ),
       lty = 1,
       lwd = 5,
       col = c("#1b9e77", "#d95f02"),
       horiz = TRUE,
       bty = 'n'
      
      )


dev.off()

```R
library(RColorBrewer)
brewer.pal(n = 5, name = "Dark2")
```


<ol class=list-inline>
	<li>'#1B9E77'</li>
	<li>'#D95F02'</li>
	<li>'#7570B3'</li>
	<li>'#E7298A'</li>
	<li>'#66A61E'</li>
</ol>




```R

```

---
#  Calculate correlation for CNV profiles
---




```R

```
