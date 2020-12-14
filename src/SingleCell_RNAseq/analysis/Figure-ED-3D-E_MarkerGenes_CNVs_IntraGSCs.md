
---
# Plot effect of CNVs on marker genes
---

/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MarkerGenes_CNVs

Reference notebook: Visualize_CNVs_Cluster_June2019


```R
setwd("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MarkerGenes_CNVs")
```

----
## 1.0 Data cleaning
----

#### 1.1 Format marker gene file
library(Seurat)

### load marker genes
#markergenes <- read.csv("TableS3_GSC_cluster_markers.csv") ### this was only top 50 fml....
markergenes <- readRDS("AllMarkers_IntraGSC.rds")
head(markergenes)
genes <- markergenes$Gene

### add cluster average expression column
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
head(BTSC@meta.data)
BTSC <- UpdateSeuratObject(BTSC)
Idents(BTSC) <- "IntraBTSC.ID"
cluster.averages <- AverageExpression(object = BTSC)
cluster.averages <- cluster.averages$RNA
saveRDS(cluster.averages, file = "IntraGSC_ClusterAvg.rds")

avg.exp <- c()
for (i in 1:nrow(markergenes)){
    avg.exp[i] <- cluster.averages[as.character(markergenes$Gene[i]),
                                   as.character(markergenes$Cluster[i])
                                  ]
}
markergenes$AverageExpression <- avg.exp
head(markergenes)

### add cluster avergae inferCNV score column
load("GlobalBTSC_CNVs_averge.RData")
head(avg.cnv.df)

avg.cnv <- c()
for (i in 1:nrow(markergenes)){
    avg.cnv[i] <- avg.cnv.df[as.character(markergenes$Gene[i]),
                                   as.character(markergenes$Cluster[i])
                                  ]
}
markergenes$AverageCNV <- avg.cnv
head(markergenes)

### add chr arm position to marker genes
gene.pos <- readRDS("GenePosition_Chrarm.rds")
chr.arm <- c()
for (i in 1:nrow(markergenes)){
    chr.arm[i] <- gene.pos[as.character(markergenes$Gene[i]), ]$arm
}
markergenes$ChromosomeArm <- chr.arm
head(markergenes)

saveRDS(markergenes, file = "Markers_Exp_CNVs.rds")
#### 1.2 Make gene position file with chromosome arm
### make a matrix of average chr arm expression across clusters
chrarms <- read.table("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/CNV_TCGA/data/chr_arms.txt", header = T)
head(chrarms)

### load start and end positions for all genes
gene.pos <- read.table("~/pughlab/projects/OICR_Brain_NucSeq/GBM/analysis/inferCNV/input/GRCh38-1.2.0_premrna_genomicPositions.txt")
colnames(gene.pos) <- c("Gene", "Chromosome", "Start", "End")
rownames(gene.pos) <- gene.pos$Gene
head(gene.pos)

### match genes to chromsome arm
gene.arm <- c()

for (i in 1:nrow(gene.pos)){

     arm.subset <- chrarms[chrarms$Chrom == gene.pos$Chromosome[i] ,]
     gene.arm[i] <- ifelse(gene.pos$End[i] > arm.subset[1 ,"End"],
           as.character(arm.subset$Idf[2]), #greater than p arm end = q
           as.character(arm.subset$Idf[1]) #NOT greater than p arm end = p
         )
}
gene.pos$arm <- gene.arm
#head(gene.pos)

saveRDS(gene.pos, file = "GenePosition_Chrarm.rds")
#### 1.3 Matrix of clusters x chr arm

Use cutoffs from Extended Data Figure 3


```R
head(gene.pos)
head(avg(avg.cnv.df))

dat <- list()
arms <- unique(gene.pos$arm)

for (i in 1:length(arms)){
    
    print(arms[i])
    genes <- rownames(gene.pos[as.character(gene.pos$arm) == arms[i], ])
    sub <- avg.cnv.df[rownames(avg.cnv.df) %in% genes, ]
    dat[[arms[i]]] <- colMeans(sub)
}

dat <- do.call(rbind, dat)
dat2 <- dat
dat[dat >= 0.17] <- 1
dat[dat <= -0.15] <- -1
dat[dat > -0.2 & dat < 1] <- 0
dat <- dat[!rownames(dat) == "NA", ]
dat[is.na(dat)] <- 0

saveRDS(dat, file = "IntraGSCcluster_chrarm_binned.rds")
```

----
## 2.0 Re-plot binary CNV heatmap
---

Use Cutoffs from WGS benchmarking to 


```R
library(pheatmap)
```


```R
dat <- readRDS("IntraGSCcluster_chrarm_binned.rds")
```


```R
pheatmap(t(dat),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         file = "BinnedChrarms_heatmap.jpeg",
         width = 6,
         height = 11.5
         )
```

----
## 3.0 Proportion marker genes in CNVs
---


```R
AllEqual <- structure(function(
	##title<< 
	## Check if all values in a vector are the same
	##description<<
	## This function is used to check if all values in a vector are equal. It can be used for example to check if a time series contains only 0 or NA values.
	
	x
	### numeric, character vector, or time series of type ts
) {
	res <- FALSE
	x <- na.omit(as.vector(x))
	if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
	return(res)
	### The function returns TRUE if all values are equal and FALSE if it contains different values.
},ex=function(){
# check if all values are equal in the following vectors:
AllEqual(1:10)
AllEqual(rep(0, 10))
AllEqual(letters)
AllEqual(rep(NA, 10))
})
```


```R
### read in marker genes
markergenes <- readRDS("Markers_Exp_CNVs.rds")
dat <- readRDS("IntraGSCcluster_chrarm_binned.rds")
genepos <- readRDS("GenePosition_Chrarm.rds")
genepos <- genepos[!is.na(genepos$arm), ]

##remove marker genes not in CNV analysis
markergenes <- markergenes[!is.na(markergenes$ChromosomeArm), ]
#head(markergenes)
```

### Figure out which regions are variable across clusters


```R
samples <- as.character(unique(markergenes$Sample))
samples
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'BT127_L'</li><li>'BT48_L'</li><li>'BT73_L'</li><li>'BT94_L'</li><li>'G523_L'</li><li>'G549_L'</li><li>'G564_L'</li><li>'G583_L'</li><li>'G620_L'</li><li>'G637_L'</li><li>'G729_L'</li><li>'G797_L'</li><li>'G800_L'</li><li>'G837_L'</li><li>'G851_L'</li><li>'G876_L'</li><li>'G885_L'</li><li>'G895_L'</li><li>'G945-I_L'</li><li>'G945-J_L'</li><li>'G945-K_L'</li><li>'G946-J_L'</li></ol>




```R
samples <- as.character(unique(markergenes$Sample))

clusters <- c()
numMarkers <- c()
numRegions <- c()
numCNV <- c()
fishers <- c()

for (i in 1:length(samples)){
    
    sub <- dat[ ,grep(samples[i], colnames(dat))]

    ### define variable regions
    var.reg <- c()
    for (j in 1:nrow(sub)){
        var.reg[j] <- AllEqual(sub[j, ])
    }
    sub <- sub[!var.reg, ]
    var.reg <- rownames(sub)
    numRegions <- append(numRegions, length(var.reg))
    if(length(var.reg) == 0){print(samples[i])} else {print(var.reg)}
    
    ### subset marker gene matrix to sample
    marker.sub <- markergenes[markergenes$Sample == samples[i], ]
    #marker.sub

    clust <- as.character(unique(marker.sub$Cluster))
    ### how many in variable regions with marker genes per cluster
    for (x in 1:length(clust)){
    
        clusters <- append(clusters, clust[x])
        d <- marker.sub[marker.sub$Cluster == clust[x], ]
        numCNV <- append(numCNV, sum(d$ChromosomeArm %in% var.reg == TRUE))
        #no <- sum(d$ChromosomeArm %in% var.reg == FALSE)
        numMarkers <- append(numMarkers, nrow(d))
        
         ### do fishers exact test
        topleft <- sum(d$ChromosomeArm %in% var.reg == TRUE)
        bottomleft <- nrow(d) - topleft
        genes.in.var <- genepos[as.character(genepos$arm) %in% var.reg, ]
        topright <- as.numeric(table(genes.in.var$Gene %in% d$Gene)["FALSE"]) 

        #not marker gene and not in var region
        genes.NOT.in.var <- genepos[!as.character(genepos$arm) %in% var.reg, ]
        bottomright <- as.numeric(table(genes.NOT.in.var$Gene %in% d$Gene)["FALSE"]) 
    
        con <- matrix(c(topleft, bottomleft, topright, bottomright), ncol = 2)
        fishers <- append(fishers, fisher.test(con, alternative = "greater")$p)
        #fishers <- append(fishers, chisq.test(con, alternative = "greater")$p)
    }

    ## calculate fishers p
    
}
```

    [1] "7p" "8q" "9p" "9q"
    [1] "9p"  "20p" "19p"
    [1] "7p"  "9p"  "19p" "22q"
    [1] "6p"  "6q"  "7p"  "18q"
    [1] "7p"  "17q" "20q"
    [1] "6q"  "17p" "17q" "19q"
    [1] "8q"  "19p"
    [1] "4p"  "7p"  "10p" "10q" "12q" "13q" "19p" "19q"
     [1] "1q"  "6p"  "6q"  "7p"  "9q"  "10p" "10q" "13q" "16q" "20p" "19p" "19q"
    [13] "21q"
    [1] "7q"  "10p" "10q" "13q" "20p" "20q" "19p" "19q"
    [1] "7p"  "7q"  "18q" "19q"
     [1] "3p"  "4p"  "4q"  "7p"  "7q"  "8q"  "9p"  "14q" "20p" "20q" "19p" "19q"
    [1] "7q"  "9q"  "13q" "19p" "19q"
     [1] "1q"  "2p"  "2q"  "4p"  "4q"  "7p"  "9q"  "11p" "13q" "16q" "18q" "20p"
    [13] "20q" "19p" "19q" "21q"
    [1] "2p"  "10p" "19q"
    [1] "7p"  "9q"  "12p" "12q" "13q" "18q" "19q" "21q"
    [1] "6p"  "6q"  "9p"  "10p" "10q" "12q"
    [1] "9q"  "18q" "19p"
    [1] "7q"  "10p" "18q"
    [1] "7q"  "17p" "19q"
    [1] "5p"  "6q"  "10p" "10q" "15q" "17p" "17q" "22q"
     [1] "4p"  "4q"  "7p"  "7q"  "9p"  "10q" "16q" "17p" "20p" "20q" "19p" "22q"



```R
plot.dat <- data.frame(clusters, numMarkers, numCNV, fishers)
#remove <- c("BT147_L", 
#            "BT67_L",
#            "BT84_L", 
#            "BT89_L", 
#            "G566_L", 
#            "G799_L", 
#            "G946-K_L"
#           ) #samples with no variable regions
#plot.dat$Sample <- gsub('.{0,3}$', '', plot.dat$clusters)
#plot.dat <- plot.dat[!plot.dat$Sample %in% remove, ]

plot.dat$propCNV <- plot.dat$numCNV / plot.dat$numMarkers
head(plot.dat)

```


<table>
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>clusters</th><th scope=col>numMarkers</th><th scope=col>numCNV</th><th scope=col>fishers</th><th scope=col>propCNV</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>BT127_L_C1</td><td>3307</td><td>339</td><td>3.155183e-06</td><td>0.10250983</td></tr>
	<tr><th scope=row>2</th><td>BT127_L_C2</td><td> 727</td><td> 84</td><td>7.349326e-04</td><td>0.11554333</td></tr>
	<tr><th scope=row>3</th><td>BT48_L_C1 </td><td> 668</td><td> 49</td><td>2.141111e-03</td><td>0.07335329</td></tr>
	<tr><th scope=row>4</th><td>BT48_L_C2 </td><td> 581</td><td> 18</td><td>9.832090e-01</td><td>0.03098107</td></tr>
	<tr><th scope=row>5</th><td>BT73_L_C1 </td><td>3226</td><td>319</td><td>3.919370e-05</td><td>0.09888407</td></tr>
	<tr><th scope=row>6</th><td>BT73_L_C2 </td><td> 264</td><td> 35</td><td>2.374066e-03</td><td>0.13257576</td></tr>
</tbody>
</table>




```R
cols <- ifelse(fishers < 0.05, "darkblue", "grey")
pdf("~/Desktop/propMarkersCNV.pdf", height = 22, width = 5)
par(oma=c(4,4,0,0))
barplot(plot.dat$propCNV, 
        names = plot.dat$clusters, 
        las = 2, 
        col=cols,
        horiz = T,
        cex.names = 1,
        xlim = c(0,0.6),
        xlab = "Proportion Marker Genes",
        #ylab = "Cluster",
        xpd = F
       )
dev.off()
```


<strong>pdf:</strong> 2

#marker gene and in var region
topleft <- plot.dat$numCNV[1] 

#marker gene and NOT in var region
bottomleft <- plot.dat$numMarkers[1] - topleft 

#not marker gene but in var region
genes.in.var <- genepos[as.character(genepos$arm) %in% var.reg, ]
topright <- as.numeric(table(genes.in.var$Gene %in% d$Gene)["FALSE"]) 

#not marker gene and not in var region
genes.NOT.in.var <- genepos[!as.character(genepos$arm) %in% var.reg, ]
bottomright <- as.numeric(table(genes.NOT.in.var$Gene %in% d$Gene)["FALSE"]) 

con <- matrix(c(topleft, bottomleft, topright, bottomright), ncol = 2)

fishers <- append(fishers, fisher.test(con, alternative = "greater")$p)