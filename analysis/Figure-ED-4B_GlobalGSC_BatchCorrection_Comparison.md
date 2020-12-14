
----
# Visualize results from batch correction comparison in GSCS
---

L.Richards


```R
options(repos='http://cran.rstudio.com/')
#install.packages("ggExtra")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggrepel)
library(dplyr)
```


```R
setwd("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/")
```

----
## 1.0 Plot UMAPs
----


```R
#load data
meta <- readRDS("Global_GSC_BatchCorrection_metadata.rds")
```


```R
colnames(meta)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'nGene'</li><li>'nUMI'</li><li>'percent.mito'</li><li>'PatientID'</li><li>'SampleID'</li><li>'Passage'</li><li>'CultureMethod'</li><li>'Pathology'</li><li>'Stage'</li><li>'Age'</li><li>'Sex'</li><li>'IDH12_Status'</li><li>'SphereFormation_Percent'</li><li>'DoublingTime_Hrs'</li><li>'OrthotopicXenoSurvival_Days'</li><li>'BulkRNA_Cluster'</li><li>'S.Score'</li><li>'G2M.Score'</li><li>'Phase'</li><li>'CC.Difference'</li><li>'Cluster.ID'</li><li>'Opt.Resolution'</li><li>'IntraBTSC.ID'</li><li>'Original_clusters'</li><li>'Original_UMAP1'</li><li>'Original_UMAP2'</li><li>'Conos_clusters'</li><li>'Conos_UMAP1'</li><li>'Conos_UMAP2'</li><li>'Liger_clusters'</li><li>'Liger_quantNorm_clusters'</li><li>'Liger_UMAP1'</li><li>'Liger_UMAP2'</li><li>'fastMNN_clusters'</li><li>'fastMNN_UMAP1'</li><li>'fastMNN_UMAP2'</li></ol>




```R
## Define sample color palette

colfunc <- colorRampPalette(c("#54278f", "#bcbddc", "#084081", "#4eb3d3", "#238b45", "#ccebc5"))
dirks <- colfunc(21)
colfunc <- colorRampPalette(c("#800026", "#fc4e2a", "#feb24c", "#ffeda0"))
weiss <- colfunc(8)
cols <- c(weiss, dirks)
length(cols)
```


29



```R
#### ORIGINAL CLUSTERING


#calculate centroids 
hc.norm.cent <- meta %>% group_by(Original_clusters) %>% select(Original_UMAP1, 
    Original_UMAP2) %>% summarize_all(median)
#hc.norm.cent



original_sample <- ggplot(meta, aes(x=Original_UMAP1, y=Original_UMAP2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="none") 

original_clusters <- ggplot(meta, aes(x=Original_UMAP1, y=Original_UMAP2, color=Original_clusters)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="none") +
                    geom_label_repel(aes(label = Original_clusters), 
                                     data = hc.norm.cent, 
                                     label.size = 0.05, 
                                     parse = T, 
                                     size = 3)

pdf("~/Desktop/OriginalClustering_NatCan_Cluster.pdf", width = 8, height= 8)
original_clusters
dev.off()

pdf("~/Desktop/OriginalClustering_NatCan_Sample.pdf", width = 8, height= 8)
original_sample
dev.off()
```

    Adding missing grouping variables: `Original_clusters`
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



```R
#### CONOS Clustering


#calculate centroids 
hc.norm.cent <- meta %>% group_by(Conos_clusters) %>% select(Conos_UMAP1, 
    Conos_UMAP2) %>% summarize_all(median)
#hc.norm.cent



Conos_sample <- ggplot(meta, aes(x=Conos_UMAP1, y=Conos_UMAP2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="bottom") + theme(legend.position="none")

Conos_clusters <- ggplot(meta, aes(x=Conos_UMAP1, y=Conos_UMAP2, color=Conos_clusters)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="none") +
                    geom_label_repel(aes(label = Conos_clusters), 
                                     data = hc.norm.cent, 
                                     label.size = 0.05, 
                                     parse = T, 
                                     size = 3)

pdf("~/Desktop/ConosClustering_NatCan_Cluster.pdf", width = 8, height= 8)
Conos_clusters
dev.off()

pdf("~/Desktop/ConosClustering_NatCan_Sample.pdf", width = 8, height= 8)
Conos_sample
dev.off()
```

    Adding missing grouping variables: `Conos_clusters`
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



```R
#### Liger Clustering


#calculate centroids 
hc.norm.cent <- meta %>% group_by(Liger_clusters) %>% select(Liger_UMAP1, 
    Liger_UMAP2) %>% summarize_all(median)
#hc.norm.cent



Liger_sample <- ggplot(meta, aes(x=Liger_UMAP1, y=Liger_UMAP2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="bottom") + theme(legend.position="none")

Liger_clusters <- ggplot(meta, aes(x=Liger_UMAP1, y=Liger_UMAP2, color=Liger_clusters)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="none") +
                    geom_label_repel(aes(label = Liger_clusters), 
                                     data = hc.norm.cent, 
                                     label.size = 0.05, 
                                     parse = T, 
                                     size = 3)

pdf("~/Desktop/LigerClustering_NatCan_Cluster.pdf", width = 8, height= 8)
Liger_clusters
dev.off()

pdf("~/Desktop/LigerClustering_NatCan_Sample.pdf", width = 8, height= 8)
Liger_sample
dev.off()
```

    Adding missing grouping variables: `Liger_clusters`
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



```R
#### fastMNN Clustering


#calculate centroids 
hc.norm.cent <- meta %>% group_by(fastMNN_clusters) %>% select(fastMNN_UMAP1, 
    fastMNN_UMAP2) %>% summarize_all(median)
#hc.norm.cent



fastMNN_sample <- ggplot(meta, aes(x=fastMNN_UMAP1, y=fastMNN_UMAP2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="bottom") + theme(legend.position="none")

fastMNN_clusters <- ggplot(meta, aes(x=fastMNN_UMAP1, y=fastMNN_UMAP2, color=fastMNN_clusters)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
                    theme(legend.position="none") +
                    geom_label_repel(aes(label = fastMNN_clusters), 
                                     data = hc.norm.cent, 
                                     label.size = 0.05, 
                                     parse = T, 
                                     size = 3)

pdf("~/Desktop/fastMNNClustering_NatCan_Cluster.pdf", width = 8, height= 8)
fastMNN_clusters
dev.off()

pdf("~/Desktop/fastMNNClustering_NatCan_Sample.pdf", width = 8, height= 8)
fastMNN_sample
dev.off()
```

    Adding missing grouping variables: `fastMNN_clusters`
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2


----
## 2.0 Plot cluster properties
----


```R
setwd("./Comparison/")
```


```R
#### load data
meta <- readRDS("./Comparison/Global_GSC_BatchCorrection_metadata.rds")
head(meta)
```


<table>
<caption>A data.frame: 6 × 36</caption>
<thead>
	<tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>percent.mito</th><th scope=col>PatientID</th><th scope=col>SampleID</th><th scope=col>Passage</th><th scope=col>CultureMethod</th><th scope=col>Pathology</th><th scope=col>Stage</th><th scope=col>Age</th><th scope=col>⋯</th><th scope=col>Conos_clusters</th><th scope=col>Conos_UMAP1</th><th scope=col>Conos_UMAP2</th><th scope=col>Liger_clusters</th><th scope=col>Liger_quantNorm_clusters</th><th scope=col>Liger_UMAP1</th><th scope=col>Liger_UMAP2</th><th scope=col>fastMNN_clusters</th><th scope=col>fastMNN_UMAP1</th><th scope=col>fastMNN_UMAP2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td> 640</td><td>  875</td><td>0.043428571</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C11</td><td> 3.91472455</td><td>-2.492949</td><td>C67</td><td>C13</td><td>-0.4835586</td><td> 2.379755</td><td>C37</td><td>6.7959390</td><td> 4.973665</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>1036</td><td> 2408</td><td>0.002076412</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C2 </td><td>-2.32181506</td><td> 1.148725</td><td>C14</td><td>C21</td><td> 0.9863441</td><td> 5.408310</td><td>C14</td><td>1.2473862</td><td>-1.664298</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>3240</td><td>10058</td><td>0.078047326</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C1 </td><td>-0.07283644</td><td>-2.741365</td><td>C56</td><td>C10</td><td> 8.9876844</td><td> 1.797849</td><td>C10</td><td>2.2776213</td><td> 1.684389</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>3337</td><td>10798</td><td>0.061863308</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C1 </td><td> 2.88304897</td><td>-4.185107</td><td>C28</td><td>C7 </td><td>-3.7655035</td><td>-4.287652</td><td>C21</td><td>0.3696045</td><td> 4.849369</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>4140</td><td>14601</td><td>0.081501267</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C4 </td><td>-3.40395693</td><td>-1.246072</td><td>C8 </td><td>C6 </td><td>-2.0646044</td><td> 5.311308</td><td>C3 </td><td>3.0244613</td><td>-2.603526</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 543</td><td>  820</td><td>0.108536585</td><td>BT127</td><td>BT127_L</td><td>NA</td><td>Sphere</td><td>GLIOBLASTOMA (GRADE IV)</td><td>PRIMARY</td><td>55</td><td>⋯</td><td>C6 </td><td> 3.92768212</td><td>-2.285280</td><td>C7 </td><td>C17</td><td> 4.2134230</td><td>-1.124725</td><td>C20</td><td>6.4704309</td><td> 4.624188</td></tr>
</tbody>
</table>



---
### Number clusters


```R
##### NUMBER CLUSTERS
##### HISTOGRAM

clusters <- c(61, 12, 78, 39)
method <- c("Original", "Conos", "Liger", "fastMNN")
clust.hist <- cbind(method, clusters)
colnames(clust.hist) <- c("method", "clusters")
clust.hist <- data.frame(clust.hist)
clust.hist$method <- factor(clust.hist$method, 
                               levels = c("Original", "Conos", "Liger", "fastMNN")
                              )

pdf("ClusterCounts.pdf", width = 5, height = 5)
cluster.plot <- ggplot(clust.hist, aes(y=clusters, x=method, fill = method)) +
                geom_bar(stat="identity") + theme_classic() +
                ylab("Number of Clusters") + xlab("") +
                theme(legend.position="none") +
                theme(text = element_text(size=20))

cluster.plot
dev.off()
```


<strong>pdf:</strong> 2


---
### Number samples per cluster


```R
#### NUMBER SAMPLES / CLUSTER (>10 cells) 
#### BOXPLOT WITH POINTS

original <- table(meta$Original_clusters, meta$SampleID)
original <- original > 10
original <- rowSums(original)
original <- data.frame(original)
colnames(original) <- "Samples"
original$method <- "Original"
original$ClusterID <- rownames(original)
print(mean(original$Samples))

conos <- table(meta$Conos_clusters, meta$SampleID)
conos <- conos > 10
conos <- rowSums(conos)
conos <- data.frame(conos)
colnames(conos) <- "Samples"
conos$method <- "Conos"
conos$ClusterID <- rownames(conos)
print(mean(conos$Samples))

liger <- table(meta$Liger_clusters, meta$SampleID)
liger <- liger > 10
liger <- rowSums(liger)
liger <- data.frame(liger)
colnames(liger) <- "Samples"
liger$method <- "Liger"
liger$ClusterID <- rownames(liger)
print(mean(liger$Samples))


fastmnn <- table(meta$fastMNN_clusters, meta$SampleID)
fastmnn <- fastmnn > 10
fastmnn <- rowSums(fastmnn)
fastmnn<- data.frame(fastmnn)
colnames(fastmnn) <- "Samples"
fastmnn$method <- "fastMNN"
fastmnn$ClusterID <- rownames(fastmnn)
print(mean(fastmnn$Samples))

plot.dat <- rbind(original, conos, liger, fastmnn)
plot.dat$method <- factor(plot.dat$method, 
                               levels = c("Original", "Conos", "Liger", "fastMNN")
                              )
head(plot.dat)
```

    [1] 1.229508
    [1] 17.41667
    [1] 19.78205
    [1] 9.769231



<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Samples</th><th scope=col>method</th><th scope=col>ClusterID</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>C1</th><td>1</td><td>Original</td><td>C1 </td></tr>
	<tr><th scope=row>C10</th><td>1</td><td>Original</td><td>C10</td></tr>
	<tr><th scope=row>C11</th><td>1</td><td>Original</td><td>C11</td></tr>
	<tr><th scope=row>C12</th><td>1</td><td>Original</td><td>C12</td></tr>
	<tr><th scope=row>C13</th><td>2</td><td>Original</td><td>C13</td></tr>
	<tr><th scope=row>C14</th><td>1</td><td>Original</td><td>C14</td></tr>
</tbody>
</table>




```R
pdf("ClusterSampleBoxplot.pdf", width = 5, height = 5)
p <- ggplot(plot.dat, aes(x=method, y=Samples, fill = method)) + 
      geom_boxplot() + theme_classic() + geom_jitter(shape=16, position=position_jitter(0.2)) +
      ylab("Samples with > 10 cells / cluster") + xlab("") + 
        theme(legend.position="none") + theme(text = element_text(size=20))
p
dev.off()
```


<strong>pdf:</strong> 2


---
### Proportion samples per cluster


```R
colfunc <- colorRampPalette(c("#54278f", "#bcbddc", "#084081", "#4eb3d3", "#238b45", "#ccebc5"))
dirks <- colfunc(21)
colfunc <- colorRampPalette(c("#800026", "#fc4e2a", "#feb24c", "#ffeda0"))
weiss <- colfunc(8)
cols <- c(weiss, dirks)
length(cols)
```


```R
original <- prop.table(table(meta$SampleID, meta$Original_clusters), margin = 2)
original <- data.matrix(original)
original <- data.frame(original)
colnames(original) <- c("Sample", "Cluster", "Prop")

pdf("Original_ClusterProportions.pdf", height = 8, width = 23)
p1 <- ggplot() + geom_bar(aes(y=Prop, x=Cluster, fill = Sample), data = original, stat = "identity") +
     scale_fill_manual(values = cols) + theme_classic() + ggtitle("Original") + 
     xlab("") + ylab("Proportion Cluster")+ 
     theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))
p1
dev.off()
```


<strong>pdf:</strong> 2



```R
conos <- prop.table(table(meta$SampleID, meta$Conos_clusters), margin = 2)
conos <- data.matrix(conos)
conos <- data.frame(conos)
colnames(conos) <- c("Sample", "Cluster", "Prop")

pdf("Conos_ClusterProportions.pdf", height = 8, width = 23)
p1 <- ggplot() + geom_bar(aes(y=Prop, x=Cluster, fill = Sample), data = conos, stat = "identity") +
     scale_fill_manual(values = cols) + theme_classic() + ggtitle("Conos") + 
     xlab("") + ylab("Proportion Cluster") + 
     theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))
p1
dev.off()
```


<strong>pdf:</strong> 2



```R
liger <- prop.table(table(meta$SampleID, meta$Liger_clusters), margin = 2)
liger<- data.matrix(liger)
liger <- data.frame(liger)
colnames(liger) <- c("Sample", "Cluster", "Prop")

pdf("Liger_ClusterProportions.pdf", height = 8, width = 23)
p1 <- ggplot() + geom_bar(aes(y=Prop, x=Cluster, fill = Sample), data = liger, stat = "identity") +
     scale_fill_manual(values = cols) + theme_classic() + ggtitle("Liger") + 
     xlab("") + ylab("Proportion Cluster") + 
     theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))
p1
dev.off()

```


<strong>pdf:</strong> 2



```R
fastmnn <- prop.table(table(meta$SampleID, meta$fastMNN_clusters), margin = 2)
fastmnn<- data.matrix(fastmnn)
fastmnn <- data.frame(fastmnn)
colnames(fastmnn) <- c("Sample", "Cluster", "Prop")

pdf("fastMNN_ClusterProportions.pdf", height = 8, width = 23)
p1 <- ggplot() + geom_bar(aes(y=Prop, x=Cluster, fill = Sample), data = fastmnn, stat = "identity") +
     scale_fill_manual(values = cols) + theme_classic() + ggtitle("fastMNN") + 
     xlab("") + ylab("Proportion Cluster") + 
     theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))
p1
dev.off()
```


<strong>pdf:</strong> 2



```R

```
