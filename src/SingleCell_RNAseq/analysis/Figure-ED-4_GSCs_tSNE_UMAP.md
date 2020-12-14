
----
# Visualize GSC Cohort
---

L.Richards
#library("Seurat", 
#        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" 
#       )


load("Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
BTSC_cellcoords <- cbind(BTSC@meta.data, 
                         BTSC@dr$tsne@cell.embeddings,
                        BTSC@dr$umap@cell.embeddings
                        )
save(BTSC_cellcoords, file = "Global_SU2C_BTSCs_CCregressed_noRibo_30PCS_umap.Rdata")

pdf("Seurat_umap.pdf")
TSNEPlot(object = BTSC, 
           do.label = TRUE, 
            pt.size = 0.2,
            group.by = "SampleID"
            
            )
DimPlot(object = BTSC, 
        reduction.use = "umap",
           do.label = TRUE, 
            pt.size = 0.2,
            group.by = "SampleID"
            
            )            
            
TSNEPlot(object = BTSC, 
           do.label = TRUE, 
            pt.size = 0.2
            
            )
            
            DimPlot(object = BTSC, 
        reduction.use = "umap",
           do.label = TRUE, 
            pt.size = 0.2)
            
           
dev.off()

```R
library(ggplot2)
```


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/data//SeuratObj/BTSCs/")
```


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo_30PCS_umap.Rdata")
```


```R
BTSC_cellcoords$PatientID <- as.character(BTSC_cellcoords$PatientID)
table(BTSC_cellcoords$SampleID)
colnames(BTSC_cellcoords)
```


    
     BT127_L  BT147_L   BT48_L   BT67_L   BT73_L   BT84_L   BT89_L   BT94_L 
        1444      862     1467     1202     1649     1796      893     1223 
      G523_L   G549_L   G564_L   G566_L   G583_L   G620_L   G637_L   G729_L 
        3745     3434     3707     2142     3502     3168     3526     2881 
      G797_L   G799_L   G800_L   G837_L   G851_L   G876_L   G885_L   G895_L 
        1160     1062     3738     4687     1501      834     1048     1308 
    G945-I_L G945-J_L G945-K_L G946-J_L G946-K_L 
        5229     5189     3681     1691     1624 



<ol class=list-inline>
	<li>'nGene'</li>
	<li>'nUMI'</li>
	<li>'percent.mito'</li>
	<li>'PatientID'</li>
	<li>'SampleID'</li>
	<li>'Passage'</li>
	<li>'CultureMethod'</li>
	<li>'Pathology'</li>
	<li>'Stage'</li>
	<li>'Age'</li>
	<li>'Sex'</li>
	<li>'IDH12_Status'</li>
	<li>'SphereFormation_Percent'</li>
	<li>'DoublingTime_Hrs'</li>
	<li>'OrthotopicXenoSurvival_Days'</li>
	<li>'BulkRNA_Cluster'</li>
	<li>'res.2'</li>
	<li>'S.Score'</li>
	<li>'G2M.Score'</li>
	<li>'Phase'</li>
	<li>'CC.Difference'</li>
	<li>'Cluster.ID'</li>
	<li>'Opt.Resolution'</li>
	<li>'IntraBTSC.ID'</li>
	<li>'tSNE_1'</li>
	<li>'tSNE_2'</li>
	<li>'UMAP1'</li>
	<li>'UMAP2'</li>
</ol>




```R
BTSC_cellcoords$BulkRNA_Cluster[is.na(BTSC_cellcoords$BulkRNA_Cluster)] <- "ND"

```

    Warning message in `[<-.factor`(`*tmp*`, is.na(BTSC_cellcoords$BulkRNA_Cluster), :
    “invalid factor level, NA generated”


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
data.frame(names(table(BTSC_cellcoords$SampleID)), cols)

```


<table>
<thead><tr><th scope=col>names.table.BTSC_cellcoords.SampleID..</th><th scope=col>cols</th></tr></thead>
<tbody>
	<tr><td>BT127_L </td><td>#800026 </td></tr>
	<tr><td>BT147_L </td><td>#B52127 </td></tr>
	<tr><td>BT48_L  </td><td>#EA4229 </td></tr>
	<tr><td>BT67_L  </td><td>#FC6A33 </td></tr>
	<tr><td>BT73_L  </td><td>#FD9542 </td></tr>
	<tr><td>BT84_L  </td><td>#FEBA57 </td></tr>
	<tr><td>BT89_L  </td><td>#FED37B </td></tr>
	<tr><td>BT94_L  </td><td>#FFEDA0 </td></tr>
	<tr><td>G523_L  </td><td>#54278F </td></tr>
	<tr><td>G549_L  </td><td>#6E4CA2 </td></tr>
	<tr><td>G564_L  </td><td>#8872B5 </td></tr>
	<tr><td>G566_L  </td><td>#A297C8 </td></tr>
	<tr><td>G583_L  </td><td>#BCBDDC </td></tr>
	<tr><td>G620_L  </td><td>#8F9DC5 </td></tr>
	<tr><td>G637_L  </td><td>#617EAE </td></tr>
	<tr><td>G729_L  </td><td>#345F97 </td></tr>
	<tr><td>G797_L  </td><td>#084081 </td></tr>
	<tr><td>G799_L  </td><td>#195C95 </td></tr>
	<tr><td>G800_L  </td><td>#2A79A9 </td></tr>
	<tr><td>G837_L  </td><td>#3C96BE </td></tr>
	<tr><td>G851_L  </td><td>#4EB3D3 </td></tr>
	<tr><td>G876_L  </td><td>#43A9AF </td></tr>
	<tr><td>G885_L  </td><td>#389F8B </td></tr>
	<tr><td>G895_L  </td><td>#2D9568 </td></tr>
	<tr><td>G945-I_L</td><td>#238B45 </td></tr>
	<tr><td>G945-J_L</td><td>#4DA365 </td></tr>
	<tr><td>G945-K_L</td><td>#77BB84 </td></tr>
	<tr><td>G946-J_L</td><td>#A1D3A5 </td></tr>
	<tr><td>G946-K_L</td><td>#CCEBC5 </td></tr>
</tbody>
</table>




```R
sample_tSNE <- ggplot(BTSC_cellcoords, aes(x=tSNE_1, y=tSNE_2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "tSNE 1", y = "tSNE 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

sample_umap <- ggplot(BTSC_cellcoords, aes(x=UMAP1, y=UMAP2, color=SampleID)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

cluster_umap <- ggplot(BTSC_cellcoords, aes(x=UMAP1, y=UMAP2, color=res.2)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))


cluster_tsne <- ggplot(BTSC_cellcoords, aes(x=tSNE_1, y=tSNE_2, color=res.2)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "tSNE 1", y = "tSNE 2") +
                   #scale_colour_manual(values = cols) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

bulkRNA_tsne <- ggplot(BTSC_cellcoords, aes(x=tSNE_1, y=tSNE_2, color=BulkRNA_Cluster)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "tSNE 1", y = "tSNE 2") +
                   scale_colour_manual(values = c("black", "darkred", "grey")) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

bulkRNA_umap <- ggplot(BTSC_cellcoords, aes(x=UMAP1, y=UMAP2, color=BulkRNA_Cluster)) + 
                   geom_point(alpha = 0.3, size = 0.8, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2") +
                   scale_colour_manual(values = c("black", "darkred", "grey")) + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))
```


```R
pdf("~/Desktop/NewCohortBTSC_tSNE_uamp.pdf", width = 8, height = 6.5)

sample_tSNE
sample_umap
cluster_tsne
cluster_umap
bulkRNA_tsne
bulkRNA_umap

dev.off()
```









    Warning message:
    “Removed 18858 rows containing missing values (geom_point).”



    Warning message:
    “Removed 18858 rows containing missing values (geom_point).”




<strong>pdf:</strong> 2


---
## Cluster distribution
---


```R
load("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo_30PCS_umap.Rdata")
```


```R
df <- prop.table(table(BTSC_cellcoords$res.2, BTSC_cellcoords$SampleID), margin = 1) 
```


```R
write.csv(df, file = "~/Desktop/BTSC_clusterProportion.csv")
```


```R

```
