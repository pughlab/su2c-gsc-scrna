
----
## Replot CRISPR Figure with new essential genes
---

L.Richards  
Previous matrix from Graham had var genes defined in a werid way with more samples, we decided to redo the anlaysis by subsetting the 

---
## 1.0 Plot gene overlaps across samples
---

Use BF > 10


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/")
getwd()
```


'/Users/laura/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR'



```R
crispr <- readRDS("TKOv3_GBM_qBF_LRichards.rds")
```


```R
xx <- rowSums(crispr > 10)
head(xx)

counts <- table(xx)
counts


```


<dl class=dl-horizontal>
	<dt>ITGB1</dt>
		<dd>7</dd>
	<dt>SCAP</dt>
		<dd>8</dd>
	<dt>EXOSC5</dt>
		<dd>8</dd>
	<dt>SUPT16H</dt>
		<dd>8</dd>
	<dt>PPP1R12A</dt>
		<dd>7</dd>
	<dt>AARS</dt>
		<dd>8</dd>
</dl>




    xx
        0     1     2     3     4     5     6     7     8     9    10    11 
    14172  1608   464   287   204   177   144   169   171   195   223   234 



```R
sum(counts[-c(1,2,3,11,12)])
```


1347



```R
pdf("~/Desktop/EssentialGene_ScreenOverlap.pdf", width = 6, height = 4)

barplot(counts[-1], 
        col = "grey",
        ylab = "Essential Genes (qBF > 10)",
        xlab = "Number of Screens"
       )

dev.off()
```


<strong>pdf:</strong> 2



```R
### how many essential genes in each screen 
aa <- colSums(crispr > 10)
mean(aa)

```


1370.18181818182


---
## 2.0 Replot the Var Gene and Corr heatmap
---


```R
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(circlize)
```


```R
var.genes <- read.table("~/Downloads/subset_SU2C_11lines_BF10_3to9.txt",
                        header = T
                       )

rownames(var.genes) <- var.genes$Gene
var.genes$Gene <- NULL
head(var.genes)
colnames(var.genes) <- c("BT67_L",
                         "G620_L",
                         "G523_L",
                         "G809_L",
                         "G583_L",
                         "G549_L",
                         "G729_L",
                         "G361_L",
                          "G691_L",
                          "G440_L",
                         "G532_L"
                        )
var.genes <- na.omit(var.genes)

scaled.dat <- scale(t(var.genes))
dat <- t(scaled.dat)
head(dat)
```


<table>
<thead><tr><th></th><th scope=col>BT67</th><th scope=col>G620</th><th scope=col>G523</th><th scope=col>G809R</th><th scope=col>G583NS</th><th scope=col>G549NS</th><th scope=col>G729</th><th scope=col>G361NS</th><th scope=col>G691R</th><th scope=col>G440</th><th scope=col>G532R</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>40.43901   </td><td>40.93868   </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td>33.67403   </td><td>37.50733   </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td>33.63656   </td><td>47.18051   </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>37.95362   </td><td>53.55119   </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td>14.84797   </td><td>28.55795   </td></tr>
	<tr><th scope=row>AARS</th><td>  0.9270344</td><td>  4.774409 </td><td>-10.0122166</td><td> 34.759483 </td><td>52.83310   </td><td>47.83145   </td><td>37.41367   </td><td>39.54539   </td><td>28.55795   </td><td>47.44578   </td><td>46.42611   </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-1.6490091 </td><td>-1.157435  </td><td>-1.319083  </td><td>-0.7377261 </td><td>0.4558661  </td><td>0.7716652  </td><td>0.5354065  </td><td>0.53125143 </td><td> 0.61899187</td><td> 0.9655297 </td><td>0.9845428  </td></tr>
	<tr><th scope=row>SCAP</th><td> 0.5183667 </td><td>-1.190939  </td><td>-1.682580  </td><td>-1.5323825 </td><td>0.7273251  </td><td>1.0041362  </td><td>0.8911478  </td><td>0.05043592 </td><td>-0.04127467</td><td> 0.5588545 </td><td>0.6969109  </td></tr>
	<tr><th scope=row>EXOSC5</th><td>-1.3470887 </td><td>-1.543471  </td><td>-1.501585  </td><td>-0.1137754 </td><td>0.9253870  </td><td>0.2497787  </td><td>0.7368359  </td><td>0.93902287 </td><td> 0.29093071</td><td> 0.3844250 </td><td>0.9795407  </td></tr>
	<tr><th scope=row>SUPT16H</th><td>-1.5049491 </td><td>-1.432480  </td><td>-1.418158  </td><td>-0.3630855 </td><td>0.4637230  </td><td>0.5671886  </td><td>0.7504430  </td><td>0.62337052 </td><td> 0.84954456</td><td> 0.3732704 </td><td>1.0911328  </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>-0.7878167 </td><td>-1.010871  </td><td>-1.515716  </td><td>-1.0741564 </td><td>0.7422234  </td><td>1.6072348  </td><td>0.6924964  </td><td>0.02810969 </td><td> 1.03810245</td><td>-0.1554960 </td><td>0.4358899  </td></tr>
	<tr><th scope=row>AARS</th><td>-1.3479086 </td><td>-1.169813  </td><td>-1.854290  </td><td> 0.2182046 </td><td>1.0548373  </td><td>0.8233094  </td><td>0.3410678  </td><td>0.43974532 </td><td>-0.06886593</td><td> 0.8054566 </td><td>0.7582557  </td></tr>
</tbody>
</table>




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
pdf("~/Desktop/Cluster_CRISPR_VarGenes_March2020.pdf", width = 7, height = 7)

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
        row_title = "Variable Genes (n=1345)",
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


dev.off()


```




<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th></tr></thead>
<tbody>
	<tr><th scope=row>BT67_L</th><td>1.00000000 </td><td> 0.38660887</td><td>0.17814257 </td><td>0.2583741  </td><td> 0.01983397</td><td>0.12551259 </td><td> 0.07658890</td><td>0.18919020 </td><td> 0.05530650</td><td>0.084204001</td><td> 0.07857563</td></tr>
	<tr><th scope=row>G620_L</th><td>0.38660887 </td><td> 1.00000000</td><td>0.12826996 </td><td>0.2278008  </td><td>-0.01174671</td><td>0.13338704 </td><td>-0.01990205</td><td>0.09477192 </td><td> 0.06723074</td><td>0.061190548</td><td>-0.01731482</td></tr>
	<tr><th scope=row>G523_L</th><td>0.17814257 </td><td> 0.12826996</td><td>1.00000000 </td><td>0.1399938  </td><td> 0.06592393</td><td>0.04094582 </td><td> 0.03449849</td><td>0.06203901 </td><td>-0.02332991</td><td>0.001131441</td><td> 0.02287645</td></tr>
	<tr><th scope=row>G809_L</th><td>0.25837407 </td><td> 0.22780078</td><td>0.13999379 </td><td>1.0000000  </td><td> 0.12575515</td><td>0.13247926 </td><td> 0.16921440</td><td>0.20869279 </td><td> 0.03474607</td><td>0.152443863</td><td> 0.17755113</td></tr>
	<tr><th scope=row>G583_L</th><td>0.01983397 </td><td>-0.01174671</td><td>0.06592393 </td><td>0.1257551  </td><td> 1.00000000</td><td>0.25786178 </td><td> 0.36892255</td><td>0.39804289 </td><td> 0.11150215</td><td>0.304134054</td><td> 0.35835544</td></tr>
	<tr><th scope=row>G549_L</th><td>0.12551259 </td><td> 0.13338704</td><td>0.04094582 </td><td>0.1324793  </td><td> 0.25786178</td><td>1.00000000 </td><td> 0.41695998</td><td>0.33320210 </td><td> 0.14681581</td><td>0.331155288</td><td> 0.39741855</td></tr>
</tbody>
</table>




<strong>pdf:</strong> 2



```R
pdf("~/Desktop/Cluster_CRISPR_Correlation_MArch2020.pdf", width = 7, height = 5)

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

corr.map
dev.off()
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th></tr></thead>
<tbody>
	<tr><th scope=row>BT67_L</th><td>1.00000000 </td><td> 0.38660887</td><td>0.17814257 </td><td>0.2583741  </td><td> 0.01983397</td><td>0.12551259 </td><td> 0.07658890</td><td>0.18919020 </td><td> 0.05530650</td><td>0.084204001</td><td> 0.07857563</td></tr>
	<tr><th scope=row>G620_L</th><td>0.38660887 </td><td> 1.00000000</td><td>0.12826996 </td><td>0.2278008  </td><td>-0.01174671</td><td>0.13338704 </td><td>-0.01990205</td><td>0.09477192 </td><td> 0.06723074</td><td>0.061190548</td><td>-0.01731482</td></tr>
	<tr><th scope=row>G523_L</th><td>0.17814257 </td><td> 0.12826996</td><td>1.00000000 </td><td>0.1399938  </td><td> 0.06592393</td><td>0.04094582 </td><td> 0.03449849</td><td>0.06203901 </td><td>-0.02332991</td><td>0.001131441</td><td> 0.02287645</td></tr>
	<tr><th scope=row>G809_L</th><td>0.25837407 </td><td> 0.22780078</td><td>0.13999379 </td><td>1.0000000  </td><td> 0.12575515</td><td>0.13247926 </td><td> 0.16921440</td><td>0.20869279 </td><td> 0.03474607</td><td>0.152443863</td><td> 0.17755113</td></tr>
	<tr><th scope=row>G583_L</th><td>0.01983397 </td><td>-0.01174671</td><td>0.06592393 </td><td>0.1257551  </td><td> 1.00000000</td><td>0.25786178 </td><td> 0.36892255</td><td>0.39804289 </td><td> 0.11150215</td><td>0.304134054</td><td> 0.35835544</td></tr>
	<tr><th scope=row>G549_L</th><td>0.12551259 </td><td> 0.13338704</td><td>0.04094582 </td><td>0.1324793  </td><td> 0.25786178</td><td>1.00000000 </td><td> 0.41695998</td><td>0.33320210 </td><td> 0.14681581</td><td>0.331155288</td><td> 0.39741855</td></tr>
</tbody>
</table>






<strong>pdf:</strong> 2


---
## 3.0 Output the appropriate supp table
---


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/")

TS8 <- readRDS("CRISPR_Avg_Diff_Matrix.rds")

head(TS8)

write.csv(TS8,
            file = "~/Desktop/SuppTable8_CRISPR_qBFs.csv", 
            col.names = T,
            row.names = T
           )
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th><th scope=col>Developmental_Average</th><th scope=col>InjuryResponse_Average</th><th scope=col>Diff</th><th scope=col>Diff_Z</th><th scope=col>FoldChange</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>40.43901   </td><td>40.93868   </td><td>-16.8872393</td><td>30.37648   </td><td>47.26371   </td><td>7.051752   </td><td> -1.798783 </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td>33.67403   </td><td>37.50733   </td><td> -8.8287541</td><td>32.77162   </td><td>41.60037   </td><td>6.206774   </td><td> -3.711919 </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td>33.63656   </td><td>47.18051   </td><td> -0.7493837</td><td>39.18887   </td><td>39.93825   </td><td>5.958783   </td><td>-52.294793 </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>37.95362   </td><td>53.55119   </td><td>  4.2116815</td><td>43.98490   </td><td>39.77322   </td><td>5.934161   </td><td> 10.443548 </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td>14.84797   </td><td>28.55795   </td><td> -6.9819679</td><td>37.50056   </td><td>44.48253   </td><td>6.636796   </td><td> -5.371059 </td></tr>
	<tr><th scope=row>AARS</th><td>  0.9270344</td><td>  4.774409 </td><td>-10.0122166</td><td> 34.759483 </td><td>52.83310   </td><td>47.83145   </td><td>37.41367   </td><td>39.54539   </td><td>28.55795   </td><td>47.44578   </td><td>46.42611   </td><td>  7.6121776</td><td>41.23631   </td><td>33.62414   </td><td>5.016710   </td><td>  5.417151 </td></tr>
</tbody>
</table>



    Warning message in write.csv(TS8, file = "~/Desktop/SuppTable8_CRISPR_qBFs.csv", :
    “attempt to set 'col.names' ignored”


```R
list.files("./CRISPR_data/")
```


<ol class=list-inline>
	<li>'CRISPR_data_Jan142020.zip'</li>
	<li>'README.txt'</li>
	<li>'TKOv3_GBM_qBF_subset_for_clustering.txt'</li>
	<li>'TKOv3_GBM_qBF.txt'</li>
	<li>'TKOv3_GBM_qBF[1].txt'</li>
	<li>'TKOv3_GBM_rawBF.txt'</li>
</ol>




```R
TS7 <- read.table("./CRISPR_data/TKOv3_GBM_rawBF.txt", header = T)
colnames(TS7) <- c("Gene", "G440_L", "G523_L", "G729_L", "G532_L", "G620_L", "BT67_L", "G361_L", "G691_L", "G809_L", "G549_L", "G583_L")
head(TS7)


```


<table>
<thead><tr><th scope=col>Gene</th><th scope=col>G440_L</th><th scope=col>G523_L</th><th scope=col>G729_L</th><th scope=col>G532_L</th><th scope=col>G620_L</th><th scope=col>BT67_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G809_L</th><th scope=col>G549_L</th><th scope=col>G583_L</th></tr></thead>
<tbody>
	<tr><td>A1BG   </td><td>-10.903</td><td>-6.396 </td><td>-12.614</td><td>-12.295</td><td>-5.526 </td><td>  4.725</td><td>-11.419</td><td>-10.883</td><td>-22.567</td><td>-22.771</td><td>-13.039</td></tr>
	<tr><td>A1CF   </td><td> -9.381</td><td>-2.805 </td><td>-26.876</td><td>-21.977</td><td>-7.734 </td><td>-16.442</td><td>-48.510</td><td>-14.388</td><td>-15.055</td><td> -7.267</td><td>-22.249</td></tr>
	<tr><td>A2M    </td><td> -3.824</td><td> 0.840 </td><td>-17.443</td><td>  2.003</td><td>-7.163 </td><td>  3.343</td><td>-45.626</td><td>-12.568</td><td>-19.626</td><td>-19.704</td><td>-23.258</td></tr>
	<tr><td>A2ML1  </td><td> -7.578</td><td>-7.185 </td><td>-31.760</td><td>-18.741</td><td>-5.462 </td><td>-15.391</td><td>-58.621</td><td> -8.590</td><td> -7.929</td><td>-17.880</td><td> -5.968</td></tr>
	<tr><td>A3GALT2</td><td> -4.173</td><td>-8.872 </td><td>-24.616</td><td>-25.149</td><td>-8.357 </td><td> -7.890</td><td>-45.236</td><td> -7.248</td><td> -6.438</td><td>-16.663</td><td>-23.153</td></tr>
	<tr><td>A4GALT </td><td> -8.087</td><td>-6.702 </td><td>-21.224</td><td>-25.741</td><td> 2.888 </td><td> -8.438</td><td>-17.055</td><td>-12.026</td><td> -9.667</td><td>-20.814</td><td>  2.463</td></tr>
</tbody>
</table>



write.csv(TS7,
            file = "~/Desktop/SuppTable7_CRISPR_rawBFs.csv", 
            col.names = T,
            row.names = T
           )

----
# 4.0 Compare Essential Genes / sample with MacLeod et al. 2019
----
Download Norm BFs from MacLeod et al., 2019 Table S2

- 1) Bar plot of number of essential genes / sample in ours vs MacLeod et al.
- 2) Overlap of essential genes in our cohort vs MacLeod


```R
### load in Grahams paper

published <- read.csv("~/Downloads/mmc3 (4).csv")
head(published)

### load in house data

setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/")

inhouse <- readRDS("CRISPR_Avg_Diff_Matrix.rds")
head(inhouse)

```


<table>
<thead><tr><th scope=col>X</th><th scope=col>G549NS</th><th scope=col>G583NS</th><th scope=col>G432NS</th><th scope=col>G440NS</th><th scope=col>G472NS</th><th scope=col>G477NS</th><th scope=col>G510NS</th><th scope=col>G361NS</th><th scope=col>HF7450</th><th scope=col>⋯</th><th scope=col>HAP1</th><th scope=col>HCT116</th><th scope=col>HELA</th><th scope=col>HPAFii</th><th scope=col>PATU8988s</th><th scope=col>RPE1</th><th scope=col>SUM149PT</th><th scope=col>DLD1</th><th scope=col>Z.score.GBM.vs..NSC</th><th scope=col>Z.score.GBM.vs..other.ex.NSC</th></tr></thead>
<tbody>
	<tr><td>JUN       </td><td> 97.47800 </td><td> 42.044667</td><td> 11.78667 </td><td>103.160556</td><td>  1.817111</td><td>50.94022  </td><td> 88.054667</td><td>-33.29122 </td><td>-68.575500</td><td>⋯         </td><td>-48.39150 </td><td>-41.470389</td><td>-64.07833 </td><td>-50.006667</td><td> -94.45667</td><td>-49.51211 </td><td>-26.58050 </td><td> -54.27397</td><td> 5.0850316</td><td>8.717187  </td></tr>
	<tr><td>SQLE      </td><td> 40.54856 </td><td>  7.615778</td><td> 40.00150 </td><td> -1.862833</td><td> 50.716444</td><td>71.12239  </td><td> 87.796778</td><td> 21.56833 </td><td> -9.068111</td><td>⋯         </td><td>-46.61517 </td><td>-51.270444</td><td>-37.67761 </td><td>-27.136278</td><td>  13.98039</td><td>-69.50706 </td><td>-30.12206 </td><td> -49.79594</td><td> 0.2994659</td><td>6.786813  </td></tr>
	<tr><td>SOX9      </td><td> 14.74061 </td><td> 40.486500</td><td> 50.22339 </td><td> 19.915278</td><td> 80.265667</td><td>94.42917  </td><td>-52.826722</td><td> 47.17117 </td><td>-25.908944</td><td>⋯         </td><td>-57.52233 </td><td> -2.514278</td><td>-39.75783 </td><td> 54.981389</td><td>-142.63917</td><td>-34.56822 </td><td>-53.57028 </td><td> -22.32247</td><td> 2.8833023</td><td>6.529596  </td></tr>
	<tr><td>SOX2      </td><td>-13.38967 </td><td>-33.731694</td><td> 28.86683 </td><td> 52.199833</td><td> 12.679611</td><td>32.71478  </td><td>  6.038389</td><td>-30.29633 </td><td> 61.616111</td><td>⋯         </td><td>-21.86828 </td><td>-46.147778</td><td>-49.18978 </td><td>-37.887000</td><td> -95.03006</td><td>-30.92989 </td><td>-46.55550 </td><td>-105.81556</td><td>-0.1768105</td><td>5.385187  </td></tr>
	<tr><td>CDK6      </td><td> 96.12300 </td><td>-39.797639</td><td>-29.69244 </td><td> 63.371778</td><td>106.376444</td><td>73.68922  </td><td> 65.610556</td><td> 92.00217 </td><td> 88.862500</td><td>⋯         </td><td>-18.70139 </td><td> 17.623500</td><td>-28.29022 </td><td>  9.250722</td><td>  29.10783</td><td> 23.36900 </td><td>-33.08367 </td><td> -45.63611</td><td>-1.7524959</td><td>5.225759  </td></tr>
	<tr><td>SOCS3     </td><td> 42.76294 </td><td> -9.191500</td><td> 48.31506 </td><td>-23.201667</td><td> 97.873389</td><td>54.55506  </td><td> 28.427556</td><td>-22.75894 </td><td> -4.405389</td><td>⋯         </td><td>-39.86256 </td><td>-30.093889</td><td>-31.54206 </td><td>  1.001722</td><td> -32.65178</td><td>-38.05142 </td><td>-46.71478 </td><td> -23.58161</td><td> 0.5142259</td><td>5.051967  </td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th><th scope=col>Developmental_Average</th><th scope=col>InjuryResponse_Average</th><th scope=col>Diff</th><th scope=col>Diff_Z</th><th scope=col>FoldChange</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>40.43901   </td><td>40.93868   </td><td>-16.8872393</td><td>30.37648   </td><td>47.26371   </td><td>7.051752   </td><td> -1.798783 </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td>33.67403   </td><td>37.50733   </td><td> -8.8287541</td><td>32.77162   </td><td>41.60037   </td><td>6.206774   </td><td> -3.711919 </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td>33.63656   </td><td>47.18051   </td><td> -0.7493837</td><td>39.18887   </td><td>39.93825   </td><td>5.958783   </td><td>-52.294793 </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>37.95362   </td><td>53.55119   </td><td>  4.2116815</td><td>43.98490   </td><td>39.77322   </td><td>5.934161   </td><td> 10.443548 </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td>14.84797   </td><td>28.55795   </td><td> -6.9819679</td><td>37.50056   </td><td>44.48253   </td><td>6.636796   </td><td> -5.371059 </td></tr>
	<tr><th scope=row>AARS</th><td>  0.9270344</td><td>  4.774409 </td><td>-10.0122166</td><td> 34.759483 </td><td>52.83310   </td><td>47.83145   </td><td>37.41367   </td><td>39.54539   </td><td>28.55795   </td><td>47.44578   </td><td>46.42611   </td><td>  7.6121776</td><td>41.23631   </td><td>33.62414   </td><td>5.016710   </td><td>  5.417151 </td></tr>
</tbody>
</table>




```R
inhouse <- inhouse[ ,c(1:11)]
head(inhouse)

rownames(published) <- published$X
published$X <- NULL
published <- published[ ,c(1:8)]
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th></tr></thead>
<tbody>
	<tr><th scope=row>ITGB1</th><td>-28.2717311</td><td>-15.353046 </td><td>-19.6011966</td><td> -4.322983 </td><td>27.04493   </td><td>35.34420   </td><td>29.13527   </td><td>29.02607   </td><td>31.33191   </td><td>40.43901   </td><td>40.93868   </td></tr>
	<tr><th scope=row>SCAP</th><td> 32.5498409</td><td>-14.911093 </td><td>-28.5620940</td><td>-24.391670 </td><td>38.35182   </td><td>46.03781   </td><td>42.90056   </td><td>19.55718   </td><td>17.01073   </td><td>33.67403   </td><td>37.50733   </td></tr>
	<tr><th scope=row>EXOSC5</th><td> -5.7701043</td><td>-10.239477 </td><td> -9.2862183</td><td> 22.298264 </td><td>45.94805   </td><td>30.57221   </td><td>41.65691   </td><td>46.25838   </td><td>31.50877   </td><td>33.63656   </td><td>47.18051   </td></tr>
	<tr><th scope=row>SUPT16H</th><td> -2.8559548</td><td> -1.281364 </td><td> -0.9701818</td><td> 21.954226 </td><td>39.91896   </td><td>42.16704   </td><td>46.14875   </td><td>43.38775   </td><td>48.30201   </td><td>37.95362   </td><td>53.55119   </td></tr>
	<tr><th scope=row>PPP1R12A</th><td>  0.1890148</td><td> -4.982018 </td><td>-16.6857273</td><td> -6.449141 </td><td>35.65962   </td><td>55.71300   </td><td>34.50681   </td><td>19.10447   </td><td>42.51892   </td><td>14.84797   </td><td>28.55795   </td></tr>
	<tr><th scope=row>AARS</th><td>  0.9270344</td><td>  4.774409 </td><td>-10.0122166</td><td> 34.759483 </td><td>52.83310   </td><td>47.83145   </td><td>37.41367   </td><td>39.54539   </td><td>28.55795   </td><td>47.44578   </td><td>46.42611   </td></tr>
</tbody>
</table>




```R
common.inhouse <- inhouse[rowSums(inhouse >= 10) == 11, ]
head(common.inhouse)
dim(common.inhouse)
```


<table>
<thead><tr><th></th><th scope=col>BT67_L</th><th scope=col>G620_L</th><th scope=col>G523_L</th><th scope=col>G809_L</th><th scope=col>G583_L</th><th scope=col>G549_L</th><th scope=col>G729_L</th><th scope=col>G361_L</th><th scope=col>G691_L</th><th scope=col>G440_L</th><th scope=col>G532_L</th></tr></thead>
<tbody>
	<tr><th scope=row>SAMM50</th><td>28.00084</td><td>16.82670</td><td>11.95773</td><td>13.83357</td><td>26.27241</td><td>48.05793</td><td>38.15273</td><td>40.51389</td><td>51.87056</td><td>41.08864</td><td>37.95362</td></tr>
	<tr><th scope=row>CCT3</th><td>18.24408</td><td>30.40479</td><td>32.93801</td><td>23.21481</td><td>35.15263</td><td>45.58241</td><td>57.04537</td><td>47.09073</td><td>53.55119</td><td>54.94900</td><td>49.20811</td></tr>
	<tr><th scope=row>RPL14</th><td>29.46006</td><td>20.09341</td><td>22.16101</td><td>30.09420</td><td>40.36555</td><td>49.31820</td><td>44.69244</td><td>52.24693</td><td>38.70517</td><td>50.03037</td><td>54.94900</td></tr>
	<tr><th scope=row>PSMA4</th><td>28.88840</td><td>31.44402</td><td>22.29826</td><td>19.30917</td><td>46.58458</td><td>46.64892</td><td>44.90366</td><td>49.65728</td><td>29.83860</td><td>43.15864</td><td>55.43283</td></tr>
	<tr><th scope=row>RUVBL1</th><td>37.18315</td><td>13.81443</td><td>48.90119</td><td>20.65362</td><td>44.76249</td><td>36.11393</td><td>54.94900</td><td>55.71300</td><td>49.00430</td><td>50.13046</td><td>50.89800</td></tr>
	<tr><th scope=row>DDX46</th><td>30.22008</td><td>34.67258</td><td>34.24049</td><td>23.46250</td><td>59.82782</td><td>41.90282</td><td>48.16158</td><td>54.03386</td><td>49.74000</td><td>44.69244</td><td>46.70927</td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>234</li>
	<li>11</li>
</ol>




```R
common.published <- published[rowSums(published >= 10) == 8, ]
head(common.published)
dim(common.published)
```


<table>
<thead><tr><th></th><th scope=col>G549NS</th><th scope=col>G583NS</th><th scope=col>G432NS</th><th scope=col>G440NS</th><th scope=col>G472NS</th><th scope=col>G477NS</th><th scope=col>G510NS</th><th scope=col>G361NS</th></tr></thead>
<tbody>
	<tr><th scope=row>CENPW</th><td>89.26472 </td><td> 82.88722</td><td>89.95572 </td><td>71.39600 </td><td>103.58300</td><td>77.48883 </td><td>126.75622</td><td>59.59550 </td></tr>
	<tr><th scope=row>XRN1</th><td>10.16317 </td><td> 63.92772</td><td>63.48311 </td><td>98.34533 </td><td> 53.86806</td><td>52.63394 </td><td> 59.51606</td><td>24.85633 </td></tr>
	<tr><th scope=row>TRAIP</th><td>64.12600 </td><td>106.37644</td><td>56.24794 </td><td>40.28667 </td><td> 55.05606</td><td>54.15322 </td><td> 47.10150</td><td>39.73922 </td></tr>
	<tr><th scope=row>CHMP4B</th><td>24.96000 </td><td> 14.10011</td><td>37.50794 </td><td>54.72017 </td><td> 24.81394</td><td>63.92772 </td><td> 72.11039</td><td>84.49544 </td></tr>
	<tr><th scope=row>NOP14</th><td>32.35411 </td><td> 44.65794</td><td>75.70100 </td><td>59.66339 </td><td> 45.29983</td><td>56.12433 </td><td> 42.76294</td><td>44.42550 </td></tr>
	<tr><th scope=row>SMC4</th><td>35.14761 </td><td> 66.96567</td><td>75.56267 </td><td>88.86250 </td><td> 56.82439</td><td>66.76494 </td><td> 55.35189</td><td>29.49583 </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>262</li>
	<li>8</li>
</ol>


