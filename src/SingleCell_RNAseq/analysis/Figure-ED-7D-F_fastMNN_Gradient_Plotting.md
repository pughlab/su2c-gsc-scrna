
---
# Plot the fastMNN gradient
---


```R
library(ggplot2)
library(ggpubr)
setwd("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/fastMNN/")
```


```R
### load data
dat <- readRDS("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/fastMNN/fastMNN_DevIR_AUCell_GSCs.rds")
```

### Plot gradient


```R
IR <-  dat[ ,c("InjuryResponse_GSC_AUC", "SampleID")]
colnames(IR)[1] <- "Score"
IR$signature <- "InjuryResponse"
IR$Score <- scale(IR$Score, center = TRUE, scale = TRUE)
IR$SampleID <- with(IR, reorder(SampleID, Score, median))



Dev <- dat[ ,c("Developmental_GSC_AUC", "SampleID")]
colnames(Dev)[1] <- "Score"
Dev$signature <- "Developmental"
Dev$Score <- scale(Dev$Score, center = TRUE, scale = TRUE)
Dev$SampleID <- with(IR, reorder(SampleID, Score, median))




dat <- rbind(Dev, IR)
head(dat)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Score</th><th scope=col>SampleID</th><th scope=col>signature</th></tr>
	<tr><th></th><th scope=col>&lt;dbl[,1]&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td> 2.7477693</td><td>BT127_L</td><td>Developmental</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>-0.4377524</td><td>BT127_L</td><td>Developmental</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td> 0.8535018</td><td>BT127_L</td><td>Developmental</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td> 0.0631048</td><td>BT127_L</td><td>Developmental</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>-1.3552608</td><td>BT127_L</td><td>Developmental</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td> 2.5307117</td><td>BT127_L</td><td>Developmental</td></tr>
</tbody>
</table>




```R
pdf("fastMNN_Corrected_Gradient.pdf", width = 10, height = 3)

p <- ggplot(dat, aes(x=SampleID, y=Score, fill = signature)) +theme_classic() +
      geom_violin()  + ylab("Relative AUC Score (z-score)") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
scale_fill_manual(values=c("red", "black")) + scale_x_discrete()  
p

dev.off()
```


<strong>pdf:</strong> 2


### Correlate medians between methods


```R
### load dat for original
dat <- readRDS("fastMNN_DevIR_AUCell_GSCs.rds")
orig <- readRDS("Original_DevIR_AUCell_GSCs.rds")
dat <- cbind(dat, orig)
head(dat)
```


<table>
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>Developmental_GSC_AUC</th><th scope=col>InjuryResponse_GSC_AUC</th><th scope=col>SampleID</th><th scope=col>Original_Developmental_GSC_AUC</th><th scope=col>Original_InjuryResponse_GSC_AUC</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>BT127_L_AAACCTGCACGGACAA</th><td>0.3543462</td><td>0.08754009</td><td>BT127_L</td><td>0.1877557</td><td>0.1704331</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGCATCCGGGT</th><td>0.2138513</td><td>0.14121799</td><td>BT127_L</td><td>0.1262983</td><td>0.1621965</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGGTACAGTTC</th><td>0.2708010</td><td>0.16961135</td><td>BT127_L</td><td>0.1389073</td><td>0.1508963</td></tr>
	<tr><th scope=row>BT127_L_AAACCTGTCTACGAGT</th><td>0.2359412</td><td>0.19577983</td><td>BT127_L</td><td>0.1516215</td><td>0.1646188</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGAGTGGTAAT</th><td>0.1733853</td><td>0.16569369</td><td>BT127_L</td><td>0.1326554</td><td>0.1643399</td></tr>
	<tr><th scope=row>BT127_L_AAACGGGCAGGACGTA</th><td>0.3447731</td><td>0.14392559</td><td>BT127_L</td><td>0.1626407</td><td>0.2095603</td></tr>
</tbody>
</table>




```R
Dev <- aggregate(Developmental_GSC_AUC~SampleID, dat ,median)
Dev$original <- aggregate(Original_Developmental_GSC_AUC~SampleID, dat ,median)$Original_Developmental_GSC_AUC
head(Dev)

IR <- aggregate(InjuryResponse_GSC_AUC~SampleID, dat ,median)
IR$original <- aggregate(Original_InjuryResponse_GSC_AUC~SampleID, dat ,median)$Original_InjuryResponse_GSC_AUC
head(IR)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>SampleID</th><th scope=col>Developmental_GSC_AUC</th><th scope=col>original</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>BT127_L</td><td>0.2609854</td><td>0.1566859</td></tr>
	<tr><th scope=row>2</th><td>BT147_L</td><td>0.2332861</td><td>0.1449265</td></tr>
	<tr><th scope=row>3</th><td>BT48_L </td><td>0.2287109</td><td>0.1686587</td></tr>
	<tr><th scope=row>4</th><td>BT67_L </td><td>0.2404907</td><td>0.1513844</td></tr>
	<tr><th scope=row>5</th><td>BT73_L </td><td>0.2727298</td><td>0.1165985</td></tr>
	<tr><th scope=row>6</th><td>BT84_L </td><td>0.2913912</td><td>0.1625484</td></tr>
</tbody>
</table>




<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>SampleID</th><th scope=col>InjuryResponse_GSC_AUC</th><th scope=col>original</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>BT127_L</td><td>0.1809427</td><td>0.1789142</td></tr>
	<tr><th scope=row>2</th><td>BT147_L</td><td>0.1584216</td><td>0.1411933</td></tr>
	<tr><th scope=row>3</th><td>BT48_L </td><td>0.1620228</td><td>0.1393815</td></tr>
	<tr><th scope=row>4</th><td>BT67_L </td><td>0.1613491</td><td>0.1494588</td></tr>
	<tr><th scope=row>5</th><td>BT73_L </td><td>0.1474356</td><td>0.1600811</td></tr>
	<tr><th scope=row>6</th><td>BT84_L </td><td>0.1580783</td><td>0.1555337</td></tr>
</tbody>
</table>




```R


a <- ggscatter(Dev,
          x= "original",
          y = "Developmental_GSC_AUC",
          xlab = "Original (Median Raw AUC Score)",
          ylab = "Corrected fastMNN Matrix (Median Raw AUC Score)",
          main = "Developmental",
          shape = 21,
          size = 3,
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson",
                                #label.x = 0,
                                label.sep = "\n")
          )

cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")


b <- ggscatter(IR,
          x= "original",
          y = "InjuryResponse_GSC_AUC",
          xlab = "Original (Median Raw AUC Score)",
          ylab = "Corrected fastMNN Matrix (Median Raw AUC Score)",
          main = "Injury Response",
          shape = 21,
          size = 3,
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson",
                                #label.x = 0,
                                label.sep = "\n")
          )

cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
```


```R
pdf("MedianCorrelation_fastMNN_Orig.pdf", width = 10, height = 5 )
ggarrange(a, b)
dev.off()
```

    `geom_smooth()` using formula 'y ~ x'

    `geom_smooth()` using formula 'y ~ x'




<strong>pdf:</strong> 2



```R

```
