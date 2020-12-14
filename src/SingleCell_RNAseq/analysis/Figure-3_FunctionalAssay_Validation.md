
---
# Investigate functional differences between GSC states
---
L.Richards

/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/FunctionalAssay_Validation


```R
options(repos='http://cran.rstudio.com/')
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(survminer) ##v0.4.8
library(survival) ##v3.2-3
```


```R
### set working dir and read in data
setwd("~/Desktop/H4H/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/FunctionalAssay_Validation")
dat <- read.csv("SU2C_scRNA_FunctionalAssays_Aug2020.csv")
colnames(dat)[8] <- "SPHERE_FORMATION_CAPACITY"
```

----
## 1.0 Survival Analysis
---
https://www.datacamp.com/community/tutorials/survival-analysis-R  
http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization  
https://rstudio.com/wp-content/uploads/2015/01/survminer-1.png

---
### 1.1 All samples, bulk RNAseq classification
---

Survival analysis of all lines with xenograft information and bin with bulk RNAsequencing classification
Take median xeno infromation for mice


```R
##subset data to only those with bulkRNA classification and xenograft
sub <- dat[!is.na(dat$BULK_RNASEQ_CLASSIFICATION), ]
remove <- c("N.A.", "N.D.", "NoTumours", "Underway")
sub <- sub[!sub$MEDIAN_ORTHOTOPIC_XENOGRAFT_SURVIVAL_DAYS %in% remove, ]
sub$STATUS <- 1
sub$MEDIAN_ORTHOTOPIC_XENOGRAFT_SURVIVAL_DAYS <- as.numeric(as.character(sub$MEDIAN_ORTHOTOPIC_XENOGRAFT_SURVIVAL_DAYS))
table(sub$BULK_RNASEQ_CLASSIFICATION)
```



     Developmental InjuryResponse
                23             14



```R
surv_object <- Surv(time = sub$MEDIAN_ORTHOTOPIC_XENOGRAFT_SURVIVAL_DAYS,
                    event = sub$STATUS
                   )
```


```R
fit1 <- survfit(surv_object ~ BULK_RNASEQ_CLASSIFICATION,
                data = sub
               )
#summary(fit1)

```


```R
pdf("Dev_IR_AllSamples_KM.pdf", width = 4, height = 6)
ggsurvplot(fit1,
           data = sub,
           pval = TRUE,
           legend.labs = c("Developmental", "InjuryResponse"),
           palette = c("red", "black"),
           risk.table = T,
           xlab = "Time (days)",
           ggtheme = theme_bw()
          )
```

    Warning message:
    “Vectorized input to `element_text()` is not officially supported.
    Results may be unexpected or may change in future versions of ggplot2.”



![png](output_8_1.png)


---
## 2.0 Compare SFC between Dev and IR bulk RNA classifications
---


```R
remove <- c("N.A.", "N.D.", "NoTumours", "Underway") #remove incomplete values
sub <- dat[!dat$SPHERE_FORMATION_CAPACITY %in% remove, ]
sub$SPHERE_FORMATION_CAPACITY <- as.numeric(as.character(sub$SPHERE_FORMATION_CAPACITY))
sub <- sub[!is.na(sub$BULK_RNASEQ_CLASSIFICATION), ]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>60</li><li>16</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>55</li><li>16</li></ol>




```R
table(sub$BULK_RNASEQ_CLASSIFICATION)
```



     Developmental InjuryResponse
                30             25



```R
### remove outlier in Developmetnal group
my_comparisons <- list( c("Developmental", "InjuryResponse")
                      )
p <- ggboxplot(sub2, x = "BULK_RNASEQ_CLASSIFICATION", y = "SPHERE_FORMATION_CAPACITY",
          color = "BULK_RNASEQ_CLASSIFICATION", palette = c("black", "darkblue", "red"),
          add = "jitter") + stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons p-value
p
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>54</li><li>16</li></ol>




![png](output_12_1.png)



```R
### test significance of tumour formation between groups
fisher.test(matrix(c(23,0,11,3), nrow = 2, ncol = 2))
```



    	Fisher's Exact Test for Count Data

    data:  matrix(c(23, 0, 11, 3), nrow = 2, ncol = 2)
    p-value = 0.04685
    alternative hypothesis: true odds ratio is not equal to 1
    95 percent confidence interval:
     0.7268299       Inf
    sample estimates:
    odds ratio
           Inf
