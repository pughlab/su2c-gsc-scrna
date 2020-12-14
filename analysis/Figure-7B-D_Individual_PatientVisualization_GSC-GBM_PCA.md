
----
# Visualize each patient
---

L.Richards


```R
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(RColorBrewer)
library(gridExtra)
library(ks)
```

    Warning message:
    “package ‘ggplot2’ was built under R version 3.4.4”Warning message:
    “package ‘ggpubr’ was built under R version 3.4.4”Loading required package: magrittr
    Warning message:
    “package ‘ggExtra’ was built under R version 3.4.4”

### Read in and format data


```R
pc.file <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/BTSC_TumourCell_G800Removed_PCA_allGenes.Rdata"
meta.file <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/BTSC_TumourCell_G800Removed_PCA_metdata.Rdata"
aucell <- "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/G800_Removed/AUCell/BTSC_Tumour_G800Removed_AUCellScores.Rdata"

load(pc.file)  
load(meta.file)
load(aucell)

file.prefix <- "G800L_Removed"

a <- strsplit(rownames(BTSC_Tumour_AUCell), "_")
a <- matrix(unlist(a), ncol = 4, byrow = TRUE)
#head(a)

BTSC_Tumour_AUCell$SampleID <- paste(a[,2], a[,3], sep = "_")
BTSC_Tumour_AUCell$SampleType <- ifelse(grepl("_L", rownames(BTSC_Tumour_AUCell)), "GSC", "Tumour")
#head(BTSC_Tumour_AUCell)

dff <- data.frame(pc@cell.embeddings)
dff <- cbind(dff, BTSC_Tumour_AUCell)

#### filter sigs to be only ones we care about

keep <- c("PC1", "PC2", "SampleType", "SampleID", "RNA.GSC.c1_AUC", "RNA.GSC.c2_AUC")
df <- dff[ ,keep]
head(df)

```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>SampleType</th><th scope=col>SampleID</th><th scope=col>RNA.GSC.c1_AUC</th><th scope=col>RNA.GSC.c2_AUC</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>-12.7289371</td><td>-0.5538633 </td><td>GSC        </td><td>BT127_L    </td><td>0.1833756  </td><td>0.1691589  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>  0.9425536</td><td> 1.6079858 </td><td>GSC        </td><td>BT127_L    </td><td>0.1225561  </td><td>0.1561081  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>  4.2237210</td><td>17.3566342 </td><td>GSC        </td><td>BT127_L    </td><td>0.1403795  </td><td>0.1532729  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>-14.3646978</td><td>12.1673063 </td><td>GSC        </td><td>BT127_L    </td><td>0.1531699  </td><td>0.1659836  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>  7.0178336</td><td>23.3358358 </td><td>GSC        </td><td>BT127_L    </td><td>0.1330963  </td><td>0.1653948  </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>-14.7980172</td><td>-8.7067006 </td><td>GSC        </td><td>BT127_L    </td><td>0.1597614  </td><td>0.2143644  </td></tr>
</tbody>
</table>




```R
#### calculate C1-C2 difference

df$Pro_Immes_Diff <- df$RNA.GSC.c1_AUC - df$RNA.GSC.c2_AUC

### normalize the scores between 0-1
### then subtract the normalized scores

df$RNA_C1_norm <- (df$RNA.GSC.c1_AUC - min(df$RNA.GSC.c1_AUC)) / (max(df$RNA.GSC.c1_AUC) - min(df$RNA.GSC.c1_AUC))
df$RNA_C2_norm <- (df$RNA.GSC.c2_AUC - min(df$RNA.GSC.c2_AUC)) / (max(df$RNA.GSC.c2_AUC) - min(df$RNA.GSC.c2_AUC))
df$Diff <- df$RNA_C1_norm - df$RNA_C2_norm

head(df)
```


<table>
<thead><tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>SampleType</th><th scope=col>SampleID</th><th scope=col>RNA.GSC.c1_AUC</th><th scope=col>RNA.GSC.c2_AUC</th><th scope=col>Pro_Immes_Diff</th><th scope=col>RNA_C1_norm</th><th scope=col>RNA_C2_norm</th><th scope=col>Diff</th></tr></thead>
<tbody>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCACGGACAA</th><td>-12.7289371</td><td>-0.5538633 </td><td>GSC        </td><td>BT127_L    </td><td>0.1833756  </td><td>0.1691589  </td><td> 0.01421673</td><td>0.5502658  </td><td>0.2204938  </td><td>0.32977199 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGCATCCGGGT</th><td>  0.9425536</td><td> 1.6079858 </td><td>GSC        </td><td>BT127_L    </td><td>0.1225561  </td><td>0.1561081  </td><td>-0.03355195</td><td>0.3172655  </td><td>0.1786540  </td><td>0.13861148 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGGTACAGTTC</th><td>  4.2237210</td><td>17.3566342 </td><td>GSC        </td><td>BT127_L    </td><td>0.1403795  </td><td>0.1532729  </td><td>-0.01289338</td><td>0.3855471  </td><td>0.1695646  </td><td>0.21598251 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACCTGTCTACGAGT</th><td>-14.3646978</td><td>12.1673063 </td><td>GSC        </td><td>BT127_L    </td><td>0.1531699  </td><td>0.1659836  </td><td>-0.01281371</td><td>0.4345473  </td><td>0.2103142  </td><td>0.22423312 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGAGTGGTAAT</th><td>  7.0178336</td><td>23.3358358 </td><td>GSC        </td><td>BT127_L    </td><td>0.1330963  </td><td>0.1653948  </td><td>-0.03229856</td><td>0.3576448  </td><td>0.2084265  </td><td>0.14921831 </td></tr>
	<tr><th scope=row>BTSC_BT127_L_AAACGGGCAGGACGTA</th><td>-14.7980172</td><td>-8.7067006 </td><td>GSC        </td><td>BT127_L    </td><td>0.1597614  </td><td>0.2143644  </td><td>-0.05460305</td><td>0.4597993  </td><td>0.3654192  </td><td>0.09438006 </td></tr>
</tbody>
</table>



### Calculate Contours


```R
###Calcualte contours

set.seed(1001)

#define 95% and 50% contour levels for GSCs
GSC <- data.frame(x=df[df$SampleType == "GSC", "PC1"], 
                  y=df[df$SampleType == "GSC", "PC2"]
                 )

kd <- ks::kde(GSC, compute.cont=TRUE)
GSC_contour_95 <- with(kd, contourLines(x=eval.points[[1]], 
                                        y=eval.points[[2]],
                                        z=estimate, 
                                        levels=cont["1%"])[[1]]
                      )
GSC_contour_95 <- data.frame(GSC_contour_95)

GSC_contour_50 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
GSC_contour_50 <- data.frame(GSC_contour_50)


###### Tumours

set.seed(1001)

#define 95% and 50% contour levels for Tumour
Tumour <- data.frame(x=df[df$SampleType == "Tumour", "PC1"], 
                  y=df[df$SampleType == "Tumour", "PC2"]
                 )

kd <- ks::kde(Tumour, compute.cont=TRUE)
Tumour_contour_95 <- with(kd, contourLines(x=eval.points[[1]], 
                                        y=eval.points[[2]],
                                        z=estimate, 
                                        levels=cont["1%"])[[1]]
                      )
Tumour_contour_95 <- data.frame(Tumour_contour_95)

Tumour_contour_50 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["5%"])[[1]])
#Tumour_contour_50 <- data.frame(Tumour_contour_50)
```

### Plot data


```R
### write function to plot sample counts

plot_sample_hex <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "black")

    
### plot patient counts
p <- ggplot(subset, aes_string(x=x, y=y, z="SampleID")) +
 geom_hex(bins=100) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) 


 p + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    }

#plot score


plot_sample_hex_AUC <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "red")

q <- ggplot(subset, aes_string(x=x, y=y, z="Diff")) +
 stat_summary_hex(bins=100, fun = "median") +
#scale_fill_gradientn("Median Raw \nAUC Score", 
#                     colours = c(rev(brewer.pal(n = 8, name = "RdYlBu"))), 
#                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
#                    ) +
    
scale_fill_gradient2("Median Raw \nAUC Score",
                     low = "blue",
                     mid = "white",
                     high = "red",
                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
                  ) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) 

 q + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    
}




```


```R
sigs <- sort(unique(df$SampleID)) #in this case sigs is patient 
sigs
```


<ol class=list-inline>
	<li>'BT127_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
	<li>'BT67_L'</li>
	<li>'BT73_L'</li>
	<li>'BT84_L'</li>
	<li>'BT89_L'</li>
	<li>'BT94_L'</li>
	<li>'G1003-A_T'</li>
	<li>'G1003-B_T'</li>
	<li>'G1003-C_T'</li>
	<li>'G1003-D_T'</li>
	<li>'G523_L'</li>
	<li>'G549_L'</li>
	<li>'G564_L'</li>
	<li>'G566_L'</li>
	<li>'G583_L'</li>
	<li>'G620_L'</li>
	<li>'G620_T'</li>
	<li>'G637_L'</li>
	<li>'G729_L'</li>
	<li>'G797_L'</li>
	<li>'G799_L'</li>
	<li>'G837_L'</li>
	<li>'G851_L'</li>
	<li>'G876_L'</li>
	<li>'G885_L'</li>
	<li>'G895_L'</li>
	<li>'G910-A_T'</li>
	<li>'G910-B_T'</li>
	<li>'G910-C_T'</li>
	<li>'G910-D_T'</li>
	<li>'G910-E_T'</li>
	<li>'G945-I_L'</li>
	<li>'G945-I_T'</li>
	<li>'G945-J_L'</li>
	<li>'G945-J_T'</li>
	<li>'G945-K_L'</li>
	<li>'G945-K_T'</li>
	<li>'G946-I_T'</li>
	<li>'G946-J_L'</li>
	<li>'G946-J_T'</li>
	<li>'G946-K_L'</li>
	<li>'G946-K_T'</li>
	<li>'G967-A_T'</li>
	<li>'G967-B_T'</li>
	<li>'G967-C_T'</li>
	<li>'G967-D_T'</li>
	<li>'G983-A_T'</li>
	<li>'G983-B_T'</li>
	<li>'G983-C_T'</li>
</ol>




```R
sigs <- sort(unique(df$SampleID)) #in this case sigs is patient 
sigs

p <- list()
q <- list()

for(i in 1:length(sigs)){

sig <- sigs[i]
print(sig)
subset <- df[df$SampleID == sigs[i], ]
p[[i]] <- plot_sample_hex(subset, "PC1", "PC2", sigs[i])
q[[i]] <- plot_sample_hex_AUC(subset, "PC1", "PC2", sigs[i])
                                        
}

```


<ol class=list-inline>
	<li>'BT127_L'</li>
	<li>'BT147_L'</li>
	<li>'BT48_L'</li>
	<li>'BT67_L'</li>
	<li>'BT73_L'</li>
	<li>'BT84_L'</li>
	<li>'BT89_L'</li>
	<li>'BT94_L'</li>
	<li>'G1003-A_T'</li>
	<li>'G1003-B_T'</li>
	<li>'G1003-C_T'</li>
	<li>'G1003-D_T'</li>
	<li>'G523_L'</li>
	<li>'G549_L'</li>
	<li>'G564_L'</li>
	<li>'G566_L'</li>
	<li>'G583_L'</li>
	<li>'G620_L'</li>
	<li>'G620_T'</li>
	<li>'G637_L'</li>
	<li>'G729_L'</li>
	<li>'G797_L'</li>
	<li>'G799_L'</li>
	<li>'G837_L'</li>
	<li>'G851_L'</li>
	<li>'G876_L'</li>
	<li>'G885_L'</li>
	<li>'G895_L'</li>
	<li>'G910-A_T'</li>
	<li>'G910-B_T'</li>
	<li>'G910-C_T'</li>
	<li>'G910-D_T'</li>
	<li>'G910-E_T'</li>
	<li>'G945-I_L'</li>
	<li>'G945-I_T'</li>
	<li>'G945-J_L'</li>
	<li>'G945-J_T'</li>
	<li>'G945-K_L'</li>
	<li>'G945-K_T'</li>
	<li>'G946-I_T'</li>
	<li>'G946-J_L'</li>
	<li>'G946-J_T'</li>
	<li>'G946-K_L'</li>
	<li>'G946-K_T'</li>
	<li>'G967-A_T'</li>
	<li>'G967-B_T'</li>
	<li>'G967-C_T'</li>
	<li>'G967-D_T'</li>
	<li>'G983-A_T'</li>
	<li>'G983-B_T'</li>
	<li>'G983-C_T'</li>
</ol>



    [1] "BT127_L"
    [1] "BT147_L"
    [1] "BT48_L"
    [1] "BT67_L"
    [1] "BT73_L"
    [1] "BT84_L"
    [1] "BT89_L"
    [1] "BT94_L"
    [1] "G1003-A_T"
    [1] "G1003-B_T"
    [1] "G1003-C_T"
    [1] "G1003-D_T"
    [1] "G523_L"
    [1] "G549_L"
    [1] "G564_L"
    [1] "G566_L"
    [1] "G583_L"
    [1] "G620_L"
    [1] "G620_T"
    [1] "G637_L"
    [1] "G729_L"
    [1] "G797_L"
    [1] "G799_L"
    [1] "G837_L"
    [1] "G851_L"
    [1] "G876_L"
    [1] "G885_L"
    [1] "G895_L"
    [1] "G910-A_T"
    [1] "G910-B_T"
    [1] "G910-C_T"
    [1] "G910-D_T"
    [1] "G910-E_T"
    [1] "G945-I_L"
    [1] "G945-I_T"
    [1] "G945-J_L"
    [1] "G945-J_T"
    [1] "G945-K_L"
    [1] "G945-K_T"
    [1] "G946-I_T"
    [1] "G946-J_L"
    [1] "G946-J_T"
    [1] "G946-K_L"
    [1] "G946-K_T"
    [1] "G967-A_T"
    [1] "G967-B_T"
    [1] "G967-C_T"
    [1] "G967-D_T"
    [1] "G983-A_T"
    [1] "G983-B_T"
    [1] "G983-C_T"



```R
#ml <- marrangeGrob(q, nrow=3, ncol=3)
#ggsave("~/Desktop/GSCs_Tumours_patient_positions_score.pdf", ml, height = 15, width = 15, dpi = 72, limitsize = FALSE)
```


```R
pdf("~/Desktop/Gradient_patient.pdf", width = 11, height = 5 )

for(i in 1:length(q)){
    
    plist <- list(p[[i]], q[[i]])
    grid.arrange(grobs = plist, ncol = 2)
    
}

dev.off()
```

    Warning message:
    “package ‘hexbin’ was built under R version 3.4.3”


<strong>pdf:</strong> 2



```R
ml <- marrangeGrob(q, nrow=3, ncol=3)
ggsave("~/Desktop/GSCs_Tumours_patient_positions_score_March2020_BlRd.pdf", 
       ml, 
       height = 12, 
       width = 15, 
       dpi = 72, 
       limitsize = FALSE
      )
```

----
## Plot pairs together
---


```R
### write function to plot sample counts

plot_sample_hex <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "black")

    
### plot patient counts
p <- ggplot(subset, aes_string(x=x, y=y, z="SampleID")) +
 geom_hex(bins=100) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) 


 p + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    }

#plot score


plot_sample_hex_AUC <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "black")

q <- ggplot(subset, aes_string(x=x, y=y, z="Diff")) +
 stat_summary_hex(bins=100, fun = "median") +
#scale_fill_gradientn("Median Raw \nAUC Score", 
#                     colours = c(rev(brewer.pal(n = 8, name = "RdYlBu"))), 
#                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
#                    ) +
    
scale_fill_gradient2("Median Raw \nAUC Score",
                     low = "blue",
                     mid = "white",
                     high = "red",
                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
                  ) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) 

 q + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) 
    
    
}




```


```R
sigs <- sort(unique(df$SampleID)) #in this case sigs is patient 

sigs <- c("G620_L", "G620_T", "G945-I_L", "G945-I_T")
sigs
```


<ol class=list-inline>
	<li>'G620_L'</li>
	<li>'G620_T'</li>
	<li>'G945-I_L'</li>
	<li>'G945-I_T'</li>
</ol>




```R
p <- list()
q <- list()


for(i in 1:length(sigs)){

sig <- sigs[i]
print(sig)
subset <- df[df$SampleID == sigs[i], ]
p[[i]] <- plot_sample_hex(subset, "PC1", "PC2", sigs[i])
q[[i]] <- plot_sample_hex_AUC(subset, "PC1", "PC2", sigs[i])
                                        
}
```

    [1] "G620_L"
    [1] "G620_T"
    [1] "G945-I_L"
    [1] "G945-I_T"



```R
pdf("~/Desktop/GSC_pairs_LineTumour.pdf", height = 3.5, width = 12)

ggarrange(q[[1]],
          q[[2]],
          q[[3]],
          q[[4]],
          ncol = 4,
          nrow = 1,
          common.legend = T,
          legend = "bottom"
    )

dev.off()
```




<strong>pdf:</strong> 2


----
## Plot examples
---


```R
sigs <- c("G1003-B_T", "G983-C_T", "BT89_L", "G729_L")
sigs
```


<ol class=list-inline>
	<li>'G1003-B_T'</li>
	<li>'G983-C_T'</li>
	<li>'BT89_L'</li>
	<li>'G729_L'</li>
</ol>




```R
p <- list()
q <- list()


for(i in 1:length(sigs)){

sig <- sigs[i]
print(sig)
subset <- df[df$SampleID == sigs[i], ]
p[[i]] <- plot_sample_hex(subset, "PC1", "PC2", sigs[i])
q[[i]] <- plot_sample_hex_AUC(subset, "PC1", "PC2", sigs[i])
                                        
}
```

    [1] "G1003-B_T"
    [1] "G983-C_T"
    [1] "BT89_L"
    [1] "G729_L"



```R
pdf("~/Desktop/GSC_pairs_LineTumour_exs.pdf", height = 7, width = 7)

ggarrange(q[[1]],
          q[[2]],
          q[[3]],
          q[[4]],
          ncol = 2,
          nrow = 2,
          common.legend = T,
          legend = "bottom"
    )

dev.off()
```




<strong>pdf:</strong> 2


----
## Plot all the rest
---


```R
### write function to plot sample counts

plot_sample_hex <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "black")

    
### plot patient counts
p <- ggplot(subset, aes_string(x=x, y=y, z="SampleID")) +
 geom_hex(bins=100) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) 


 p + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) 
    
    }

#plot score


plot_sample_hex_AUC <- function(dat, x, y, sig){ 

plot.title <- paste0(sig, " (", nrow(subset), " cells)")
    
    
title.col <- ifelse(unique(subset$SampleType) == "GSC", "black", "black")

q <- ggplot(subset, aes_string(x=x, y=y, z="Diff")) +
 stat_summary_hex(bins=100, fun = "median") +
#scale_fill_gradientn("Median Raw \nAUC Score", 
#                     colours = c(rev(brewer.pal(n = 8, name = "RdYlBu"))), 
#                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
#                    ) +
    
scale_fill_gradient2("Median Raw \nAUC Score",
                     low = "blue",
                     mid = "white",
                     high = "red",
                     limits=c(min(df[ ,"Diff"]),max(df[ ,"Diff"]))
                  ) +
 theme_classic(base_size=8) + ggtitle(plot.title) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(colour = title.col)) +  theme(legend.title = element_blank(), legend.position = "none") +
    xlab("") + ylab("")

 q + geom_path(data=Tumour_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "#225ea8", 
              lwd = 0.5
             ) +
 geom_path(data=GSC_contour_95, inherit.aes = FALSE,
              aes(x,y), 
              col = "black", 
              lwd = 0.5
             ) 
    
    
}




```


```R
sigs.exclude <- c("G1003-B_T", "G983-C_T", "BT89_L", "G729_L", "G620_L", "G620_T", "G945-I_L", "G945-I_T")
sigs <- sort(unique(df$SampleID)) #in this case sigs is patient 
sigs <- sigs[!sigs %in% sigs.exclude]
length(sigs)
```


43



```R
p <- list()
q <- list()


for(i in 1:length(sigs)){

sig <- sigs[i]
print(sig)
subset <- df[df$SampleID == sigs[i], ]
p[[i]] <- plot_sample_hex(subset, "PC1", "PC2", sigs[i])
q[[i]] <- plot_sample_hex_AUC(subset, "PC1", "PC2", sigs[i])
                                        
}
```

    [1] "BT127_L"
    [1] "BT147_L"
    [1] "BT48_L"
    [1] "BT67_L"
    [1] "BT73_L"
    [1] "BT84_L"
    [1] "BT94_L"
    [1] "G1003-A_T"
    [1] "G1003-C_T"
    [1] "G1003-D_T"
    [1] "G523_L"
    [1] "G549_L"
    [1] "G564_L"
    [1] "G566_L"
    [1] "G583_L"
    [1] "G637_L"
    [1] "G797_L"
    [1] "G799_L"
    [1] "G837_L"
    [1] "G851_L"
    [1] "G876_L"
    [1] "G885_L"
    [1] "G895_L"
    [1] "G910-A_T"
    [1] "G910-B_T"
    [1] "G910-C_T"
    [1] "G910-D_T"
    [1] "G910-E_T"
    [1] "G945-J_L"
    [1] "G945-J_T"
    [1] "G945-K_L"
    [1] "G945-K_T"
    [1] "G946-I_T"
    [1] "G946-J_L"
    [1] "G946-J_T"
    [1] "G946-K_L"
    [1] "G946-K_T"
    [1] "G967-A_T"
    [1] "G967-B_T"
    [1] "G967-C_T"
    [1] "G967-D_T"
    [1] "G983-A_T"
    [1] "G983-B_T"



```R


ml <- marrangeGrob(q, nrow=9, ncol=5)
ggsave("~/Desktop/GSC_pairs_LineTumour_TheRest.pdf", 
       ml, 
       height = 14, 
       width = 10, 
       dpi = 72, 
       limitsize = FALSE
      )


```


```R
?ggsave
```


```R

```
