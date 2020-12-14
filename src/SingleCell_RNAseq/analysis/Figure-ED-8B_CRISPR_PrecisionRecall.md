
---
# Plot Precision Recall Curves for CRISPR Data
----


From Graham: "Precision-recall curves are often included to show screen performance. I can either plot these for you, or provide the necessary files if you want to do it yourself to keep the aesthetics consistent."

Figure 2 of BAGEL paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1015-8

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/PrecisionRecall

----
## 1.0 Plot Precision-Recall curves
----


```R
setwd("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CRISPR/PrecisionRecall")

files <- list.files("./PR_files_for_SU2C_manuscript/", pattern = ".PR")
files

names <- gsub(".PR", "_L", files)
names[7] <- "G691_L"



```


<ol class=list-inline>
	<li>'BT67.PR'</li>
	<li>'G361.PR'</li>
	<li>'G523.PR'</li>
	<li>'G549.PR'</li>
	<li>'G583.PR'</li>
	<li>'G620.PR'</li>
	<li>'G691R.PR'</li>
	<li>'G729.PR'</li>
	<li>'G809.PR'</li>
</ol>




```R
pdf("~/Desktop/CRISPR_PR_Curves.pdf", height = 5, width = 5)

i <- 1
#cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','black')
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','black','#a65628','#f781bf','#999999')

PR <- read.table(paste0("./PR_files_for_SU2C_manuscript/", files[i]),
                 header = T, 
                 sep = "\t"
                )
    
    plot(PR$Recall,
     PR$Precision,
     type = "l",
     lwd = "2",
     #main = gsub(".PR", "_L", files[i]),
     xlab = "Recall",
     ylab = "Precision",
     col = cols[i]
    )

for(i in 2:length(files)){
    
    
   PR <- read.table(paste0("./PR_files_for_SU2C_manuscript/", files[i]),
                 header = T, 
                 sep = "\t"
                )
    
 lines(PR$Recall,
       PR$Precision,
       col = cols[i],
       lwd = 2,
       #type = "S"
      )
    
}

legend("bottomleft",
       legend = names,
       col = cols,
       lty = 1,
       lwd = 2
      
      )


dev.off()
```


<strong>pdf:</strong> 2

