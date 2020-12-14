
---
# Regressing out cell cycle effect in BTSCs
---

L.Richards  
Analysis Dir: /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression
#!/bin/bash
#
#$ -cwd

module load R/3.5.0
module load umap

Rscript CellCycleRegression_Alternative.R $1
### Original BTSC Cohort


```R
library("optparse") 
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))



args <- commandArgs(TRUE)
sample <- args[1]
print(sample)

print("")
print("-------")
print(sample)
print(Sys.time())
print("-------")

print("Step 1: Load QC'ed input data and gene signatures")
print(Sys.time())

file.path <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/input/regress_mito_nUMI/"
file <- paste0(sample, "_Seurat.RData")
print(file)

load.file <- paste0(file.path, file)
print(load.file)

BTSC <- local({
        load(load.file)
        stopifnot(length(ls())==1)
        environment()[[ls()]]
})

print(dim(BTSC@data))

cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

print("Step 2: Assign Cell Cycle Scores")
print(Sys.time())

BTSC <- CellCycleScoring(object = BTSC, 
                 s.genes = s.genes, 
                 g2m.genes = g2m.genes, 
                 set.ident = TRUE
                )
print(table(BTSC@ident))

print("Step 3: Regress CC difference, nUMI and percent.mito")
print(Sys.time())

BTSC@meta.data$CC.Difference <- BTSC@meta.data$S.Score - BTSC@meta.data$G2M.Score

BTSC <- ScaleData(object = BTSC, 
                  vars.to.regress = c("CC.Difference", "nUMI", "percent.mito"), 
                  display.progress = TRUE
                 )


print("Step 4: Run PCA on all genes EXCEPT ribo")
print(Sys.time())

ribo.genes <- c(rownames(BTSC@data)[grep("^RP[1:9]", rownames(BTSC@data))],
                rownames(BTSC@data)[grep("^RP[L,S]", rownames(BTSC@data))]
               ) #667 genes


BTSC <- RunPCA(BTSC, 
               pc.genes = rownames(BTSC@data)[!rownames(BTSC@data) %in% ribo.genes],
               pcs.compute = 100,
               genes.print = 10
              )

print("Step 5: Output scree plot")
print(Sys.time())

scree.name <- paste0(sample, "_CCregressed_noRibo_PCscreeplot.pdf")
print(scree.name)

pdf(scree.name)
PCElbowPlot(object =BTSC, num.pc = 30)
dev.off()

print("Save data")
print(Sys.time())

save.file <- paste0(sample, "_CCregressed_Seurat.Rdata")
print(save.file)
save(BTSC, file = save.file)

print("End of Analysis")
print(Sys.time())
```

### New BTSC Samples
#!/bin/bash
#
#$ -cwd

module load R/3.5.0
module load umap

Rscript CellCyle_Regression_new.R $1

```R
library("optparse") 
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))



args <- commandArgs(TRUE)
sample <- args[1]
print(sample)

print("")
print("-------")
print(sample)
print(Sys.time())
print("-------")

print("Step 1: Load QC'ed input data and gene signatures")
print(Sys.time())

file.path <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/SeuratProcessing/"
file <- paste0(sample, "_scClustViz_input.Rdata")
print(file)

load.file <- paste0(file.path, file)
print(load.file)

BTSC <- local({
        load(load.file)
        stopifnot(length(ls())==1)
        environment()[[ls()]]
})

print(dim(BTSC@data))

cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

print("Step 2: Assign Cell Cycle Scores")
print(Sys.time())

BTSC <- CellCycleScoring(object = BTSC, 
                 s.genes = s.genes, 
                 g2m.genes = g2m.genes, 
                 set.ident = TRUE
                )
print(table(BTSC@ident))

print("Step 3: Regress CC difference, nUMI and percent.mito")
print(Sys.time())

BTSC@meta.data$CC.Difference <- BTSC@meta.data$S.Score - BTSC@meta.data$G2M.Score

BTSC <- ScaleData(object = BTSC, 
                  vars.to.regress = c("CC.Difference", "nUMI", "percent.mito"), 
                  display.progress = TRUE
                 )


print("Step 4: Run PCA on all genes EXCEPT ribo")
print(Sys.time())

ribo.genes <- c(rownames(BTSC@data)[grep("^RP[1:9]", rownames(BTSC@data))],
                rownames(BTSC@data)[grep("^RP[L,S]", rownames(BTSC@data))]
               ) #667 genes


BTSC <- RunPCA(BTSC, 
               pc.genes = rownames(BTSC@data)[!rownames(BTSC@data) %in% ribo.genes],
               pcs.compute = 100,
               genes.print = 10
              )

print("Step 5: Output scree plot")
print(Sys.time())

scree.name <- paste0(sample, "_CCregressed_noRibo_PCscreeplot.pdf")
print(scree.name)

pdf(scree.name)
PCElbowPlot(object =BTSC, num.pc = 30)
dev.off()

print("Save data")
print(Sys.time())

save.file <- paste0(sample, "_CCregressed_Seurat.Rdata")
print(save.file)
save(BTSC, file = save.file)

print("End of Analysis")
print(Sys.time())
```

### Pass to scClustViz
#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash

#### Example usage: qsub run_scClustViz_CCregressed.sh G549_L 10

module load R/3.5.0
module load umap

Rscript /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/scClustViz_v1.2.1_CC.R --sample $1 --output.dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/ --input.seurat /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/$1_CCregressed_Seurat.Rdata --pseudocount 1 --maxPCt $2 --perplexity 30 --max.resolution 1

```R
###########################################################
#          Script to run scClustViz on Seurat Input       #
#                         L.Richards                      #
#                          June 2019                      #
###########################################################

suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))

suppressMessages(library(scales))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(cluster))
suppressMessages(library(pbapply))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(scClustViz))

print(sessionInfo())

###################
#  PARSE OPTIONS  #
###################

option_list <- list(

make_option("--sample",
            type = "character",
            default = NULL,
            help = "Name of sample to append to output files",
            metavar= "character"),

make_option("--output.dir",
            type = "character",
            default = NULL,
            help = "Path to directory to output results",
            metavar= "character"),

make_option("--input.seurat",
            type = "character",
            default = NULL,
            help = "Path to seurat object saved as .RData",
            metavar= "character"),

make_option("--pseudocount",
            type = "integer",
            default = 1,
            help = "Pseudocount added to all log-normalized values in your input data. Most methods use a pseudocount of 1 to eliminate log(0) errors. Default =1",
            metavar= "integer"),

make_option("--maxPCt",
            type = "integer",
            default = 10,
            help = "Significant PCs to use in analysis based on scree plot. Default = 10",
            metavar= "integer"),

make_option("--perplexity",
            type = "integer",
            default = 30,
            help = "Perplexity to run tSNE. Script will output variety of tSNEs with varying perplex to assess shape. Default = 30.",
            metavar= "integer"),

make_option("--max.resolution",
            type = "integer",
            default = 1,
            help = "Resolution max, increase by 0.1 increments. ",
            metavar= "integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dataName <- opt$sample
output <-  opt$output.dir
dataRData <- opt$input.seurat
pseudocount <- opt$pseudocount
maxPCt <- opt$maxPCt
perplex <- opt$perplexity
max.res <- opt$max.resolution


dataPath <- paste(output, dataName, sep ="") #output.dir
PCuse <- seq(1,maxPCt)
resVal <- resVal <- seq(from = 0.1, to = max.res, by = 0.1)

########### DEVELOPMENT ##############


#dataName <- "G523_L"
#output <-  "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/output/"
#dataPath <- paste(output, dataName, sep ="") #output.dir
#dataRData <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/input/G523_L_Seurat.RData"  #input data
#exponent <- exp(1)
#pseudocount <- 1
#logFCthresh <- 1
#WRSTalpha <- 0.01
#maxPCt <- 10
#PCuse <- seq(1,maxPCt)
#perplex <- 30

#scClust_viz_prep <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/output/scripts/run_scClustViz_clusterDE.R"


#####################################

print("Running scClustViz with the following options:")
print(opt)

print("$dataPath")
dataPath <- paste(output, dataName, sep ="") #output.dir
print(dataPath)

print("$PCuse")
PCuse <- seq(1,maxPCt)
print(PCuse)

print("$resVal")
print(resVal)

########
# MAIN #
########

print("")
print("#----Setting up results directory----#")
print("")


setwd(output)
dir.create(dataName)
setwd(dataPath)
print(getwd())

print("")
print("#----Reading in Seurat Object----#")
Sys.time()
print("")

#name object BTSC

BTSC <- local({
        load(dataRData)
        stopifnot(length(ls())==1)
        environment()[[ls()]]
})


print(dim(BTSC@data))
print(BTSC@data[1:5, 1:5])

print(paste("#----Run tSNE with perplexity =", perplex, "----#", sep = ""))
Sys.time()

BTSC <- RunTSNE(BTSC,
                dims.use = PCuse,
                perplexity = perplex,
                reduction.use = "pca",
                do.fast = TRUE
               )

print(paste("#----Run UMAP", "----#", sep = ""))
Sys.time()

BTSC <- RunUMAP(BTSC,
                dims.use = PCuse,
                reduction.use = "pca"
               )

print(paste("#----Run Diffusion Map", "----#", sep = ""))
Sys.time()

BTSC <- RunDiffusion(BTSC,
                     dims.use = PCuse,
                     reduction.use = "pca"
                    )

print("")
print("#----Sequentially cluster and perform DE analysis for a series of resolutions----#")
Sys.time()
print("")

for (i in 1:length(resVal)){

    print(resVal[i])
    print(Sys.time())
    BTSC <- FindClusters(BTSC,
                     reduction.type="pca",
                     dims.use=PCuse,
                     k.param=30,
                     print.output=F,
                     resolution=resVal[i],
                     algorithm=3,
                     n.start=100,
                     n.iter=100,
                     save.SNN=T
                    )
}

print("")
print("#----Save results from sequential clustering----#")
Sys.time()
print("")

filename <- paste(dataName, "_ClustVizInput.RData", sep = "")
print(filename)
save(BTSC, file=filename)

print("")
print("#----Run scClustViz----#")
Sys.time()
print("")

your_cluster_columns <- grepl("res[.0-9]+$",
                              names(getMD(BTSC)))

your_cluster_results <- BTSC@meta.data[ ,your_cluster_columns]

### check is there is a clustering solution with only one option

print("Check is any clustering solutions have only one cluster")
print("We need to remove them because they throw an error")

filtercluster <- c()

for (i in 1:ncol(your_cluster_results)){
    
    filtercluster[i] <- length(names(table(your_cluster_results[i]))) > 1
    
}

print(filtercluster)

your_cluster_results <- your_cluster_results[ ,filtercluster]


sCVdata_list <- CalcAllSCV(
  inD=BTSC,
  clusterDF=your_cluster_results,
  #assayType=NULL, #specify assay slot of data
  DRforClust="pca",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=TRUE, 
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

print("")
print("#----Save scClustViz Results----#")
Sys.time()
print("")

save.file <- paste0(dataName, "_scClustViz_CCregressed_output.Rdata")

save(BTSC,
     sCVdata_list,
     file=save.file)

print("")
print("#----END OF ANALYSIS----#")
Sys.time()
print("")
```

----
## Visualize scClustViz Results
----


```R
library(scClustViz)

sample <- "G583_L"
load.file <- paste0("~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/",
               sample, "/", sample, "_scClustViz_CCregressed_output.Rdata"
              )
print(load.file)

runShiny(filePath=load.file,
         outPath="./",
         #imageFileType="jpeg"
)
```

---
## Extract optimal clustering solution
---
Want to create a new seurat object with the resolution we want

Metadata
> Opt.Res = optimal resolution for IntraBTSC clustering  
> Cluster.ID = clustering ID from seurat plus 1 (ie. 0 becomes C1)  

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/opt_solution


```R
suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))

#load in matrix with solution info

opt.clust <- read.csv("newClustSolutions.csv")


```


```R
for (i in 1:nrow(opt.clust)){

print("#---Set up sample parameters---#")
print(Sys.time())

sample <- opt.clust$Sample[i]
res <- opt.clust$Resolution[i]
file.path <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/"
load.file <- paste(file.path, sample, "/", sample, "_scClustViz_CCregressed_output.Rdata", sep ="")
print(sample)
print(res)
print(load.file)

print("#---Load Seurat Object---#")
print(Sys.time())
load(load.file)


print("#---Remove non-optimal clustering solution from metadata---#")
print(Sys.time())

res.name <- paste("res.", res, sep = "")
print("Optimal resolution is:")
print(res.name)

BTSC@meta.data <- BTSC@meta.data[ ,c(1:9, grep(res.name, colnames(BTSC@meta.data)))]
BTSC@meta.data$Opt.Resolution <- res.name

#change cluster names

colnames(BTSC@meta.data)[grep(res.name, colnames(BTSC@meta.data))] <- "Cluster.ID"
BTSC@meta.data$Cluster.ID <- paste0("C", as.numeric(BTSC@meta.data$Cluster.ID)+1)

print(head(BTSC@meta.data))
    
BTSC@ident <- factor(BTSC@meta.data$Cluster.ID)
print(BTSC@ident[1:10])

print("#---Saving file---#")
print(Sys.time())

save.file <- paste(sample, "_", res.name, ".RData", sep ="")
print(save.file)

save(BTSC, file = save.file)

}
```

----
## Plot dimensionality reduction with optimal clusters for BTSCs
----


```R
suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))

path <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/opt_solution/"

files <- list.files(path, pattern = "RData")
print(files)

samples <- gsub('.{14}$', '', files)
samples


for (i in 1:length(samples)){
    
    
    sample <- samples[i]
    print("----------")
    print(sample)
    
    load(files[i])
    
    dat <- cbind(BTSC@dr$tsne@cell.embeddings, 
             BTSC@dr$umap@cell.embeddings,
             BTSC@dr$dm@cell.embeddings,
             BTSC@meta.data
            )
    head(dat)
    
    plot.title <- paste(sample, "     ", nrow(dat), "cells")
    print(plot.title)
    
    
    sample_tSNE <- ggplot(dat, aes(x=tSNE_1, y=tSNE_2, color=Cluster.ID)) + 
                   geom_point(alpha = 0.6, size = 1.5, pch = 16) +  
                   labs(x = "tSNE 1", y = "tSNE 2", title = plot.title) +
                   scale_colour_brewer(palette = "Dark2") + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

    sample_umap <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color=Cluster.ID)) + 
                   geom_point(alpha = 0.6, size = 1.5, pch = 16) +  
                   labs(x = "UMAP 1", y = "UMAP 2", title = plot.title) +
                   scale_colour_brewer(palette = "Dark2") + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

    sample_dm <- ggplot(dat, aes(x=DM1, y=DM2, color=Cluster.ID)) + 
                   geom_point(alpha = 0.6, size = 1.5, pch = 16) +  
                   labs(x = "Diffusion Map 1", y = "Diffusion Map 2", title = plot.title) +
                   scale_colour_brewer(palette = "Dark2") + 
                   theme_bw() + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1)))

plot.name <- paste0(sample, "_tSNE_UMAP_DM.pdf")

pdf(plot.name, height = 5, width = 5.5)
print(sample_tSNE)
print(sample_umap)
print(sample_dm)
dev.off()
    
}
```

---
## Make meta dataframe with cluster intra cluster info
---

This will be added to the big dataset for downstream analysis.


```R
suppressMessages(library("Seurat", lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" ))

path <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/opt_solution/"

files <- list.files(path, pattern = "RData")
print(files)

samples <- gsub('.{14}$', '', files)
samples

i <- 1
load(files[i])

intra.meta <- data.frame(BTSC@meta.data)
rownames(intra.meta) <- paste0(samples[i], "_", rownames(intra.meta))


for (i in 2:length(samples)){
    
    sample <- samples[i]
    print("----------")
    print(sample)
    
    print(files[i])
    load(files[i])
    
    meta <- data.frame(BTSC@meta.data)
    rownames(meta) <- paste0(sample, "_", rownames(meta))
    print(head(meta))
    
    intra.meta <- rbind(intra.meta, meta)
    
}

save(intra.meta, file = "BTSC_IntrCluster_metadata.RData")
```


```R
#append this to global object and save

#load big BTSC object
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")

#load metadata
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/CellCyle_Regression/scClustViz/opt_solution/BTSC_IntrCluster_metadata.RData")
```


```R
intra.meta <- intra.meta[ ,c("Cluster.ID", "Opt.Resolution")]

BTSC <- AddMetaData(BTSC,
                   metadata = intra.meta
                   )

BTSC@meta.data$IntraBTSC.ID <- paste0(BTSC@meta.data$SampleID, "_", BTSC@meta.data$Cluster.ID)
```


```R
#Average by sample
BTSC <- SetIdent(BTSC, ident.use = BTSC@meta.data$SampleID)
BTSC_SampleAvg <- AverageExpression(BTSC,
                                    show.progress = TRUE
                                   )
 
BTSC_SampleAvg_scaled <- AverageExpression(BTSC,
                                           show.progress = TRUE,
                                           use.scale = TRUE
                                          )
save(BTSC_SampleAvg, file = "BTSC_SampleAveraged.Rdata")
save(BTSC_SampleAvg_scaled, file = "BTSC_SampleAveraged_scaled.Rdata")
     
     
#Average by IntraBTSC Cluster
BTSC <- SetIdent(BTSC, ident.use = BTSC@meta.data$IntraBTSC.ID)
BTSC_IntraClusterAvg <- AverageExpression(BTSC,
                                    show.progress = TRUE
                                   )

BTSC_IntraClusterAvg_scaled <- AverageExpression(BTSC,
                                    show.progress = TRUE, use.scale = TRUE
                                   )

save(BTSC_IntraClusterAvg, file = "BTSC_IntraClusterAvgeraged.Rdata")
save(BTSC_IntraClusterAvg_scaled, file = "BTSC_IntraClusterAvgeraged_scaled.Rdata")

save(BTSC, 
     file = "~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata"
    )
```
############################
Reprocessed SU2C BTSC scRNAseq Cohort
June 2019
L.Richards
############################

"Global_SU2C_BTSCs_CCregressed_noRibo.Rdata" 
- Seurat object with all cells across samples
- @data slot is LogNormalized counts
- @scale.data slot is gene centered data with nUMI, percent.mito and CC.difference regressed 
- Cell cycle regression method details here (see "Alternate Workflow")
    - https://satijalab.org/seurat/v2.4/cell_cycle_vignette.html
    
"BTSC_SampleAveraged.Rdata"
- LogNormalized expression data averaged by SampleID

"BTSC_SampleAveraged_scaled.Rdata"
- Scaled expression data averaged by SampleID

"BTSC_IntraClusterAvgeraged.Rdata"
- LogNormalized expression data averaged by IntraBTSC cluster ID

"BTSC_IntraClusterAvgeraged_scaled.Rdata"
- Scaled expression data averaged by  by IntraBTSC cluster ID
----
## Compare cluster number between scClustViz runs
----


```R
dat <- read.csv("~/Downloads/scClustVizComparison.csv")
head(dat)
```


```R
plot(dat$scClustViz_v0.5.0,
     dat$scClustViz_v1.2.1,
     xlab = "scClustViz_v0.5.0",
     ylab = "scClustViz_v1.2.1",
     main = "Cluster Number",
     xlim = c(1,8),
     ylim = c(1,8)
    )
 text(dat$scClustViz_v0.5.0, 
      dat$scClustViz_v1.2.1, 
      labels=dat$Sample, 
      cex= 0.7,
      pos=1
     )


plot(dat$scClustViz_v0.5.0,
     dat$scClustViz_v1.2.1_CC,
     xlab = "scClustViz_v0.5.0",
     ylab = "scClustViz_v1.2.1 Cell Cycle Regressed",
     main = "Cluster Number",
     xlim = c(1,8),
     ylim = c(1,8)
    )
 text(dat$scClustViz_v0.5.0, 
      dat$scClustViz_v1.2.1_CC, 
      labels=dat$Sample, 
      cex= 0.7,
      pos=1
     )

plot(dat$scClustViz_v1.2.1,
     dat$scClustViz_v1.2.1_CC,
     xlab = "scClustViz_v1.2.1",
     ylab = "scClustViz_v1.2.1 Cell Cycle Regressed",
     main = "Cluster Number",
     xlim = c(1,8),
     ylim = c(1,8)
    )
 text(dat$scClustViz_v1.2.1, 
      dat$scClustViz_v1.2.1_CC, 
      labels=dat$Sample, 
      cex= 0.7,
      pos=1
     )
```
