
---
# Run Seurat Processing Pipeline
---

**L.Richards**  
Pugh Lab  

Process and run new BTSC samples with Seurat. 
Preparing for input into scClustViz:
> run up to screen plot to determine number of PCs


/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/SeuratProcessing

---
## 1.0 runSeurat_preprocess.sh 
---

Bash script to run 
#!/bin/bash
#
#$ -cwd

#first argument is sample prefix

module load R/3.5.0
module load umap

Rscript Seurat_preprocess.R $1

```R
qsub runSeurat_preprocess.sh G945-J_L
qsub runSeurat_preprocess.sh G945-I_L
qsub runSeurat_preprocess.sh G945-K_L
qsub runSeurat_preprocess.sh G946-J_L
qsub runSeurat_preprocess.sh G946-K_L
```


```R
############################################
# Run Seurat Pipeline on a single sample  
# L.Richards
# June 2019
############################################

#### Load packages
 
suppressMessages(library("Seurat", 
                         lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/"),
                 library("optparse")
                )


args <- commandArgs(TRUE)
sample <- args[1]
print(sample)

#DEVELOPMENT
#sample <- "G946-K_L"


print("")
print("-------")
print(sample)
print(Sys.time())
print("-------")

print("Step 1: Load QC'ed input data")
print(Sys.time())

file.path <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/QC/"
file <- paste0(sample, "_QCFiltered.RData")
print(file)

load.file <- paste0(file.path, file)
print(load.file)
load(load.file)
print(dim(BTSC@data))

print("Step 2: Normalize and Scale Data")
print(Sys.time())

BTSC <-   NormalizeData(BTSC, 
                          assay.type = "RNA",
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000,
                          display.progress = TRUE
                         )
mito.genes <- grep(pattern = "^MT-", 
                   x = rownames(x = BTSC@raw.data), 
                   value = TRUE)

percent.mito <- Matrix::colSums(BTSC@raw.data[mito.genes, ])/Matrix::colSums(BTSC@raw.data)

BTSC <- AddMetaData(object = BTSC, 
                          metadata = percent.mito, 
                          col.name = "percent.mito")
BTSC <- ScaleData(object = BTSC, vars.to.regress = c('percent.mito', "nUMI"))

print("Step 3: Find Variable Genes")
print(Sys.time())

BTSC <- FindVariableGenes(object = BTSC, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, 
                            x.high.cutoff = 3, 
                            y.cutoff = 0.5
                           )


print("Step 4: Run PCA on all genes")
print(Sys.time())

BTSC <- RunPCA(object = BTSC, 
                 pc.genes = rownames(BTSC@data), #all genes
                 do.print = FALSE, 
                 pcs.compute = 100
                )


print("Step 5: Output scree plot")
print(Sys.time())

scree.name <- paste0(sample, "_PCscreeplot.pdf")
print(scree.name)

pdf(scree.name)
PCElbowPlot(object =BTSC, num.pc = 25)
dev.off()

print("Step 6: Save file")
print(Sys.time())

save.file <- paste0(sample, "_scClustViz_input.Rdata")
print(save.file)
save(BTSC, file = save.file)

print("END")
print(Sys.time())

```
