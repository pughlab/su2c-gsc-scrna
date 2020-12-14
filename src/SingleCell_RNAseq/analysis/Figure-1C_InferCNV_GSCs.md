
---
# Run InferCNV on GSCs
---

**Parent Analysis Directory:** 

Use [InferCNV](https://github.com/broadinstitute/inferCNV) to assess CNV heterogeneity in GSC populations. 

Import Seurat objects from: '/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/input/regress_mito_nUMI'

---
# 1.0 Pilot Run
---

Use G637_L as the pilot run, AUCell showed that this sample had heterogeneity in chr10 deletion. Some cells appeared WT and one cluster appeared to have a DEL. 

**InferCNV parameters:**  
> Seurat LogNormalization with scale factor of 100,000  
> cutoff to 0.5  
> noise filter to 0.1  
> vis.bound to -1,1  


```R
library(infercnv)
library(Seurat)

load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/input/regress_mito_nUMI/G637_L_Seurat.RData")

print("#---Renormalize with scale.factor = 100,000---#")

cnv.dat <- LogNormalize(data = eb1S@raw.data,
                        scale.factor = 100000
                        )
print(dim(cnv.dat))
print(cnv.dat[1:5, 1:5])


#remove double rows

write.table(as.matrix(cnv.dat), 
            file = "G637_L_inferCNV_DGE.txt", 
            sep = "\t", 
            col.names = T, 
            row.names = T, 
            quote = F
           )
```


```R
xx <- read.table("~/pughlab/projects/BTSCs_scRNAseq/inferCNV/Cellranger_1.2.0_GRCh38_genPos.txt")
xu <- xx[!duplicated(xx$V1), ]

chrs <-c(1:22)
xu <- xu[xu$V2 %in% chrs, ]

write.table(xx, 
            "GenePos_GRCh38.txt",
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F
           )
```


```R
#remove KI270734.1., MT, and X from contig list in gene pos file


python /mnt/work1/users/pughlab/bin/inferCNV/src/gtf_to_position_file.py 
--attribute_name gene_name 
/mnt/work1/data/commondata/cellranger/1.2.0/hg19/genes/genes.gtf 
k562_gen_pos.txt
```


```R
module load R/3.2.2

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
--noise_filter 0.1 \
--output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/output \
--vis_bound_threshold " -1,1" \
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/G637_inferCNV_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/GenePos_GRCh38.txt
```


```R
#convert pdf to jpeg 

convert infercnv.pdf infercnv.jpg 
```

---
# 2.0 Use normal oligodendrocytes from brain tumours
---

Extract out "normal" brain cells to use as a reference for the BTSCs.   
The following samples have a high number of "normal" brain cells....

>910_A  
910_B  
910_C  
910_E  
945_J  
946_I  
946_J  
983_A    
983_B   
983_C     
1003_C    

This is a kind of annoying method becasuse you have to normalize the samples together before being able to call CNVs.


```R
library(Seurat)
```


```R
setwd("~/Desktop/Samwise/projects/Multiregional_GBM_scRNAseq/Manuscript/scClustViz/regress_nUMI_mito/input/")

G910_A <- load("910_A_Seurat.RData")
G910_A <- get(G910_A)
G910_A

G910_B <- load("910_B_Seurat.RData")
G910_B <- get(G910_B)
G910_B

G910_C <- load("910_C_Seurat.RData")
G910_C <- get(G910_C)
G910_C

G910_E <- load("910_E_Seurat.RData")
G910_E <- get(G910_E)
G910_E

G945_J <- load("945_J_Seurat.RData")
G945_J  <- get(G945_J )
G945_J 

G946_I <- load("946_I_Seurat.RData")
G946_I <- get(G946_I)
G946_I

G946_J <- load("946_J_Seurat.RData")
G946_J <- get(G946_J)
G946_J

G1003_C <- load("1003_C_Seurat.RData")
G1003_C  <- get(G1003_C )
G1003_C 
```


    An object of class seurat in project 910_B 
     15945 genes across 1610 samples.



    An object of class seurat in project 910_C 
     14691 genes across 893 samples.



    An object of class seurat in project 910_E 
     16082 genes across 964 samples.



    An object of class seurat in project 945_J 
     14829 genes across 1468 samples.



    An object of class seurat in project 946_I 
     13255 genes across 1388 samples.



    An object of class seurat in project 946_J 
     12674 genes across 1135 samples.



    An object of class seurat in project 1003_C 
     17529 genes across 1322 samples.



```R
MasterGBM <- MergeSeurat(object1 = G910_A, 
                         object2 = G910_B, 
                         add.cell.id1 = "910_A", 
                         add.cell.id2 = "910_B", 
                         project = "MultiRegionGBM"
                        )
```

    Warning message:
    “package ‘bindrcpp’ was built under R version 3.4.4”


```R
MasterGBM <- AddSamples(object = MasterGBM, 
                        new.data = G910_E@raw.data, 
                        add.cell.id = "910_E")
```


```R
MasterGBM <- AddSamples(object = MasterGBM, 
                        new.data = G945_J@raw.data, 
                        add.cell.id = "945_J")
```


```R
MasterGBM <- AddSamples(object = MasterGBM, 
                        new.data = G946_I@raw.data, 
                        add.cell.id = "946_I")
```


```R
MasterGBM <- AddSamples(object = MasterGBM, 
                        new.data = G946_J@raw.data, 
                        add.cell.id = "946_J")
```


```R
MasterGBM <- AddSamples(object = MasterGBM, 
                        new.data = G1003_C@raw.data, 
                        add.cell.id = "1003_C")
```


```R
mito.genes <- grep(pattern = "^MT-", x = rownames(x = MasterGBM@data), value = TRUE)
percent.mito <- Matrix::colSums(MasterGBM@raw.data[mito.genes, ])/Matrix::colSums(MasterGBM@raw.data)
```


```R
MasterGBM <- AddMetaData(object = MasterGBM, metadata = percent.mito, col.name = "percent.mito")
```


```R
MasterGBM <- NormalizeData(object = MasterGBM, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
```


```R
MasterGBM <- FindVariableGenes(object = MasterGBM, 
                               mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, 
                               x.high.cutoff = 3, 
                               y.cutoff = 0.5
                              )
length(x = MasterGBM@var.genes)
```


4295



![png](output_19_1.png)



```R
MasterGBM <- ScaleData(object = MasterGBM, vars.to.regress = c("nUMI", "percent.mito"))
```

    Regressing out: nUMI, percent.mito


    
    Time Elapsed:  1.48238751888275 mins

    Scaling data matrix



```R
MasterGBM <- RunPCA(object = MasterGBM, 
                    pc.genes = MasterGBM@var.genes, 
                    do.print = TRUE, 
                    pcs.print = 1:5, 
                    enes.print = 5)
```

    [1] "PC1"
     [1] "GPM6B"    "TUBA1A"   "CLU"      "PTN"      "MAP1B"    "MARCKSL1"
     [7] "SOX2"     "TUBB2B"   "NOVA1"    "FEZ1"     "UCHL1"    "S100B"   
    [13] "DDR1"     "CKB"      "CNN3"     "TSC22D4"  "PTPRZ1"   "SPARC"   
    [19] "MT3"      "CALD1"    "DNER"     "C1orf61"  "TTYH1"    "SCRG1"   
    [25] "SCD5"     "APP"      "PFN2"     "PCDH9"    "IGFBP2"   "PMP2"    
    [1] ""
     [1] "C1QC"     "C1QA"     "HLA-DPA1" "RGS1"     "HLA-DPB1" "HLA-DRB1"
     [7] "CD14"     "S100A9"   "MSR1"     "HLA-DQA1" "S100A4"   "IER3"    
    [13] "HMOX1"    "CD163"    "CCL4"     "REL"      "GPR183"   "CD83"    
    [19] "MAFB"     "HAMP"     "S100A8"   "CXCL8"    "C5AR1"    "FCGBP"   
    [25] "HLA-DRB5" "NR4A2"    "PLIN2"    "IFI30"    "CCL3L3"   "RGS2"    
    [1] ""
    [1] ""
    [1] "PC2"
     [1] "MAG"      "TF"       "PPP1R14A" "CLDN11"   "ERMN"     "KLK6"    
     [7] "EDIL3"    "MOG"      "APOD"     "PLP1"     "TUBB4A"   "NKX6-2"  
    [13] "CNDP1"    "ENPP2"    "SPOCK3"   "PTGDS"    "SEPT4"    "CAPN3"   
    [19] "TMEM144"  "CNTN2"    "HAPLN2"   "MOBP"     "ABCA2"    "CARNS1"  
    [25] "UGT8"     "SLAIN1"   "MAP7"     "GJB1"     "MAL"      "TMEM125" 
    [1] ""
     [1] "PTPRZ1"  "IGFBP2"  "FHL1"    "GPM6A"   "EGFR"    "IGFBP5"  "IGFBP7" 
     [8] "CNN3"    "ATP1B2"  "AQP4"    "FXYD6"   "VIM"     "AGT"     "FABP7"  
    [15] "ITM2C"   "CPE"     "SOX9"    "GAP43"   "SOX2"    "F3"      "PIFO"   
    [22] "RARRES2" "BCAN"    "NTRK2"   "OCIAD2"  "CYR61"   "MGST1"   "LANCL2" 
    [29] "NR2F1"   "SOCS2"  
    [1] ""
    [1] ""
    [1] "PC3"
     [1] "CAPS"          "C9orf24"       "PIFO"          "C1orf194"     
     [5] "FAM183A"       "TCTEX1D1"      "RSPH1"         "CFAP126"      
     [9] "TPPP3"         "C7orf57"       "IGFBP7-AS1"    "ZMYND10"      
    [13] "SPEF1"         "AQP4"          "AK1"           "LRRIQ1"       
    [17] "KIF9"          "SUSD2"         "GJA1"          "ANXA1"        
    [21] "AC015936.3"    "C9orf116"      "ROPN1L"        "IGFBP7"       
    [25] "DNAH9"         "FOXJ1"         "MORN5"         "NTRK2"        
    [29] "FOLR1"         "RP11-356K23.1"
    [1] ""
     [1] "SLC35E3"    "NUP107"     "B4GALNT1"   "MDM2"       "MARCH9"    
     [6] "CDK4"       "MIAT"       "SOX11"      "MIR181A1HG" "ELAVL4"    
    [11] "MARS"       "STMN2"      "DCTN2"      "ASIC4"      "OS9"       
    [16] "DLL3"       "MIRLET7BHG" "GRIA2"      "KIF5A"      "SEZ6L"     
    [21] "TSPAN31"    "MEG3"       "DCX"        "AGAP2-AS1"  "DTX3"      
    [26] "METTL1"     "MLLT11"     "TSFM"       "GPC2"       "MRPS17"    
    [1] ""
    [1] ""
    [1] "PC4"
     [1] "TOP2A"    "KIAA0101" "NUSAP1"   "UBE2C"    "BIRC5"    "PBK"     
     [7] "UBE2T"    "CENPF"    "MKI67"    "PTTG1"    "CDK1"     "TYMS"    
    [13] "MAD2L1"   "RRM2"     "FBLN1"    "TPX2"     "NUF2"     "AURKB"   
    [19] "PRC1"     "CDKN3"    "FAM64A"   "CCNB2"    "GTSE1"    "SGOL1"   
    [25] "IL13RA2"  "CDC20"    "SPC25"    "CENPA"    "RAD51AP1" "CCNB1"   
    [1] ""
     [1] "CYR61"     "VEGFA"     "CAV1"      "IGFBP3"    "SUSD2"     "ADM"      
     [7] "AKAP12"    "NDRG1"     "NRN1"      "KISS1R"    "DOK5"      "APLN"     
    [13] "ANGPTL4"   "CHI3L1"    "GADD45B"   "CTGF"      "EGR1"      "HILPDA"   
    [19] "DNAJB1"    "PPP1R15A"  "FGFR1"     "JUNB"      "EPAS1"     "LANCL2"   
    [25] "MEG3"      "NR2E1"     "TNFRSF12A" "FOSB"      "LGALS3"    "IER2"     
    [1] ""
    [1] ""
    [1] "PC5"
     [1] "C1QC"     "C1QA"     "CTSL"     "CD14"     "HLA-DPA1" "HLA-DRB1"
     [7] "S100A9"   "S100A8"   "HMOX1"    "HAMP"     "HLA-DPB1" "MSR1"    
    [13] "CD163"    "HLA-DQA1" "APOC2"    "GSN"      "MAFB"     "LGALS1"  
    [19] "FCGBP"    "C5AR1"    "GPNMB"    "CSTB"     "IFI27"    "SGK1"    
    [25] "NUPR1"    "IL1B"     "CXCL8"    "CPM"      "HLA-DRB5" "IFI30"   
    [1] ""
     [1] "IL32"       "CD2"        "CD3E"       "CD3D"       "TRAC"      
     [6] "CD52"       "GZMA"       "TRBC2"      "LCK"        "LTB"       
    [11] "CCL5"       "KLRB1"      "TRBC1"      "CD7"        "CD3G"      
    [16] "CD247"      "ACAP1"      "CD96"       "IL2RG"      "GZMK"      
    [21] "CTSW"       "AC092580.4" "NKG7"       "CST7"       "ITM2A"     
    [26] "IL7R"       "GZMH"       "GZMM"       "CD27"       "TUBA4A"    
    [1] ""
    [1] ""



```R
MasterGBM <- ProjectPCA(object = MasterGBM, do.print = FALSE)
```


```R
PCElbowPlot(object = MasterGBM)
```




![png](output_23_1.png)



```R
MasterGBM <- FindClusters(object = MasterGBM, 
                     reduction.type = "pca", 
                     dims.use = 1:15, 
                     resolution = 0.6, 
                     print.output = 0, 
                     save.SNN = TRUE
                    )
```


```R
MasterGBM <- RunTSNE(object = MasterGBM, 
                     dims.use = 1:15, 
                     do.fast = TRUE
                    )
```


```R
TSNEPlot(object = MasterGBM, pt.size = 0.4, do.label = T)
```


![png](output_26_0.png)



```R
FeaturePlot(object = MasterGBM, 
            features.plot = c("EGFR", "MOG", "MAG", "CD2", "CD3E", "CD14", "ITGAM"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
```


![png](output_27_0.png)



```R
#cluster 8 is normal brain cells
#

NormalBrain <- WhichCells(object = MasterGBM, ident = 8)
length(NormalBrain)
```


512



```R
NormalBrain.DGE <- MasterGBM@raw.data[ ,NormalBrain]
```


```R
write.table(as.matrix(NormalBrain.DGE),
            file = "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/Oligo_DGE.txt",
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = F
           )
```


```R
save(MasterGBM, file = "~/Desktop/Samwise/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/MultiRegion_Seurat.RData")
```

-----
## 3.0 Prepare InferCNV input files for each GSC line 
---

Run on Samwise


```R
##########################
# USER DEFINED VARIABLES #
##########################

path.to.reference <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/Oligo_DGE.txt"
tumour.dir <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript/scClustViz/input/regress_mito_nUMI/"

##########
# SCRIPT #
##########

print("")
print("#---Prepare input expression file for inferCNV---#")
print("")

suppressMessages(library(infercnv))
suppressMessages(library(Seurat))

print("")
print("#---Read in reference cell DGE---#")
print("")

print(path.to.reference)
Normal.rawDat <- read.table(path.to.reference)
print(dim(Normal.rawDat))
print(Normal.rawDat[1:5, 1:5])

files <- list.files(tumour.dir, pattern = "Seurat.RData")
samples <- gsub("_Seurat.RData", "", files)

for (i in 1:length(files)){
    
    print("")
    print("#---Read in tumour cells--#")
    print("")
    
    load.file <- paste(tumour.dir, files[i], sep = "")
    
    print("")
    print(samples[i])
    print(load.file)
    print("")

    load(load.file)
    
    print("")
    print("#---Combine reference and tumour cells into Seurat object--#")
    print("")


    CNV <- AddSamples(object = eb1S, 
                        new.data = Normal.rawDat, 
                        add.cell.id = "Normal")
    table(CNV@meta.data$orig.ident)
    
    print("")
    print("#---Renormalize with scale.factor = 100,000---#")
    print("")

    cnv.dat <- LogNormalize(data = CNV@raw.data,
                        scale.factor = 100000
                        )
    print(dim(cnv.dat))
    print(cnv.dat[1:5, 1:5])
    
    
    print("")
    print("#---Write out expression file---#")
    print("")

    DGE.name <- paste(samples[i], "_refOligos_DGE.txt", sep = "")
    print(DGE.name)

    write.table(as.matrix(cnv.dat), 
            file = DGE.name, 
            sep = "\t", 
            col.names = T, 
            row.names = T, 
            quote = F
           )
}


print("")
print("#---Write out comma delinited reference cell barcode file---#")
print("")

normal.BCs <- rownames(CNV@meta.data)[CNV@meta.data$orig.ident == "Normal"]

write.table(matrix(as.character(normal.BCs),nrow=1), 
            file = "Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

print("")
print("#---Input generation for inferCNV complete---#")
print("")

print(sessionInfo())
```

---
## 4.0 Run inferCNV
---

Make this into a bash script.  
Run for each sample. 
### Run 'InferCNV' on scRNA-seq data 
### L.Richards
### Novemeber 2018

#!/bin/bash
#
#$ -cwd

### Example Usage:
### qsub runInferCNV.sh SAMPLE_NAME
### qsub runInferCNV.sh

module load R/3.2.2

mkdir $1

/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/inferCNV/scripts/inferCNV.R --cutoff 0.5 \
--noise_filter 0.1 \
--output_dir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/output/$1 \
--vis_bound_threshold " -1,1" \
--log_file $1.log.txt \
--ref /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/input/Reference_Barcodes.txt \
/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/input/$1_refOligos_DGE.txt /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/input/GenePos_GRCh38.txt

convert infercnv.pdf infercnv.jpg 
---
## 5.0 Run on all cells at once....
---


```R
##########################
# USER DEFINED VARIABLES #
##########################

path.to.reference <- "~/pughlab/projects/BTSCs_scRNAseq/Manuscript/InferCNV/NormalBrainRef/Oligo_DGE.txt"
tumour.dir <- "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/PCA/Trans_programs/regress_mito_nUMI/merged_data/BTSC_25lines_regressMitoUMI_Seurat.RData"

##########
# SCRIPT #
##########

print("")
print("#---Prepare input expression file for inferCNV---#")
print("")

suppressMessages(library(infercnv))
suppressMessages(library(Seurat))

print("")
print("#---Read in reference cell DGE---#")
print("")

print(path.to.reference)
Normal.rawDat <- read.table(path.to.reference)
print(dim(Normal.rawDat))
print(Normal.rawDat[1:5, 1:5])


    
    print("")
    print("#---Read in tumour cells--#")
    print("")
    
    load.file <- tumour.dir
    
    print("")
   
    print(load.file)
    print("")

    load(load.file)
    
    print("")
    print("#---Combine reference and tumour cells into Seurat object--#")
    print("")
    
    print(dim(BTSC@data))

    CNV <- AddSamples(object = BTSC, 
                        new.data = Normal.rawDat, 
                        add.cell.id = "Normal")
    table(CNV@meta.data$orig.ident)
    
    print("")
    print("#---Renormalize with scale.factor = 100,000---#")
    print("")

    cnv.dat <- LogNormalize(data = CNV@raw.data,
                        scale.factor = 100000
                        )
    print(dim(cnv.dat))
    print(cnv.dat[1:5, 1:5])
    
    
    print("")
    print("#---Write out expression file---#")
    print("")

    DGE.name <- paste("All_BTSCs", "_refOligos_DGE.txt", sep = "")
    print(DGE.name)

    write.table(as.matrix(cnv.dat), 
            file = DGE.name, 
            sep = "\t", 
            col.names = T, 
            row.names = T, 
            quote = F
           )
}


print("")
print("#---Write out comma delinited reference cell barcode file---#")
print("")

normal.BCs <- rownames(CNV@meta.data)[CNV@meta.data$orig.ident == "Normal"]

write.table(matrix(as.character(normal.BCs),nrow=1), 
            file = "Reference_Barcodes.txt",
            sep=",",
            row.names=FALSE,
            col.names=FALSE,
            quote = F
           )

print("")
print("#---Input generation for inferCNV complete---#")
print("")

print(sessionInfo())
```


```R
All_BTSCs
```
