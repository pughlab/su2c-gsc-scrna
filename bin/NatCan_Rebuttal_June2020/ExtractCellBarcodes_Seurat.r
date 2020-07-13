##############################################################
#             Extract cell barcodes from Seurat Object       #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### This is part of the scRNA mutation calling pipeline
### Reference: https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/using/bamslice#configuration

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Extract cell barcodes from Seurat Object
### 2) Format into Bamslice compatible csv
### 3) Save
##############################################################

library(Seurat) #v3.0

##loads an RData file, and returns it with new name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

seurat.obj <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/G523_L_res.0.2.RData"
outFilePrefix <- "G523_L"

##############################################################
### 1) Load seurat object ##############################################################
print("")
print("********************")
print("Load Data")
print(Sys.time())

if(grepl(".rds$", tolower(seurat.obj))){
  dat <- readRDS(seurat.obj)
} else if (grepl(".rdata$", tolower(seurat.obj))){
  dat <- loadRData(seurat.obj)
}

print("Updating seurat object to v3....")
dat <- UpdateSeuratObject(dat)


##############################################################
### 2) Extract cell barcodes ##############################################################
print("")
print("********************")
print("Extract Cell Barcodes")
print(Sys.time())

cells <- colnames((GetAssayData(dat, slot = "counts")))
cells <- paste0(cells, "-1") ### add 10X GEM group to cell barcode
print(paste0(length(cells), " cells in object...."))


##############################################################
### 3) Save data as BamSlice compatible csv ##############################################################
print("")
print("********************")
print("Output csv file")
print(Sys.time())

### The barcodes CSV files will each have one barcode entry per line,
### including the GEM well suffix (see GEM wells). Each such file will ### look something like: AAACGGGTCAAAGTGA-1
fileName <- paste0(outFilePrefix, "_CellBarcodes.csv")
fileConn <- file(fileName)
writeLines(cells, fileConn)
close(fileConn)
