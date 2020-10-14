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

suppressMessages(library(Seurat)) #v3.0
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))

##loads an RData file, and returns it with new name
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--seurat.obj",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object",
                                metavar= "character"
                               ),
                     make_option("--outFilePrefix",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files",
                                metavar= "character"
                              ))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
seurat.obj <- opt$seurat.obj
outFilePrefix <- opt$outFilePrefix

#### DEVELOPMENT
#seurat.obj <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/G523_L_res.0.2.RData"
#outFilePrefix <- "G523_L"



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
cells <- as.character(cells)
print(paste0(length(cells), " cells in object...."))

###chunk cell barcodes
splitSize <- length(cells) ##output each cell into own csv
print(paste0("Splitting cells into ", splitSize, " chunks..."))
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
cells <- chunk(cells, splitSize)
#names(cells) <- paste0("Chunk", as.numeric(names(cells))+1)
names(cells) <- cells


##############################################################
### 3) Save data as BamSlice compatible csv ##############################################################
print("")
print("********************")
print("Output csv file")
print(Sys.time())
###Output 2 csv files

### (1) The barcodes CSV files will each have one barcode entry per line,
### including the GEM well suffix (see GEM wells). Each such file will ### look something like: AAACGGGTCAAAGTGA-1
### output 1 csv file per chunk of cell barcodes
print("Writing out individual cell csvs....")
files <- c()
for (i in 1:length(names(cells))){
  #print(names(cells)[i])
  fileName <- paste0(getwd(), "/", outFilePrefix,  "_", names(cells)[i], ".csv")
  files[i] <- fileName
  sub <- as.matrix(cells[[i]])
  write.table(sub,
             file = fileName,
             quote = FALSE,
             col.names = FALSE,
             row.names = FALSE,
             sep = ","
             )
}

### (2) master csv file pointing to csv above
print("Writing out config csv....")
mastercsv <- matrix(c(paste0(outFilePrefix, ".", names(cells)), files),
                    ncol = 2
                  )

colnames(mastercsv) <- c("library_id", "barcodes_csv")
fileName2 <- paste0(outFilePrefix, "_BamSlice_Config.csv")
write.table(mastercsv,
           file = fileName2,
           quote = FALSE,
           col.names = TRUE,
           row.names = FALSE,
           sep = ","
           )

print("")
print("********************")
print("End of ExtractCellBarcodes_Seurat.r")
print(Sys.time())
print("********************")