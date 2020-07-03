##############################################################
#               Run Diffusion Map SU2C data                  #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Owen was running out of memory, lets see what superhimem can do
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/DiffusionMap

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Run Diffusion Map on GSCs
### 2) Save data
##############################################################
PatientID <- sapply(strsplit(rownames(BTSC_TumourCells@meta.data),"_"), `[`, 2)
SampleID <- paste(sapply(strsplit(rownames(BTSC_TumourCells@meta.data),"_"), `[`, 2),
                  sapply(strsplit(rownames(BTSC_TumourCells@meta.data),"_"), `[`, 3),
                  sep = "_"
                )
##############################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 96:00:00
## #SBATCH --mem=600G
## #SBATCH -p superhimem
## #SBATCH -c 50
## #SBATCH -N 1
## #SBATCH --account=pughlab
##
## module load R/3.6.1
##
## Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/DiffusionMap.r --data /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata \
##  --outName GSC \
##  --downsample 10
##############################################################

library(optparse)
suppressMessages(library("Seurat",
                         lib ="/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/bin/")
                ) #v2.3.4 needed to run Diffusion

#loads an RData file, and returns it
loadRData <- function(fileName){
                    load(fileName)
                    get(ls()[ls() != "fileName"])
                  }


##############################################################
# 1) Parse options
##############################################################
print("****************************")
print("Parse options")
print(Sys.time())
print("****************************")

option_list <- list(make_option("--data",
                                type = "character",
                                default = NULL,
                                help = "path to Seurat object",
                                metavar= "character"
                               ),
                     make_option("--outName",
                                type = "character",
                                default = NULL,
                                help = "will be appended to all output files",
                                metavar= "character"
                              ),
                     make_option("--downsample",
                                type = "integer",
                                default = NULL,
                                help = "Downsample to this many cells per ID; if you dont want to downsample enter 0",
                               metavar= "interger"
                            ))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

data <- opt$data
outName <- opt$outName
downsample <- opt$downsample

print(paste0("Dataset: ",data))
print(paste0("Output file prefix: ", outName))
print(paste0("Downsample cells: ", downsample))

######## DEVELOPMENT
#data <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata"
#outName <- "GSCs"
#downsample <- 10



##############################################################
# 2) Load and Reformat Data
##############################################################
print("****************************")
print("Loading Data")
print(Sys.time())
print("****************************")

if (grepl(".rds$", data)){
  dat <- readRDS(data)
}

if (grepl(".data$", data)){
  dat <- loadRData(data)
}

### set SampleID to ident
dat <- SetAllIdent(dat, id = "SampleID")

### downsample cells
if (downsample > 0){

    print(paste0("Downsampling data to max ", downsample, " cells per sample..."))

    ### get list of cell barcodes with downsampling
    cells <- c()
    samples <- as.character(unique(dat@ident))
    for (i in 1:length(samples)){
      dat_sub <- dat@data[ ,grepl(samples[i], colnames(dat@data))]
      cells <- c(cells, sample(colnames(dat_sub), downsample, replace = F))
    }

    print(paste0("Downsampled to ", length(cells), " cells....."))
    dat <- SubsetData(dat,
                      cells.use = cells
                      )

}



##############################################################
# 3) Run Diffusion Map
##############################################################
print("****************************")
print("Run Diffusion Map")
print(Sys.time())
print("****************************")

print(paste0("Running Diffusion map on ", ncol(dat@data), " cells....."))
dat <- RunDiffusion(dat,
                    genes.use = dat@var.genes, #to save time (~2765 genes)
                    #dims.use = dat@calc.params$RunTSNE$dims.use
                    dims.use = 1:10
                   )



##############################################################
# 4) Save Results
##############################################################
print("****************************")
print("Save")
print(Sys.time())
print("****************************")

print("Plotting....")
### output plot
plot.name <- paste0(outName, "_", ncol(dat@data),"cells", "_DM.pdf")
pdf(plot.name, height = 8, width = 10)
DMPlot(dat, group.by = "SampleID")
dev.off()

print("Saving....")
### save metadata
meta <- dat@meta.data
meta <- cbind(meta, dat@dr$dm@cell.embeddings) #add DM coordinates
meta.file <- paste0(outName, "_", ncol(dat@data),"cells", "_DM_Metadata.rds")
print(meta.file)
saveRDS(meta, meta.file)
### save seurat obj
seurat.file <- paste0(outName, "_", ncol(dat@data),"cells", "_DM_Seurat.rds")
print(seurat.file)
saveRDS(dat, file = seurat.file)


##############################################################
print("****************************")
print("END")
print(Sys.time())
print("****************************")
