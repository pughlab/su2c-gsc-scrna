##############################################################
#     Run Diffusion Map SU2C data using destiny R package    #
#                         L.Richards                         #
#                         Sept 2020                          #
##############################################################
### Destiny: https://bioconductor.org/packages/release/bioc/html/destiny.html
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/DiffusionMap/destiny

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Run Diffusion Map on GSCs
### 2) Save data
##############################################################

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
## Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/DiffusionMap_destiny.r --data /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata \
##  --outName GSC_destiny \
##  --downsample 100
##############################################################

library(optparse)
library(destiny) #v3.0.1
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
data <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata"
outName <- "GSCs_desinty"
downsample <- 100


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
      if(ncol(dat_sub) > downsample){
              cells <- c(cells, sample(colnames(dat_sub), downsample, replace = F))
      } else if (ncol(dat_sub) <= downsample){ #if not enough cell just use all
          cells <- c(cells, colnames(dat_sub))
      }
    }

    print(paste0("Downsampled to ", length(cells), " cells....."))
    dat <- SubsetData(dat,
                      cells.use = cells
                      )
    print(table(dat@ident))

}


##############################################################
# 3) Run Diffusion Map
##############################################################
print("****************************")
print("Run Diffusion Map using Destiny")
print(Sys.time())
print("****************************")

print(paste0("Running Diffusion map on ", ncol(dat@data), " cells....."))
pc.genes <- dat@calc.params$RunPCA$pc.genes
exprMatrix <- as.matrix(dat@data[pc.genes, ])
dim(exprMatrix)
### run diffusion map using destiny package
distance <-
sigma <-
dm <- DiffusionMap(exprMatrix,
                   suppress_dpt = TRUE,
                   n_eigs = 5,
                   k = 30, #same as seurat clustering
                   verbose = TRUE,
                   distance = "euclidean", #can also be "cosine" and "rankcor",
                  )

pdf("DiffusionMap_2900cells_euclidean.pdf")
plot(dm)
dev.off()



##############################################################
# 4) Save Results
##############################################################
print("****************************")
print("Save")
print(Sys.time())
print("****************************")

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


print("Plotting....")
### output plot
plot.name <- paste0(outName, "_", ncol(dat@data),"cells", "_DM.pdf")
pdf(plot.name, height = 8, width = 10)
#DMPlot(dat, group.by = "SampleID")
DMPlot(dat)
dev.off()



##############################################################
print("****************************")
print("END")
print(Sys.time())
print("****************************")
