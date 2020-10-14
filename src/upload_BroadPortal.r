###########################################################
#            Generate files for upload to Broad           #
#                        scRNA portal                     #
#                         L.Richards                      #
#                         June 2020                       #
###########################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Read in seurat object (v2)
### 2) Extract and format data file
### 3) Extract and format meta.data file
### 4) Extract and format clusters files
### 5) Save data
###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
###
###
###
###
###
###########################################################

library(Seurat)
library(optparse)
library(data.table)
library(dplyr)
library(taRifx)

###########################################################
### USER DEFINE VARIABLES
###########################################################

#path to seurat object to extract data from
seurat.obj <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Hypoxia/AllTumour_AUCell_Seurat.Rdata"

#append string to output files
out.name <- "SU2C_scRNA_TumourOnly"

#path to text file with names of metadata columns to include
#each column name on own line
#if "ALL", use all columns in metadatafile
meta.columns <- "ALL"

#extract tSNE Coordinates
tsne <- FALSE

#extract UMAP Coordinates
umap <- TRUE



###########################################################
### 1) Load data
###########################################################
print("")
print("***************************")
print("Loading Seurat Object...")
print(Sys.time())
print("***************************")

#function to load and rename Rdata Object
load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

dat <- load_obj(seurat.obj)
dat <- UpdateSeuratObject(dat) #update if v2, check structure if v3

###########################################################
### 2) Extract normalized expression
###########################################################
print("")
print("***************************")
print("Normalized expression matrix")
print(Sys.time())
print("***************************")

print("Size of expression matrix....")
print(dim(dat@assays$RNA@data))

print("Extracting normalized expression matrix...")
norm.exp <- data.frame(as.matrix(dat@assays$RNA@data))
norm.exp <- tibble::rownames_to_column(norm.exp, "GENE")
print(norm.exp[1:5, 1:5])

file1 <- paste0(out.name, "_NormalizedExpression.txt")
print(paste0("Writing out expression matrix...", file1))
write.table(norm.exp,
            file=file1,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
          )


###########################################################
### 3) Extract metadata
###########################################################
print("")
print("***************************")
print("Metadata")
print(Sys.time())
print("***************************")

print("Extract Metadata matrix...")
meta <- data.frame(dat@meta.data)
meta <- remove.factors(meta)

if(meta.columns == "ALL"){

  include <- colnames(meta)
  print("Using all columns in metadata...")
  print(cat(include ,sep = "\n"))


} else if (grepl(".txt", meta.columns)){

  include <- readLines(meta.columns)
  print("Using pre-specified set of metadata columns...")
  print(cat(include ,sep = "\n"))

}

#subset metadata matrix
meta <- meta[, include]
print(dim(meta))

##format cluster column to be group
## search "res." which is seurat default
cluster.idx <- grep("^res\\.", colnames(meta))
colnames(meta)[cluster.idx] <- paste0("Cluster_", colnames(meta)[cluster.idx])
meta[ ,cluster.idx] <- paste0("C", meta[ ,cluster.idx])


print("Formatting...")
#format for Broad Guidelines
classes <- as.vector(unlist(lapply(meta, class)))
TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric", "group"))
meta <- data.frame(setDT(meta, keep.rownames = "NAME"))
meta <- rbind(TYPE, meta)

file2 <- paste0(out.name, "_MetaData.txt")
print(paste0("Writing out metadata...", file2))
write.table(meta,
          file = file2,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t",
          col.names = TRUE
          )

###########################################################
### 4) Extract cluster information
###########################################################
print("")
print("***************************")
print("Cluster Coordinates")
print(Sys.time())
print("***************************")

if(umap == TRUE){

  print("Extract UMAP Coordinates....")
  umap.coords <- as.data.frame(dat@reductions$umap@cell.embeddings)
  colnames(umap.coords) <- c("X", "Y")
  umap.coords <- remove.factors(umap.coords)

  classes <- as.vector(unlist(lapply(umap.coords , class)))
  TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric", "group"))
  umap.coords  <- data.frame(setDT(umap.coords , keep.rownames = "NAME"))
  umap.coords  <- rbind(TYPE, umap.coords)

  file3 <- paste0(out.name, "_UMAP_ClusterCoordinates.txt")
  print(paste0("Writing out UMAP Coordinates....", file3))

  write.table(umap.coords,
            file=file3,
            row.names = FALSE,
            quote = FALSE,
             sep = "\t",
            #col.names = TRUE
            )
}

if(tsne == TRUE){

  print("Extract tSNE Coordinates....")
  tsne.coords <- as.data.frame(dat@reductions$tsne@cell.embeddings)
  colnames(tsne.coords) <- c("X", "Y")
  tsne.coords <- remove.factors(tsne.coords)

  classes <- as.vector(unlist(lapply(tsne.coords , class)))
  TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric", "group"))
  tsne.coords  <- data.frame(setDT(tsne.coords , keep.rownames = "NAME"))
  tsne.coords  <- rbind(TYPE, tsne.coords)

  file3 <- paste0(out.name, "_tSNE_ClusterCoordinates.txt")
  print(paste0("Writing out tSNE Coordinates....", file3))

  write.table(tsne.coords,
            file=file3,
            row.names = FALSE,
            quote = FALSE,
             sep = "\t",
            #col.names = TRUE
            )
}



#######################################################
### EXAMPLE: PUSH FILES TO PORTAL VIA GOOGLE CLOUD
###
### 1) make sure you are on the interactive build parition on H4H
### salloc -c 1 -t 1:0:0 --mem 5G -p build
###
### 2) move to where the tool is located
### cd /cluster/home/lrichard/yea
###
### 3) Make sure files to trasnfer are in home dir
### cp filename /cluster/home/lrichard/Broad_scRNA_Uploads/
###
### 4) login to google cloud
### gcloud init
###
### 5) copy file to bucket # (located on your portal study page)
### gsutil cp /cluster/home/lrichard/Broad_scRNA_Uploads/SU2C_scRNA/* gs://fc-dc9e0b19-ec18-47ca-9e6c-ae13ac5e1889
###
### 6) login to Broad single cell portal wbstire and sync study
#######################################################
