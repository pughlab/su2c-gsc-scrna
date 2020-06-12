###########################################################
#            Generate files for upload to Broad           #
#                        scRNA portal                     #
#                         L.Richards                      #
#                         April 2020                      #
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

#load old version of seurat
suppressMessages(library("Seurat",
                         lib ="/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/bin/")
                )
library("optparse")
library("data.table")
library("dplyr")
library("taRifx")

###########################################################
### USER DEFINE VARIABLES
###########################################################

#path to seurat object to extract data from
seurat.obj <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Hypoxia/AllTumour_AUCell_Seurat.Rdata"

#append string to output files
out.name <- "SU2C_scRNA_TumourOnly"

#path to text file with names of metadata columns to include
#meta.columns <-



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


###########################################################
### 2) Extract normalized expression
###########################################################
print("")
print("***************************")
print("Normalized expression matrix")
print(Sys.time())
print("***************************")

print("Size of expression matrix....")
print(dim(dat@data))

print("Extracting normalized expression matrix...")
norm.exp <- data.frame(as.matrix(dat@data))
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

meta$Patient.ID <- gsub("_L", "", meta$orig.ident)
meta$Sample.Type <- ifelse(grep("_L", meta$orig.ident), "GSC", "TUMOUR")


## define columns to keep in
include <- c("orig.ident",
             "Patient.ID",
             "Sample.Type",
             "nGene",
             "nUMI",
             "percent.mito",
             "Phase",
             "G2M.Score",
             "Cluster.ID",
             "Opt.Resolution"
            )
meta <- meta[, include]
colnames(meta)[c(1,10)] <- c("Sample.ID", "Optimal.Clustering.Resolution")


classes <- as.vector(unlist(lapply(meta, class)))
TYPE <- c("TYPE", ifelse(classes == "numeric", "numeric", "group"))
meta <- data.frame(setDT(meta, keep.rownames = "NAME"))
meta <- rbind(TYPE, meta)

file2 <- paste0(sample, "_MetaData.txt")
print(file2)

write.table(meta,
          file=file2,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t",
          col.names = TRUE
          )
