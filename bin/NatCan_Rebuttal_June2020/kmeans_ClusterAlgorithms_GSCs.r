##############################################################
#              Cluster GSCs with  kmeans                     #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/kmeans_leiden_clustering

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load invididual GSC objects
### 2) run TSCAN to cluster
### 4) Calculate sil width
### 5) Plot tSNEs
### 6) Save Data
##############################################################

##############################################################
### EXAMPLE EXECUTION ON H4H
### #!/bin/bash
### #SBATCH --mem=60G
### #SBATCH -p himem
### #SBATCH -c 30
### #SBATCH -N 1
### #SBATCH --account=pughlab
###
### module load R/3.6.1
###
### Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/kmeans_clustering_GSCs.r
###
##############################################################
options(stringsAsFactors = F)

library(Seurat)
library(BBmisc)
library(ggplot2)
library(gridExtra)
library(cluster)
library(kernlab)
library(optparse)
library(factoextra)
library(NbClust)


##############################################################
# 1) Run Clustering Tools (Seurat and Spectral)
##############################################################

file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/"
files <- list.files(file.path, pattern = "_L_res")
load.files <- paste0(file.path, files)

meta_combo <- list()
sil_list <- list()


for (i in 1:length(files)){
##############################################################
print("")
print("*****************")
print(files[i])
print("*****************")
load(load.files[i]) #load GSC Data
names(BTSC@ident) <- rownames(BTSC@meta.data)

#############################################################
### isolate parameters
meta <- BTSC@meta.data
meta <- cbind(meta, BTSC@dr$tsne@cell.embeddings[ ,1:2])
meta$PC.use <- max(BTSC@calc.params$FindClusters.res.0.5$dims.use)
colnames(meta) <- gsub("Cluster.ID", "Cluster_3_SLM", colnames(meta))

dims <- 1:unique(meta$PC.use)
res <- unique(meta$Opt.Resolution)
res <- gsub("res.", "", res)
res <- as.numeric(res)

sample <- as.character(unique(meta$orig.ident))

##############################################################
print("Determining optimal k.....")
pc.dat <- BTSC@dr$pca@cell.embeddings[ ,dims]
a <- NbClust(data = pc.dat,  distance = "euclidean",
        min.nc = 2, max.nc = 10, method = "kmeans")
k.file <- paste0(sample, "kmeans_optimalCluster.pdf")
pdf(k.file)
fviz_nbclust(a)
dev.off()

print("Running k-means.....")
print(table(a$Best.partition))
col.name <- paste0("kmeans_", length(table(a$Best.partition)))
meta[ ,col.name] <- a$Best.partition

print("Calculate sil width....")
pc.dist <- dist(pc.dat)
clusters <- as.integer(factor(meta[ ,col.name]))
sil <- cluster::silhouette(clusters,
                             dist = pc.dist,
                             do.clus.stat = TRUE
                           )
sil_list[[sample]] <- summary(sil)$clus.avg.widths

##############################################################
meta_combo[[sample]] <- meta

}

print("Saving data")
saveRDS(meta_combo, file = "kmeans_metadata.rds")
saveRDS(sil_list, file = "kmeans_silwidths.rds")
