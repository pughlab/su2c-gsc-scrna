##############################################################
#              Cluster GSCs with seurat and spectral         #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/scran_clustering


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load invididual GSC objects
### 2) Run alternate seurat clustering algorithms
### 3) Spectral clustering across 2-10 centers (k)
### 4)
### 3) Compare cluster number
### 4) Compare cluster compostion
### 5) Plot scran clusters
##############################################################
options(stringsAsFactors = F)

library("Seurat",
        lib ="/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/bin"
      ) #old version v2.3.4
library(BBmisc)
library(ggplot2)
library(gridExtra)
library(cluster)
library(kernlab)

#### write a function for UMAP plotting
plot_tSNE <- function(dat){

        plot.title <-  paste(unique(dat$orig.ident), " ", nrow(dat), "cells")

        sample_umap <- ggplot(dat, aes(x=tSNE1, y=tSNE2, color=Rphenograph_clusters)) +
                   geom_point(alpha = 0.6, size = 1.5, pch = 16) +
                   labs(x = NULL, y = NULL, title = plot.title) +
                   #scale_colour_brewer(palette = "Dark2") +
                   theme_bw() +
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                        panel.border = element_rect(linetype = "solid", fill = NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())  +
                    guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
    theme(legend.position = "none")
}


##############################################################
# 1) Run Clustering Tools (Seurat and Spectral)
##############################################################

file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/"
files <- list.files(file.path, pattern = "_L_res")
load.files <- paste0(file.path, files)

plots <- list() ##append phenograph UMAPs to
DE_scran <- list()
DE_seurat <- list()
meta_combo <- list()

for (i in 1:length(files)){

      ##############################################################
      print("")
      print("*****************")
      print(paste0(i, "/", length(files)))
      print(files[i])
      print("*****************")
      load(load.files[i]) #load GSC Data
      names(BTSC@ident) <- rownames(BTSC@meta.data)

      ##############################################################
      ###isolate parameters
      meta <- BTSC@meta.data
      meta <- cbind(meta, BTSC@dr$tsne@cell.embeddings[ ,1:2])
      meta$PC.use <- max(BTSC@calc.params$FindClusters.res.0.5$dims.use)
      colnames(meta) <- gsub("Cluster.ID", "Cluster_3_SLM", colnames(meta))

      dims <- 1:unique(meta$PC.use)
      res <- unique(meta$Opt.Resolution)
      res <- gsub("res.", "", res)
      res <- as.numeric(res)

      ##############################################################
      print("Running Seurat Original Louvain #1...")
      BTSC <- FindClusters(object = BTSC,
                           dims.use=dims,
                           reduction.type = "pca",
                           resolution=res,
                           algorithm=1,
                           force.recalc = TRUE,
                           save.SNN = TRUE,
                           print.output = FALSE
                          )
      meta$Cluster_1_Louvian <- paste0("C", as.numeric(BTSC@ident))

      ##############################################################
      print("Running Seurat Louvain multi-level refinement #2...")
      BTSC <- FindClusters(object = BTSC,
                           dims.use=dims,
                           reduction.type = "pca",
                           resolution=res,
                           algorithm=2,
                           force.recalc = TRUE,
                           save.SNN = TRUE,
                           print.output = FALSE
                          )
      meta$Cluster_2_LouvianMultiLev <- paste0("C", as.numeric(BTSC@ident))

      ##############################################################
      print("Running spectral clustering on pc matrix...")
      pc.dat <- BTSC@dr$pca@cell.embeddings[ ,dims]

      sc <- specc(pc.dat, centers=2)
      meta$Cluster_Spectral_2 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=3)
      meta$Cluster_Spectral_3 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=4)
      meta$Cluster_Spectral_4 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=5)
      meta$Cluster_Spectral_5 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=6)
      meta$Cluster_Spectral_6 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=7)
      meta$Cluster_Spectral_7 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=8)
      meta$Cluster_Spectral_8 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=9)
      meta$Cluster_Spectral_9 <- paste0("C", sc)

      sc <- specc(pc.dat, centers=10)
      meta$Cluster_Spectral_10 <- paste0("C", sc)

      ##############################################################
      print("Calculate silhouette width across methods...")

      BTSC <- SetIdent(BTSC, value = "Cluster.ID")
      seurat_clusters <- as.integer(Idents(BTSC))
      sil_1 <- cluster::silhouette(seurat_clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_seurat[[sample]] <- summary(sil_1)$clus.avg.widths


}



##############################################################
# 2) Plot and save results
##############################################################

pdf("GSC_UMAP_RphenographClustering.pdf", height = 21, width = 15)
do.call(grid.arrange, plots)
dev.off()

saveRDS(meta_combo, file = "GSC_Rphenograph_meta.rds")
saveRDS(DE_seurat, file = "GSC_Rphenograph_DE_SeuratClusters.rds")
saveRDS(DE_phenograph, file = "GSC_Rphenograph_DE_PhenographClusters.rds")


##############################################################
# 3) Calculate sil widht per cluster
##############################################################
### Ref: https://www.rdocumentation.org/packages/cluster/versions/2.1.0/topics/silhouette

options(stringsAsFactors = F)

file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/"
files <- list.files(file.path, pattern = "_L_res")
load.files <- paste0(file.path, files)

meta.file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Phenograph/"
meta.files <- list.files(meta.file.path, pattern = "phenograph_meta.rds")

sil_phenograph <- list() #averag per cluster
sil_seurat <- list() #average per cluster

for (i in 1:length(files)){

      print("")
      print("*****************")
      print(paste0(i, "/", length(files)))
      print(files[i])
      print("*****************")
      load(load.files[i]) #load GSC Data
      BTSC@meta.data$PC.use <- length(BTSC@calc.params$RunUMAP$dims.use)
      BTSC <- UpdateSeuratObject(BTSC)
      sample <- as.character(unique(BTSC@meta.data$orig.ident))
      pc.dat <- BTSC@reductions$pca@cell.embeddings[ ,1:unique(BTSC@meta.data$PC.use)] #subset out pc cell embeddings
      pc.dist <- dist(pc.dat)

      BTSC <- SetIdent(BTSC, value = "Cluster.ID")
      seurat_clusters <- as.integer(Idents(BTSC))
      sil_1 <- cluster::silhouette(seurat_clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_seurat[[sample]] <- summary(sil_1)$clus.avg.widths

      #load R phenograph ids from saved meta
      meta.load.file <- paste0(meta.file.path, meta.files[grep(sample, meta.files)])
      meta <- readRDS(meta.load.file)
      phenograph_clusters <- as.integer(meta$Rphenograph_clusters)
      sil_2 <- cluster::silhouette(phenograph_clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_phenograph[[sample]] <- summary(sil_2)$clus.avg.widths

}

saveRDS(sil_seurat, file = 'GSC_Rphenograph_Sil_SeuratClusters.rds')
saveRDS(sil_phenograph, file = 'GSC_Rphenograph_Sil_PhenographClusters.rds')
