##############################################################
#                Cluster GSCs with Phenograph                #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: https://github.com/JinmiaoChenLab/Rphenograph
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/SC3


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load invididual GSC objects
### 2) Run Rphenogrpah using same number of PCs as original pipeline
### 3) Compare cluster number
### 4) Compare cluster compostion
### 5) Plot Rphenograph
##############################################################
library(Seurat) #v3.1.5
library(AUCell) #v1.8.0
library(BBmisc)
library(ggplot2)
library(gridExtra)
library(Rphenograph)
library(cluster)

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
# 1) Run Phenograph
##############################################################

options(stringsAsFactors = F)

file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/"
files <- list.files(file.path, pattern = "_L_res")
load.files <- paste0(file.path, files)

plots <- list() ##append phenograph UMAPs to
DE_phenograph <- list()
DE_seurat <- list()
meta_combo <- list()

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

      #Run phenopgraph
      Rphenograph_out <- Rphenograph(pc.dat, k = 45)
      Rphenograph_clusters <- membership(Rphenograph_out[[2]])
      names(Rphenograph_clusters) <- rownames(BTSC@meta.data)
      BTSC@meta.data$Rphenograph_clusters <- as.factor(Rphenograph_clusters)
      BTSC@meta.data$tSNE1 <- BTSC@reductions$tsne@cell.embeddings[ ,1]
      BTSC@meta.data$tSNE2 <- BTSC@reductions$tsne@cell.embeddings[ ,2]

      meta <- BTSC@meta.data
      meta_combo[[sample]] <- meta
      meta.file <- paste0(sample, "_Rphenograph_meta.rds")
      saveRDS(meta, file = meta.file)
      plots[[i]] <- plot_tSNE(meta)

      #differential gene expression analysis
      BTSC <- SetIdent(BTSC, value = "Cluster.ID")
      DE_seurat[[sample]] <- FindAllMarkers(BTSC) #seurat clusters
      BTSC <- SetIdent(BTSC, value = "Rphenograph_clusters")
      DE_phenograph[[sample]] <- FindAllMarkers(BTSC) #phenograph clusters

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