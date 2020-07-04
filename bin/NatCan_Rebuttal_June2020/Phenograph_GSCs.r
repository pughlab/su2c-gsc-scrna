##############################################################
#                Cluster GSCs with Phenograph                #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: https://github.com/JinmiaoChenLab/Rphenograph
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/Phenograph


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

#### write a function for UMAP plotting
plot_UMAP <- function(dat){

        plot.title <-  paste(unique(dat$orig.ident), " ", nrow(dat), "cells")

        sample_umap <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color=Rphenograph_clusters)) +
                   geom_point(alpha = 0.6, size = 1.5, pch = 16) +
                   labs(x = NULL, y = NULL, title = plot.title) +
                   scale_colour_brewer(palette = "Dark2") +
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
      BTSC@meta.data$UMAP1 <- BTSC@reductions$pca@cell.embeddings[ ,1]
      BTSC@meta.data$UMAP2 <- BTSC@reductions$pca@cell.embeddings[ ,2]

      meta <- BTSC@meta.data
      meta_combo[[sample]] <- meta
      meta.file <- paste0(sample, "_Rphenograph_meta.rds")
      saveRDS(meta, file = meta.file)
      plots[[i]] <- plot_UMAP(meta)

      #differential gene expression analysis
      BTSC <- SetIdent(BTSC, value = "Cluster.ID")
      DE_seurat[[sample]] <- FindAllMarkers(BTSC) #seurat clusters
      BTSC <- SetIdent(BTSC, value = "Rphenograph_clusters")
      DE_phenograph[[sample]] <- FindAllMarkers(BTSC) #phenograph clusters

}



##############################################################
# 2) Plot and save results
##############################################################

pdf("GSC_RphenographClustering.pdf", height = 21, width = 15)
do.call(grid.arrange, plots)
dev.off()

saveRDS(meta_combo, file = "GSC_Rphenograph_meta.rds")
saveRDS(DE_seurat, file = "GSC_Rphenograph_DE_SeuratClusters.rds")
saveRDS(DE_phenograph, file = "GSC_Rphenograph_DE_PhenographClusters.rds")
