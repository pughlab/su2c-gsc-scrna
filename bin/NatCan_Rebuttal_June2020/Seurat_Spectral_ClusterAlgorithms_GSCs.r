##############################################################
#              Cluster GSCs with seurat and spectral         #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### REF: http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html
### /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/spectral_clustering

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load invididual GSC objects
### 2) Run alternate seurat clustering algorithms
### 3) Spectral clustering across 2-10 centers (k)
### 4) Calculate sil width arcross methods
### 5) Plot tSNEs
### 6) Save Data
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
plot_tSNE <- function(dat, clustername){

        plot.title <-  paste(unique(dat$orig.ident), " ", nrow(dat), "cells")

        sample_umap <- ggplot(dat, aes_string(x="tSNE_1", y="tSNE_2", color=clustername)) +
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
    return(sample_umap)
}


##############################################################
# 1) Run Clustering Tools (Seurat and Spectral)
##############################################################

file.path <- "/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/"
files <- list.files(file.path, pattern = "_L_res")
load.files <- paste0(file.path, files)

meta_combo <- list()

plots_Louvain <- list()
plots_Louvain_MLR <- list()
plots_Spectral2 <- list()
plots_Spectral3 <- list()
plots_Spectral4 <- list()
plots_Spectral5 <- list()
plots_Spectral6 <- list()
plots_Spectral7 <- list()
plots_Spectral8 <- list()
plots_Spectral9 <- list()
plots_Spectral10 <- list()

sil_Seurat <- list()
sil_Louvain <- list()
sil_Louvain_MLR <- list()
sil_Spectral2 <- list()
sil_Spectral3 <- list()
sil_Spectral4 <- list()
sil_Spectral5 <- list()
sil_Spectral6 <- list()
sil_Spectral7 <- list()
sil_Spectral8 <- list()
sil_Spectral9 <- list()
sil_Spectral10 <- list()

#for (i in 1:length(files)){
for (i in 1:2){
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

      sample <- as.character(unique(meta$orig.ident))
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
      plots_Louvain[[sample]] <- plot_tSNE(meta, "Cluster_1_Louvian")

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
      plots_Louvain_MLR[[sample]] <- plot_tSNE(meta, "Cluster_2_LouvianMultiLev")

      ##############################################################
      print("Running spectral clustering on pc matrix...")
      pc.dat <- BTSC@dr$pca@cell.embeddings[ ,dims]

      sc <- specc(pc.dat, centers=2)
      meta$Cluster_Spectral_2 <- paste0("C", sc)
      plots_Spectral2[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_2")

      sc <- specc(pc.dat, centers=3)
      meta$Cluster_Spectral_3 <- paste0("C", sc)
      plots_Spectral3[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_3")

      sc <- specc(pc.dat, centers=4)
      meta$Cluster_Spectral_4 <- paste0("C", sc)
      plots_Spectral4[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_4")

      sc <- specc(pc.dat, centers=5)
      meta$Cluster_Spectral_5 <- paste0("C", sc)
      plots_Spectral5[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_5")

      sc <- specc(pc.dat, centers=6)
      meta$Cluster_Spectral_6 <- paste0("C", sc)
      plots_Spectral6[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_6")

      sc <- specc(pc.dat, centers=7)
      meta$Cluster_Spectral_7 <- paste0("C", sc)
      plots_Spectral7[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_7")

      sc <- specc(pc.dat, centers=8)
      meta$Cluster_Spectral_8 <- paste0("C", sc)
      plots_Spectral8[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_8")

      sc <- specc(pc.dat, centers=9)
      meta$Cluster_Spectral_9 <- paste0("C", sc)
      plots_Spectral9[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_9")

      sc <- specc(pc.dat, centers=10)
      meta$Cluster_Spectral_10 <- paste0("C", sc)
      plots_Spectral10[[sample]] <- plot_tSNE(meta, "Cluster_Spectral_10")

      ######
      meta_combo[[sample]] <- meta

      ##############################################################
      print("Calculate silhouette width across methods...")

      pc.dist <- dist(pc.dat)

      ###############
      ### seurat
      clusters <- as.integer(factor(meta$Cluster_3_SLM))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Seurat[[sample]] <- summary(sil)$clus.avg.widths

      ###############
      ### Louvain
      clusters <- as.integer(factor(meta$Cluster_1_Louvian))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Louvain[[sample]] <- summary(sil)$clus.avg.widths

      ###############
      ### LouvainMLR
      clusters <- as.integer(factor(meta$Cluster_2_LouvianMultiLev))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Louvain_MLR[[sample]] <- summary(sil)$clus.avg.widths

      ###############
      ### Spectral

      clusters <- as.integer(factor(meta$Cluster_Spectral_2))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral2[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_3))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral3[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_4))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral4[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_5))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral5[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_6))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral6[[sample]] <- summary(sil)$clus.avg.widths


      clusters <- as.integer(factor(meta$Cluster_Spectral_7))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral7[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_8))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral8[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_9))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral9[[sample]] <- summary(sil)$clus.avg.widths

      clusters <- as.integer(factor(meta$Cluster_Spectral_10))
      sil <- cluster::silhouette(clusters,
                                   dist = pc.dist,
                                   do.clus.stat = TRUE
                                 )
      sil_Spectral10[[sample]] <- summary(sil)$clus.avg.widths


}



##############################################################
# 2) Plot tSNEs
##############################################################
print("*****************")
print("Plotting tSNEs across clustering methods....")
print("*****************")

pdf("GSC_tSNE_Louvain.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Louvain)
dev.off()

pdf("GSC_tSNE_LouvainMLR.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Louvain_MLR)
dev.off()

pdf("GSC_tSNE_Spectral2.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral2)
dev.off()

pdf("GSC_tSNE_Spectral3.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral3)
dev.off()

pdf("GSC_tSNE_Spectral4.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral4)
dev.off()

pdf("GSC_tSNE_Spectral5.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral5)
dev.off()

pdf("GSC_tSNE_Spectral6.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral6)
dev.off()

pdf("GSC_tSNE_Spectral7.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral7)
dev.off()

pdf("GSC_tSNE_Spectral8.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral8)
dev.off()

pdf("GSC_tSNE_Spectral9.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral9)
dev.off()

pdf("GSC_tSNE_Spectral10.pdf", height = 21, width = 15)
do.call(grid.arrange, plots_Spectral10)
dev.off()


##############################################################
# 3) Save results
##############################################################
print("*****************")
print("Saving metadata....")
print("*****************")
meta_2 <- do.call(rbind, meta_combo)
saveRDS(meta_2, file = "GSC_Louvain_Spectral_meta.rds")


print("*****************")
print("Saving sil widths....")
print("*****************")
saveRDS(sil_Seurat, file = 'GSC_Sil_SLM.rds')
saveRDS(sil_Louvain, file = 'GSC_Sil_Louvain.rds')
saveRDS(sil_LouvainMLR, file = 'GSC_Sil_LouvainMLR.rds')
saveRDS(sil_Spectral2, file = 'GSC_Sil_Spectral2.rds')
saveRDS(sil_Spectral3, file = 'GSC_Sil_Spectral3.rds')
saveRDS(sil_Spectral4, file = 'GSC_Sil_Spectral4.rds')
saveRDS(sil_Spectral5, file = 'GSC_Sil_Spectral5.rds')
saveRDS(sil_Spectral6, file = 'GSC_Sil_Spectral6.rds')
saveRDS(sil_Spectral7, file = 'GSC_Sil_Spectral7.rds')
saveRDS(sil_Spectral8, file = 'GSC_Sil_Spectral8.rds')
saveRDS(sil_Spectral9, file = 'GSC_Sil_Spectral9.rds')
saveRDS(sil_Spectral10, file = 'GSC_Sil_Spectral10.rds')
