###########################################################
#             Run Batch Correction Tools on GSCs          #
#                         L.Richards                      #
#                         June 2020                       #
###########################################################
### Reference: https://github.com/satijalab/seurat-wrappers
### Working dir: /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/BatchCorrection/

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Load Global GSC data (corresponds to Ext Data Fig 4.)
### 2) Run CONOS wrapper in Seurat
### 3) Run LIGER in Seurat
### 4) Run fastMNN in Seurat
###########################################################


#########################################
# 1) Install + Load packages
#########################################
#library(devtools)
#devtools::install_github("RcppCore/Rcpp")
#devtools::install_github("satijalab/seurat-wrappers")
#devtools::install_github("hms-dbmi/conos")
#devtools::install_github('MacoskoLab/liger')

library(Seurat) #v3.1.5
library(SeuratWrappers) #v0.2.0
library(conos) #v1.3.0
library(liger)

#########################################
# 2) Load GSC object
#########################################

### file copied from /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata to H4H
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
table(BTSC@meta.data$SampleID) #29 samples, 69393 cells
table(BTSC@meta.data$res.2) #61 clusters
### convert object to seurat v3 to be compatible with wrappers
### used 30 PCs in original
BTSC <- UpdateSeuratObject(BTSC)

### Plot default UMAP with current paper parameters
### 2x2 grid where you color by cluster ID and
pdf("OriginalClustering_res.2.pdf", height = 7, width = 18)

DimPlot(BTSC,
        group.by = c("SampleID", "res.2"),
        reduction = "umap",
        pt.size = 0.1,
        label = T,
        ncol = 2
      ) + NoLegend()

DimPlot(BTSC,
        group.by = c("SampleID", "res.2"),
        reduction = "tsne",
        pt.size = 0.1,
        label = T,
        ncol = 2
        ) + NoLegend()

dev.off()

### Plot default tSNE with current paper parameters


#########################################
# 3) Run CONOS
#########################################

### Ref: http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/conos.html
### https://github.com/hms-dbmi/conos

### 3.1) split up by sample
BTSC.panel <- SplitObject(BTSC, split.by = "SampleID")

### 3.2) re-processs each sample
### we want to mimic the preporcessing to be as close to ours as possible
### run PCA on all genes except ribosomal
### regress out CC.difference, nUMI and percent.mito
for (i in 1:length(BTSC.panel)){

  print(names(BTSC.panel[i]))

  ribo.genes <- c(rownames(BTSC.panel[[i]]@assays$RNA@data)[grep("^RP[1:9]", rownames(BTSC.panel[[i]]@assays$RNA@data))],
                rownames(BTSC.panel[[i]]@assays$RNA@data)[grep("^RP[L,S]", rownames(BTSC.panel[[i]]@assays$RNA@data))]
               ) #667 genes

  pca.genes <- rownames(BTSC.panel[[i]]@assays$RNA@data)[!rownames(BTSC.panel[[i]]@assays$RNA@data) %in% ribo.genes]
  print(length(pca.genes))

  BTSC.panel[[i]] <- NormalizeData(BTSC.panel[[i]]) %>% ScaleData(vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"), verbose = T) %>%
        RunPCA(pc.genes = pca.genes, verbose = FALSE)

}

### 3.3) Create conos object
### This step also clusters each sample in its own
BTSC.con <- Conos$new(BTSC.panel)

### 3.4) Build Graph with PCA space
### There are other options but PCA is default
### This step is very slow.....
### "Next we will build the joint graph that encompasses all the samples. We do
### that by pairwise projecting samples onto a common space and establishing
### kNN of mNN pairs between the samples. We then append within-sample kNN #
### neighbours to the graph to ensure that all the cell are included in the
## graph.""
BTSC.con$buildGraph(k = 15,
                    k.self = 5,
                    space = "PCA",
                    ncomps = 30,  ### 30  for global dataset to match original processing
                    n.odgenes = 2000, ### over dispersed genes
                    matching.method = "mNN",
                    metric = "angular",
                    score.component.variance = TRUE,
                    verbose = TRUE
                  )

### 3.5) Identify global clusters
### Can also cluster with different methods
### use default resolution of 1
BTSC.con$findCommunities(method=leiden.community, resolution=1)

### 3.6) Embed graph
### options are LargeVis or umap
### LargeVis is default, but use umap to be consistent with original processing
### use same min.dist and spread as Seurat defaults
### both seurat and CONOS implement UMAP with uwot package
BTSC.con$embedGraph(method = "UMAP", min.dist = 0.3, spread = 1)

### 3.7) Convert conos object back to seurat object
BTSC.conos <- as.Seurat(BTSC.con)
#ident is conos clusters

### 3.8) Plot results
pdf("CONOS_GSC_umaps.pdf", height = 7, width = 18)
DimPlot(BTSC.conos,
        reduction = "UMAP",
        group.by = c("SampleID", "ident"),
        ncol = 2
      )
dev.off()

### fraction of samples in each cluster
pdf("CONOS_GSC_fractionSamples.pdf", height = 18, width = 18)
plotClusterBarplots(BTSC.con, legend.height = 0.1)
dev.off()

### 3.9) Save results
saveRDS(BTSC.con, file = "Global_SU2C_GSCs_CONOS.rds")
saveRDS(BTSC.conos, file = "Global_SU2C_GSCs_Seurat.CONOS.rds")
saveRDS(BTSC.panel, file = "Global_SU2C_GSCs_Seurat_samplePanel.rds")



#########################################
# 4) Run LIGER
#########################################
### Ref: http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/liger.html
### https://github.com/MacoskoLab/liger

library(liger)

### load GSC data and update to seurat v3 object
load("/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")
BTSC <- UpdateSeuratObject(BTSC)

### 4.1) Scale data within each sample
### Regress out sample factors as original analysis
BTSC.liger <- ScaleData(BTSC,
                        split.by = "SampleID",
                        do.center = FALSE,
                        vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                         verbose = T
                      )

### 4.2) Joint Matrix Factorization (very time consuming)
### lambda = larger values penalize dataset-specific effects more strongly (ie.
### alignment should increase as lambda increases). Default is 5.
BTSC.liger <- RunOptimizeALS(BTSC.liger,
                              k = 20,
                              lambda = 5, #default value, effects alignemnt
                              split.by = "SampleID"
                            )
saveRDS(BTSC.liger, file = "Global_SU2C_GSCs_Seurat_LIGER.rds")
### 4.3)  Quantile Normalization and Joint Clustering
### This step concerns me because of the ref_dataset paramter use in cancer
### "Name of dataset to use as a "reference" for normalization.By default, the dataset with the largest number of cells is used."
### It does not seem appropiate to use the largest GSC dataset as the refernece ### since we know they have different transcriptomic and genetic properties
BTSC.liger <- RunQuantileNorm(BTSC.liger,
                              split.by = "SampleID",
                              ref_datset = NULL
                            )

### 4.4) FindNeighbours using NMF from last step
BTSC.liger <- FindNeighbors(BTSC.liger,
                            reduction = "iNMF",
                            dims = 1:30
                          )
### 4.5) Find Clusters
BTSC.liger <- FindClusters(BTSC.liger,
                            resolution = 2 #use same resolution to match original protocol
                          )

### 4.6) Dimensional reduction wiht UMAP
BTSC.liger <- RunUMAP(BTSC.liger,
                      dims = 1:ncol(BTSC.liger[["iNMF"]]),
                      reduction = "iNMF"
                    )

### 4.7) Plotting
DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)

### 4.8) Save Data
saveRDS(BTSC.liger, file = "Global_SU2C_GSCs_Seurat_LIGER.rds")




#########################################
# 5) Run fastMNN
#########################################

### Ref: http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html
### https://bioconductor.org/packages/release/bioc/html/batchelor.html




#########################################
# 6) Run Batchelor package options
#########################################

### Ref: http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html
### https://bioconductor.org/packages/release/bioc/html/batchelor.html


##########################
# 6.1) fastMNN
##########################
### already run above with seurat wrapper

##########################
# 6.2) Batch rescaling
##########################
### While this method is fast and simple, it makes the strong assumption that
### the population composition of each batch is the same. This is usually not
### the case for scRNA-seq experiments in real systems that exhibit biological
### variation. Thus, rescaleBatches() is best suited for merging technical
### replicates of the same sample, e.g., that have been sequenced separately.




##########################
# 6.3) Multi-batch Normalization
##########################
### Differences in sequencing depth between batches are an obvious cause for
### batch-to-batch differences. These can be removed by multiBatchNorm(), which
### downscales all batches to match the coverage of the least-sequenced batch.
### This function returns a list of SingleCellExperiment objects with log-
### transformed normalized expression values that can be directly used for
### further correction.


##########################
# 6.4) Multi-batch PCA
##########################
