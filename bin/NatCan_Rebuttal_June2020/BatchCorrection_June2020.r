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

### 3.2) re-processs each sample### we want to mimic the preporcessing to be as close to ours as possible
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
### "Next we will build the joint graph that encompasses all the samples. We do
### that by pairwise projecting samples onto a common space and establishing
### kNN of mNN pairs between the samples. We then append within-sample kNN #
### neighbours to the graph to ensure that all the cell are included in the
## graph.""
BTSC.con$buildGraph(k = 15,
                    k.self = 5,
                    space = "PCA",
                    ncomps = 30,  ### 30  for global dataset to match original processing
                    n.odgenes = 2000,
                    matching.method = "mNN",
                    metric = "angular",
                    score.component.variance = TRUE,
                    verbose = TRUE
                  )

### 3.5) Identify global clusters
BTSC.con$findCommunities(method=leiden.community, resolution=1)

### 3.6) Embed graph
BTSC.con$embedGraph()

### 3.7) Convert conos object back to seurat object
BTSC.conos <- as.Seurat(BTSC.con)

### 3.8) Plot results
DimPlot(ifnb,
        reduction = "largeVis",
        group.by = c("stim", "ident", "seurat_annotations"),
        ncol = 3
      )

#fraction of samples in each cluster
plotClusterBarplots(con, legend.height = 0.1)

### 3.9) Save results
saveRDS(BTSC.con, file = "Global_SU2C_GSCs_CONOS.rds")
saveRDS(BTSC.conos, file = "Global_SU2C_GSCs_Seurat.CONOS.rds")
saveRDS(BTSC.panel, file = "Global_SU2C_GSCs_Seurat_samplePanel.rds")



#########################################
# 4) Run LIGER
#########################################

### Ref: http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/liger.html
### https://github.com/MacoskoLab/liger
