### idea to scale GSC, Tumour and Nuclei as batches and then run PCA
###

library(Seurat)

seurat <- readRDS("GSCs_Tumour_LiveNuclei_fastMNN_Seurat.rds")
seurat@meta.data$SampleTech <- paste(seurat@meta.data$SampleType, seurat@meta.data$Technology, sep = "_")

ribo.genes <- c(rownames(seurat@assays$RNA@data)[grep("^RP[1:9]", rownames(seurat@assays$RNA@data))],
              rownames(seurat@assays$RNA@data)[grep("^RP[L,S]", rownames(seurat@assays$RNA@data))]
             ) #667 genes
pca.genes <- rownames(seurat@assays$RNA@data)[!rownames(seurat@assays$RNA@data) %in% ribo.genes]

split <- ScaleData(seurat,
                  split.by = "SampleTech",
                  do.center = TRUE,
                  vars.to.regress = c("CC.Difference", "nCount_RNA", "percent.mito"),
                  verbose = T,
                  features = rownames(seurat@assays$RNA@scale.data)
                  )
print(dim(split@assays$RNA@scale.data))
split.file <-  "Nuclei_SampleTech_SplitData_unprocessed.rds"
print("Saving split intermediate file that takes way to long to generate.....")
print(split.file)
saveRDS(split, file = split.file)

print("Running PCA...")
split <- RunPCA(split,
               features = pca.genes,
               npcs = 10,
               verbose = FALSE,
               #reduction.key = "ScaledPC_"
              )

saveRDS(split, file = "Nuclei_SampleTech_SplitData_PCA.rds")
