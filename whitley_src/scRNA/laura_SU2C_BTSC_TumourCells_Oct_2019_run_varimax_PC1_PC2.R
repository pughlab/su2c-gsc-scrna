###############################################################################
## Add Varimax Rotation to Nov 2019 GSC + Tumour combined data

library(stats)
library(matrixStats)

preproc_dir <- "~/projects/su2c_v2/data/preprocessed/scRNA/GSCs_Tumour_combined"
inp_file <- 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L.rds'
output_dir <- '~/projects/su2c_v2/data/preprocessed/scRNA/GSCs_Tumour_combined_w_varimax_PC1_PC2'
output_file <- 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L_w_varimax_PC1_PC2.rds'
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

seurat_obj <- readRDS(file = file.path(preproc_dir, inp_file))

gene_loadings <- seurat_obj@dr$pca@gene.loadings[,c('PC1', 'PC2')]

varimax_rot <- varimax(gene_loadings)
gc(full = TRUE)
varimax_loadings <- varimax_rot$loadings
varimax_loadings <- varimax_loadings[1:nrow(varimax_loadings), 1:ncol(varimax_loadings)]
colnames(varimax_loadings) <- paste0('VM', 1:ncol(varimax_loadings))
gc(full = TRUE)
data_mat <- t(seurat_obj@scale.data[rownames(varimax_loadings), ])
rm(seurat_obj)
gc(full = TRUE)
new_embeddings <- data_mat %*% varimax_loadings
rm(data_mat)
gc(full = TRUE)
seurat_obj <- readRDS(file = file.path(preproc_dir, inp_file))
dim_sdev <- sqrt(matrixStats::colVars(new_embeddings))

rownames(new_embeddings) <- colnames(seurat_obj@scale.data)
colnames(new_embeddings) <- colnames(varimax_loadings)

# a quick hack to modify dimensionality reduction object and add back to seurat object
dim_red_obj <- seurat_obj@dr$pca
dim_red_obj@gene.loadings <- varimax_loadings
dim_red_obj@cell.embeddings <- new_embeddings
dim_red_obj@sdev <- dim_sdev
dim_red_obj@key <- "VM"
dim_red_obj@jackstraw <- NULL
dim_red_obj@misc <- NULL
seurat_obj@dr$varimax <- dim_red_obj

saveRDS(seurat_obj, file = file.path(output_dir, output_file))
sessionInfo()
