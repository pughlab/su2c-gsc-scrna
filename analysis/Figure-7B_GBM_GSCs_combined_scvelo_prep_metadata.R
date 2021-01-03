library(velocyto.R)
library(Seurat)

## LOAD DATA

top_dir <- '/home/owenwhitley/projects/su2c_v2'
data_dir <- file.path(top_dir, 'data')
output_dir <- file.path(top_dir, 'results/scRNA/scvelo_GSCs_Tumour_combined_no_G800_L')
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# contains Tumour + GSC data, labeled with classification as tumour or line
GSC_GBM_seurat_dir <- file.path(top_dir, 'results/scRNA/laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic')
dir(GSC_GBM_seurat_dir)
seurat_obj <- readRDS(file = file.path(GSC_GBM_seurat_dir, "logistic_regression_T_v_L_seurat.rds"))



## Data Cleaning

# clean up cell identifiers by just extracting barcode. do for sce object and loom data.
# from sce object, remove any barcodes that were present in more than 1 sample
# (we've come across at least 1 of such edge cases)
cell_name_regex <- '(G|BT)[0-9]+-*[A-Z]*_[TL]_[ATCG]*$'
cell_names <- rownames(seurat_obj@meta.data)
sample_barcode_id <- regmatches(rownames(seurat_obj@meta.data), regexpr(cell_name_regex, cell_names))
dup_barcodes <- sample_barcode_id[which(duplicated(sample_barcode_id))]
valid_columns <-  !(vapply(cell_names, FUN.VALUE = logical(1), FUN = function(x) {
  regmatches(x, regexpr('[TCAG]*$', x)) %in% dup_barcodes
}))
if (!all(valid_columns)) {
  print(paste('subsetting', sum(valid_columns), 'of', length(valid_columns), 'cells'))
  seurat_obj <- Seurat::SubsetData(seurat_obj, cells.use = cell_names[valid_columns], subset.raw = TRUE)
  
}
# rename seurat object's cells to barcode names. Get all data to be used, as RenameCells runs into memory issues
pca_coords <- Seurat::GetDimReduction(seurat_obj, reduction.type = 'pca', slot = 'cell.embeddings')
tsne_coords <- Seurat::GetDimReduction(seurat_obj, reduction.type = 'tsne', slot = 'cell.embeddings')
meta_data <- seurat_obj@meta.data

all(rownames(pca_coords) == rownames(meta_data)) && all(rownames(tsne_coords) == rownames(meta_data))

gene_names <- rownames(seurat_obj@data)

rm(seurat_obj)
gc(full = TRUE)

new_cols_seurat <- regmatches(rownames(meta_data), regexpr(cell_name_regex, rownames(meta_data)))
stopifnot(!any(duplicated(new_cols_seurat)))



rownames(pca_coords) <- rownames(tsne_coords) <- rownames(meta_data) <- new_cols_seurat

write.csv(meta_data, file = file.path(output_dir, 'seurat_meta_data.csv'))
write.csv(pca_coords, file = file.path(output_dir, 'seurat_pca_coords.csv'))
write.csv(tsne_coords, file = file.path(output_dir, 'seurat_tsne_coords.csv'))
write.csv(gene_names, file = file.path(output_dir, 'seurat_gene_names.csv'))

Sys.time()
sessionInfo()
