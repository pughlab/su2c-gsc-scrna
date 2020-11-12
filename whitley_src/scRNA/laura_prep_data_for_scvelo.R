###############################################################################
library(knitr)
library(rmarkdown)
source("laura_SU2C_GBM_GSCs_combined_explore_helpers.R")
diffusion_map_data_dir <- '~/projects/su2c_v2/data/raw/scRNA/Data/laura_Diffusion_Map_July_2020/'
tumour_gsc_dir <- file.path(diffusion_map_data_dir, 'Tumour_GSCs')
gsc_dir <- file.path(diffusion_map_data_dir, 'GSCs')
dir(gsc_dir)
GSC_Tumour_seurat <- readRDS(file.path(tumour_gsc_dir, "GSC_Tumour_500cells_20343cells_DM_Seurat.rds"))
GSC_Tumour_seurat <- calc_feat_diff(GSC_Tumour_seurat, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff', scale = TRUE)
GSC_Tumour_seurat <- calc_feat_diff(GSC_Tumour_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff', scale = TRUE)
GSC_Tumour_seurat@meta.data$Dev_IR_Diff_zscore <- scale(GSC_Tumour_seurat@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
GSC_Tumour_seurat@meta.data$Stem_Astro_Diff_zscore <- scale(GSC_Tumour_seurat@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)

scvelo_plot_data <- Seurat::FetchData(GSC_Tumour_seurat, vars.all = c(colnames(GSC_Tumour_seurat@meta.data), 'DM1', 'DM2'))

scvelo_plot_data$unique_id <- sub('(BTSC|GBM)_', '', rownames(scvelo_plot_data))

output_dir <- '~/projects/su2c_v2/data/preprocessed/scRNA/velocyto_GBM_GSCs_split_by_sample_no_G800_L'
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(scvelo_plot_data, file = file.path(output_dir, 'scvelo_plot_data_GBM_GSCs.csv'))

rm(GSC_Tumour_seurat)
gc(full = TRUE)

GSC_seurat <- readRDS(file.path(gsc_dir, "GSC_500cells_14000cells_DM_Seurat.rds"   ))
GSC_seurat <- calc_feat_diff(GSC_seurat, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff', scale = TRUE)
GSC_seurat <- calc_feat_diff(GSC_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff', scale = TRUE)
GSC_seurat@meta.data$Dev_IR_Diff_zscore <- scale(GSC_seurat@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
GSC_seurat@meta.data$Stem_Astro_Diff_zscore <- scale(GSC_seurat@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)

scvelo_plot_data <- Seurat::FetchData(GSC_seurat, vars.all = c(colnames(GSC_seurat@meta.data), 'DM1', 'DM2'))

scvelo_plot_data$unique_id <- rownames(scvelo_plot_data)

write.csv(scvelo_plot_data, file = file.path(output_dir, 'scvelo_plot_data_GSCs.csv'))


BTSC_GBM_full_data_dir <- file.path('~/projects/su2c_v2/data/raw/scRNA/Data/SU2C_BTSC_TumourCells_combined_Oct_2019')
BTSC_GBM_file <- 'BTSC_TumourCells_G800Removed_AUCell_Seurat.Rdata'
load(file.path(BTSC_GBM_full_data_dir, BTSC_GBM_file))
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff', scale = TRUE)
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff', scale = TRUE)
BTSC_TumourCells@meta.data$Dev_IR_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
BTSC_TumourCells@meta.data$Stem_Astro_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)

scvelo_plot_data <- Seurat::FetchData(BTSC_TumourCells, vars.all = c(colnames(BTSC_TumourCells@meta.data), 'PC1', 'PC2'))

scvelo_plot_data$unique_id <- sub('(BTSC|GBM)_', '', rownames(scvelo_plot_data))

write.csv(scvelo_plot_data, file = file.path(output_dir, 'scvelo_plot_data_GBM_GSCs_PCA_full.csv'))
