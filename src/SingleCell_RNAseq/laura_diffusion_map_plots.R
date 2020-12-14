###############################################################################
library(knitr)
library(rmarkdown)
source("laura_SU2C_GBM_GSCs_combined_explore_helpers.R")
data_dir <- '~/projects/su2c_v2/data/raw/scRNA/Data/laura_Diffusion_Map_July_2020/'
tumour_gsc_dir <- file.path(data_dir, 'Tumour_GSCs')
gsc_dir <- file.path(data_dir, 'GSCs')
dir(gsc_dir)
GSC_Tumour_seurat <- readRDS(file.path(tumour_gsc_dir, "GSC_Tumour_500cells_20343cells_DM_Seurat.rds"))
GSC_Tumour_seurat <- calc_feat_diff(GSC_Tumour_seurat, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff')
GSC_Tumour_seurat <- calc_feat_diff(GSC_Tumour_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
GSC_Tumour_seurat@meta.data$Dev_IR_Diff_zscore <- scale(GSC_Tumour_seurat@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
GSC_Tumour_seurat@meta.data$Stem_Astro_Diff_zscore <- scale(GSC_Tumour_seurat@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)
reports_dir <- '~/projects/su2c_v2/scripts/scRNA/reports/laura_rebuttal_July_2020_Diffusion_Map'

feat_regex <- 'Neftel|orig.ident|cahoy_(astrocyte|oligodendrocyte)|zscore'
feats_plot <- colnames(GSC_Tumour_seurat@meta.data)[grep(feat_regex, colnames(GSC_Tumour_seurat@meta.data))]
feats_plot <- c(feats_plot, 'PC1', 'PC2')

rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_DM.html',
                  params = list(inp_data = GSC_Tumour_seurat,
                                feats_plot = feats_plot,
                                dr_use = 'dm',
                                dims_plot = list(dm = list(c(1, 2)))
                                ),
                  )

rm(GSC_Tumour_seurat)

gc(full = TRUE)

GSC_seurat <- readRDS(file.path(gsc_dir, "GSC_500cells_14000cells_DM_Seurat.rds"))
GSC_seurat <- calc_feat_diff(GSC_seurat, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff')
GSC_seurat <- calc_feat_diff(GSC_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
GSC_seurat@meta.data$Dev_IR_Diff_zscore <- scale(GSC_seurat@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
GSC_seurat@meta.data$Stem_Astro_Diff_zscore <- scale(GSC_seurat@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_DM.html',
                  params = list(inp_data = GSC_seurat,
                                feats_plot = feats_plot,
                                dr_use = 'dm',
                                dims_plot = list(dm = list(c(1, 2)))
                  ),
)

sessionInfo()
