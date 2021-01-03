###############################################################################
## Do Exploration of Wang 2019 data after filtering for glioma cells
library(knitr)
library(rmarkdown)
library(Seurat)
suppressPackageStartupMessages({source('whitley_scRNA_helpers.R')})
# setup
html_output_dir <- 'reports/Neftel_2019'
preproc_data_dir <- '~/projects/su2c_v2/data/preprocessed/scRNA/Neftel_2019_Filtered/'
dir(preproc_data_dir)
scRNA_file <- "Neftel_2019_scRNA_filtered.rds" 
Neftel_2019_seurat <- readRDS(file = file.path(preproc_data_dir, scRNA_file))
# calculate score difference features (needs to be recalculated)
Neftel_2019_seurat <- calc_feat_diff(Neftel_2019_seurat, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
Neftel_2019_seurat <- calc_feat_diff(Neftel_2019_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
# define features to plot
feats_plot <- colnames(Neftel_2019_seurat@meta.data)[grep('(glioma.stem.cell|Developmental|Injury_Response|cahoy_(astrocyte|oligodendrocyte|neuron))_AUC$|Diff|(tumour name)',
                                                           colnames(Neftel_2019_seurat@meta.data))]
# analysis routine
analysis_routine <- function(seurat_obj, html_dir, features, prefix) {
  
  print('Running TSNE')
  print(Sys.time())
  seurat_obj <- Seurat::RunTSNE(seurat_obj)
  # print('Running Diffusion Map')
  # seurat_obj <- Seurat::RunDiffusion(seurat_obj)
  # print(Sys.time())
  print('Finished Dim Red steps')
  print(Sys.time())
  rmarkdown::render('seurat_dr_plot_generic.Rmd',
                    output_dir = html_dir,
                    output_file = paste0(prefix, '_dim_red.html'),
                    params = list(inp_data = seurat_obj,
                                  feats_plot = features,
                                  dr_use = c('pca', 'tsne')),
                    envir = new.env()
  )
  gc(full = TRUE)
  
  rmarkdown::render('corplot_2_var.Rmd',
                    output_dir = html_dir,
                    output_file = paste0(prefix, '_corplots.html'),
                    params = list(inp_data = seurat_obj,
                                  cor_pairs = list(c('Developmental_AUC', 'Injury_Response_AUC'), 
                                                   c('glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC'), 
                                                   c('Dev_IR_Diff', 'Stem_Astro_Diff')),
                                  feats_plot = features,
                                  filter_by = 'tumour name'),
                    envir = new.env()
  )
  gc(full = TRUE)
  print('Finished making notebooks')
  print(Sys.time())
}

analysis_routine(Neftel_2019_seurat, html_output_dir, feats_plot, 'Neftel_2019_scRNA_filtered')

sessionInfo()
