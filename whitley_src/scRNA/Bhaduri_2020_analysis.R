###############################################################################
### Bhaduri Run Analyses
library(Seurat)
library(knitr)
library(rmarkdown)
source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')
top_dir <- '~/projects/su2c_v2'
html_output_dir <- './reports/Bhaduri_2020'
data_dir <- file.path(top_dir, 'data')
preproc_data_dir <- file.path(data_dir, 'preprocessed/scRNA/Bhaduri_2020')
inp_fname <- 'bhaduri_2020_tumor_seurat_scored.rds'
bhaduri_2020_seurat <- readRDS(file = file.path(preproc_data_dir, inp_fname))


feats_plot <- colnames(bhaduri_2020_seurat@meta.data)[grep('(glioma.stem.cell|Developmental|Injury_Response|cahoy_astrocyte)_AUC$|Diff|Study_ID',
                                                           colnames(bhaduri_2020_seurat@meta.data))]

rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = html_output_dir,
                  output_file = 'Bhaduri_2020_dim_red.html',
                  params = list(inp_data = bhaduri_2020_seurat,
                                feats_plot = feats_plot),
                  envir = new.env()
                  )
gc(full = TRUE)
rmarkdown::render('corplot_2_var.Rmd',
                  output_dir = html_output_dir,
                  output_file = 'Bhaduri_2020_corplots.html',
                  params = list(inp_data = bhaduri_2020_seurat,
                                cor_pairs = list(c('Developmental_AUC', 'Injury_Response_AUC'), 
                                                 c('glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC'), 
                                                 c('Dev_IR_Diff', 'Stem_Astro_Diff')),
                                feats_plot = feats_plot,
                                filter_by = 'Study_ID'),
                  envir = new.env()
                  )
gc(full = TRUE)