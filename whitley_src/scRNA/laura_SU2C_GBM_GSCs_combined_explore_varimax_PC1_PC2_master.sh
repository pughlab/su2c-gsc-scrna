#!/bin/bash

# Overview: run scRNA exploratory analyses for GBM + GSCs with varimax rotation on PC1, PC2

# Run Varimax Rotation

Rscript laura_SU2C_BTSC_TumourCells_Oct_2019_run_varimax_PC1_PC2.R > \
  laura_SU2C_BTSC_TumourCells_Oct_2019_run_varimax_PC1_PC2_out.txt \
  2>laura_SU2C_BTSC_TumourCells_Oct_2019_run_varimax_PC1_PC2_stderr.txt


# setup preprocessed data dir + html output directory
preproc_dir="~/projects/su2c_v2/data/preprocessed/scRNA/GSCs_Tumour_combined_w_varimax_PC1_PC2"
master_html_output_dir="./reports/laura_SU2C_GSCs_Tumour_combined_varimax_PC1_PC2"

set -e

# Run Exploration of Reduced Dims

Rscript -e "suppressPackageStartupMessages({source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')});\
            seurat_obj <- readRDS(file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L_w_varimax_PC1_PC2.rds'));\
            seurat_obj <- calc_feat_diff(seurat_obj, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff');\
            seurat_obj <- calc_feat_diff(seurat_obj, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff');\
            seurat_obj@meta.data[['Dev_IR_Diff_zscore']] <- scale(seurat_obj@meta.data[['Dev_IR_Diff']], center = TRUE, scale = TRUE);\
            seurat_obj@meta.data[['Stem_Astro_Diff_zscore']] <- scale(seurat_obj@meta.data[['Stem_Astro_Diff']], center = TRUE, scale = TRUE);\
            auc_feats <- colnames(seurat_obj@meta.data)[grep('(glioma.stem.cell|RNA.GSC.c[1-2]|cahoy_astrocyte|Neftel)', colnames(seurat_obj@meta.data))];\
            auc_feats <- auc_feats[grep('AUC$', auc_feats)];\
            rmarkdown::render('./seurat_dr_plot_generic.Rmd',\
                              output_file = 'GSC_Tumour_combined_w_varimax_explore_dims.html',\
                              output_dir = '${master_html_output_dir}',\
                              params = list(inp_data = seurat_obj,\
                                            dr_use = 'varimax',\
                                            dims_plot = list(varimax = list(c(1,2))),\
                                            feats_plot = c('Dev_IR_Diff_zscore', 'Dev_IR_Diff', 'Stem_Astro_Diff_zscore', 'Stem_Astro_Diff', auc_feats, 'SampleType'),\
                                            pt_size = 0.5,\
                                            split_by = 'SampleType',\
                                            grid_dist = 5.0)\
                              )"

# Pathway Analysis on reduced dims                              
Rscript -e "rmarkdown::render('laura_GSCs_GBM_combined_VM1_VM2.Rmd', output_dir = '${master_html_output_dir}')" 
# master_output_dir='/home/owenwhitley/projects/su2c_v2/results/scRNA/GSC_Tumour_combined_w_varimax_PC1_PC2'
# if [[ ! -d $master_output_dir ]]; then
#   mkdir $master_output_dir
# fi
# # Run GSEA on first 5 dims
# Rscript -e "rmarkdown::render('./seurat_dr_GSEA_generic.Rmd',\
#                               output_file = 'GSC_Tumour_combined_w_varimax_GSEA.html',\
#                               output_dir = '${master_html_output_dir}',\
#                               params = list(inp_data = file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L_w_varimax_PC1_PC2.rds'),\
#                                             output_dir = file.path('${master_output_dir}', 'GSEA_results'),\
#                                             dims_use = 1:2,\
#                                             dr_use = 'varimax',\
#                                             n_cores = 1)\
#                               )"

# # Run logistic regression notebook using 2 PCs
# results_dir_PC1_PC2='../../results/scRNA/laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic'
# if [ ! -d $results_dir_PC1_PC2 ];
# then
#   mkdir $results_dir_PC1_PC2
# fi
# 
# Rscript -e "BTSC_Tumour_combined_obj <- readRDS(file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L.rds'));\
# rmarkdown::render('./laura_SU2C_GBM_GSCs_combined_logistic.Rmd',\
#                   output_file = 'laura_SU2C_GBM_GSCs_combined_logistic.html',\
#                   output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic'),\
#                   params = list(rerun_logistic = TRUE,\
#                                 results_dir = '${results_dir_PC1_PC2}',\
#                                 file_prefix = 'logistic_regression_T_v_L',\
#                                 inp_data = BTSC_Tumour_combined_obj,\
#                                 n_splits = 30L,\
#                                 coef_thresh = 0.20,\
#                                 heatmap_fontsize = 10.0,\
#                                 base_seed = 32L,\
#                                 reduction_use = 'pca',\
#                                 other_reductions_plot = c('pca', 'tsne'),\
#                                 feats_plot = c('SampleType', 'pred_class', 'misclassified', 'SampleID'),\
#                                 num_dims = 2),\
#                   envir = new.env())\
# "
# 
# 
# ## correlation plot
# Rscript -e "BTSC_Tumour_combined_obj <- readRDS(file.path('${results_dir_PC1_PC2}', 'logistic_regression_T_v_L_seurat.rds'));\
# rmarkdown::render('./corplot_2_var.Rmd',\
#                   output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L'),\
#                   output_file = 'corplot_2_var.html',\
#                   params = list(inp_data = BTSC_Tumour_combined_obj,\
#                                 feats_plot = c('pred_class', 'SampleType'),\
#                                 cor_pairs = list(c('RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC')),\
#                                 filter_by = 'SampleType'),\
#                   envir = new.env());\
# gc(full = TRUE)\
# "
# 
# ## correlation plot, tumour cells only
# Rscript -e "library(Seurat);\
# seurat_obj <- readRDS(file.path('${results_dir_PC1_PC2}', 'logistic_regression_T_v_L_seurat.rds'));\
# tumour_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$SampleType == 'Tumour'];\
# seurat_obj <- Seurat::SubsetData(seurat_obj, cells.use = tumour_cells, subset.raw = TRUE);\
# rmarkdown::render('./corplot_2_var.Rmd',\
#                   output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L'),\
#                   output_file = 'corplot_2_var_tumour_only.html',\
#                   params = list(inp_data = seurat_obj,\
#                                 feats_plot = c('pred_class', 'SampleType'),\
#                                 cor_pairs = list(c('RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC'), c('glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC')),\
#                                 filter_by = 'pred_class'),\
#                   envir = new.env());\
# gc(full = TRUE)\
# "
# ## show tumour + line cells' feature scores, subcategorized by whether or not cells were properly classified
# Rscript -e "BTSC_Tumour_combined_obj <- readRDS(file.path('${results_dir_PC1_PC2}', 'logistic_regression_T_v_L_seurat.rds'));\
# boxplot_feats <- c('RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC');\
# rmarkdown::render('./scRNA_make_boxplots.Rmd',\
#                   params = list(inp_data = BTSC_Tumour_combined_obj,\
#                                 split_by = 'SampleType',\
#                                 color_by = 'pred_class',\
#                                 feats_plot = boxplot_feats),\
#                   output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L'),\
#                   envir = new.env());\
# gc(full = TRUE)"
#