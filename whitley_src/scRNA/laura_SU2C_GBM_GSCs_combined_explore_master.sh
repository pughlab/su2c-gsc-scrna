#!/bin/bash

# Overview: run scRNA exploratory analyses for GBM + GSCs

# setup preprocessed data dir + html output directory
preproc_dir="~/projects/su2c_v2/data/preprocessed/scRNA/GSCs_Tumour_combined"
master_html_output_dir="./reports/laura_SU2C_GSCs_Tumour_combined"

set -e

# # explore PCs with and without G800_L
Rscript -e "\
rmarkdown::render('./laura_SU2C_GBM_GSCs_combined_annotate_PCs.Rmd',\
                  output_file = 'laura_GBM_GSCs_annotate_PCs_no_G800_L.html',\
                  output_dir = file.path('${master_html_output_dir}', 'laura_GBM_GSCs_annotate_PCs_no_G800_L'),\
                  params = list(inp_file = file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L.rds'),\
                                output_dir = '../../results/scRNA/laura_GBM_GSCs_annotate_PCs_no_G800_L'),\
                  envir = new.env())\
"


Rscript -e "\
rmarkdown::render('./laura_SU2C_GBM_GSCs_combined_annotate_PCs.Rmd',\
                  output_file = 'laura_GBM_GSCs_annotate_PCs.html',\
                  output_dir = file.path('${master_html_output_dir}', 'laura_GBM_GSCs_annotate_PCs'),\
                  params = list(inp_file = file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019.rds'),\
                                output_dir = '../../results/scRNA/laura_GBM_GSCs_annotate_PCs'),\
                  envir = new.env())"
# Run logistic regression notebook using 2 PCs
results_dir_PC1_PC2='../../results/scRNA/laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic'
if [ ! -d $results_dir_PC1_PC2 ];
then
  mkdir $results_dir_PC1_PC2
fi

Rscript -e "BTSC_Tumour_combined_obj <- readRDS(file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L.rds'));\
rmarkdown::render('./laura_SU2C_GBM_GSCs_combined_logistic.Rmd',\
                  output_file = 'laura_SU2C_GBM_GSCs_combined_logistic.html',\
                  output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic'),\
                  params = list(rerun_logistic = TRUE,\
                                results_dir = '${results_dir_PC1_PC2}',\
                                file_prefix = 'logistic_regression_T_v_L',\
                                inp_data = BTSC_Tumour_combined_obj,\
                                n_splits = 30L,\
                                coef_thresh = 0.20,\
                                heatmap_fontsize = 10.0,\
                                base_seed = 32L,\
                                reduction_use = 'pca',\
                                other_reductions_plot = c('pca', 'tsne'),\
                                feats_plot = c('SampleType', 'pred_class', 'misclassified', 'SampleID'),\
                                num_dims = 2),\
                  envir = new.env())\
"

# Do the same but using 10 PCs
results_dir_10_PCs='../../results/scRNA/laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic_10_PCs'
if [ ! -d $results_dir_10_PCs ];
then
  mkdir $results_dir_10_PCs
fi

Rscript -e " BTSC_Tumour_combined_obj <- readRDS(file.path('${preproc_dir}', 'laura_BTSC_Tumour_combined_Oct_2019_no_G800_L.rds'));\
rmarkdown::render('./laura_SU2C_GBM_GSCs_combined_logistic.Rmd',\
                  output_file = 'laura_SU2C_GBM_GSCs_combined_logistic_10_PCs.html',\
                  output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L_logistic_10_PCs'),\
                  params = list(rerun_logistic = TRUE,\
                                results_dir = '${results_dir_10_PCs}',\
                                file_prefix = 'logistic_regression_T_v_L',\
                                inp_data = BTSC_Tumour_combined_obj,\
                                n_splits = 30L,\
                                coef_thresh = 0.20,\
                                heatmap_fontsize = 10.0,\
                                base_seed = 32L,\
                                reduction_use = 'pca',\
                                other_reductions_plot = c('pca', 'tsne'),\
                                feats_plot = c('SampleType', 'pred_class', 'misclassified', 'SampleID'),\
                                num_dims = 10),\
                  envir = new.env())\
"
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
# ## correlation plot, tumour cells only
# ## TODO: This Rscript call fails for some reason with a failure in cor.test (not enough samples)
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
## show tumour + line cells' feature scores, subcategorized by whether or not cells were properly classified
Rscript -e "BTSC_Tumour_combined_obj <- readRDS(file.path('${results_dir_PC1_PC2}', 'logistic_regression_T_v_L_seurat.rds'));\
boxplot_feats <- c('RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC');\
rmarkdown::render('./scRNA_make_boxplots.Rmd',\
                  params = list(inp_data = BTSC_Tumour_combined_obj,\
                                split_by = 'SampleType',\
                                color_by = 'pred_class',\
                                feats_plot = boxplot_feats),\
                  output_dir = file.path('${master_html_output_dir}', 'laura_SU2C_GSCs_GBM_combined_Oct_2019_no_G800_L'),\
                  envir = new.env());\
gc(full = TRUE)\
"

# # Find genes correlated with Developmental (RNA GSC cluster 1) and Injury Response (RNA GSC cluster 2) signature scores (AUC) in GSC lines + tumour cells
# # scripts here relies on output of previous analyses run from this script (uses tumour cells with logistic regression classification label). don't need that classification now but may need later.
# Rscript -e "rmarkdown::render('laura_SU2C_GSCs_AUCell_gene_correlation.Rmd', output_dir = 'reports/AUCell_genes')" > laura_SU2C_GSCs_gene_correlation_out.txt 2 > laura_SU2C_GSCs_AUCell_gene_correlation_stderr.txt
# Rscript -e "rmarkdown::render('laura_SU2C_GBM_AUCell_gene_correlation.Rmd', output_dir = 'reports/AUCell_genes')" > laura_SU2C_GBM_AUCell_gene_correlation_out.txt 2 > laura_SU2C_GBM_AUCell_gene_correlation_stderr.txt
# 
# explore tumour + GSC lines on their own. GBM exploration relies on class labels from logistic regression script
Rscript -e "rmarkdown::render('laura_SU2C_GBM_Nov_2019_explore.Rmd', output_dir = 'reports/laura_SU2C_GSCs_GBM_Nov_2019_explore', envir = new.env());\
gc(full = TRUE)"
Rscript -e "rmarkdown::render('laura_SU2C_GSCs_Nov_2019_explore.Rmd', output_dir = 'reports/laura_SU2C_GSCs_GBM_Nov_2019_explore', envir = new.env());\
gc(full = TRUE)"
