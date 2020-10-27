###############################################################################
## Do Exploration of Wang 2019 data after filtering for glioma cells
library(knitr)
library(rmarkdown)
library(Seurat)
suppressPackageStartupMessages({source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')})
# setup
html_output_dir <- 'reports/Wang_2019_filtered_analyses'
preproc_data_dir <- '~/projects/su2c_v2/data/preprocessed/scRNA/Wang_2019_CancerDiscov_Tumour_Filtered/'
dir(preproc_data_dir)
scRNA_file <- "Wang_2019_scRNA_filtered.rds"
Wang_2019_scRNA <- readRDS(file = file.path(preproc_data_dir, scRNA_file))
snRNA_file <- 'Wang_2019_snRNA_filtered.rds'
Wang_2019_snRNA <- readRDS(file = file.path(preproc_data_dir, snRNA_file))
# calculate score difference features
Wang_2019_scRNA <- calc_feat_diff(Wang_2019_scRNA, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
Wang_2019_scRNA <- calc_feat_diff(Wang_2019_scRNA, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
Wang_2019_snRNA <- calc_feat_diff(Wang_2019_snRNA, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
Wang_2019_snRNA <- calc_feat_diff(Wang_2019_snRNA, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
# define features to plot
feats_plot <- colnames(Wang_2019_scRNA@meta.data)[grep('(glioma.stem.cell|Developmental|Injury_Response|cahoy_(astrocyte|oligodendrocyte|neuron))_AUC$|Diff|Study_ID|nUMI',
                                                           colnames(Wang_2019_scRNA@meta.data))]
feats_plot <- c(feats_plot, 'SampleID')
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
                                  feats_plot = features),
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
                                  filter_by = 'SampleID'),
                    envir = new.env()
  )
  gc(full = TRUE)
  print('Finished making notebooks')
  print(Sys.time())
}

analysis_routine(Wang_2019_scRNA, html_output_dir, feats_plot, 'Wang_2019_scRNA_filtered')
analysis_routine(Wang_2019_snRNA, html_output_dir, feats_plot, 'Wang_2019_snRNA_filtered')
# GSM4119525 is an outlier. rerun without this sample
cells_use <- rownames(Wang_2019_snRNA@meta.data)[Wang_2019_snRNA@meta.data$SampleID != 'GSM4119525']
data_use <- Wang_2019_snRNA@raw.data[,cells_use]
meta_data <- Wang_2019_snRNA@meta.data[cells_use,]
rm(Wang_2019_snRNA)
gc(full = TRUE)
Wang_2019_snRNA_no_GSM4119525 <- seurat_subroutine(data_use, meta_data)
# Load Genesets
top_dir <- '~/projects/su2c_v2'
gene_sets_dir <- file.path(top_dir, '/data/preprocessed/GeneSets')
gene_sets_file <- 'genesets_and_info.rds'
genesets_and_info <- readRDS(file = file.path(gene_sets_dir, gene_sets_file))
genesets <- genesets_and_info$gene_set_list
# rename RNA.GSC.c1, RNA.GSC.c2 to Developmental, Injury Response, repectively
genesets$Developmental <- genesets$RNA.GSC.c1
genesets$RNA.GSC.c1 <- NULL
genesets$Injury_Response <- genesets$RNA.GSC.c2
genesets$RNA.GSC.c2 <- NULL
# calculate scores + score diffs
output_dir <- '~/projects/su2c_v2/data/preprocessed/scRNA/Wang_2019_CancerDiscov_Tumour_Filtered_no_GSM411925/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
Wang_2019_snRNA_no_GSM4119525 <- scoring_subroutine(Wang_2019_snRNA_no_GSM4119525, genesets, preproc_data_dir = output_dir, prefix = 'Wang_2019_snRNA_no_GSM4119525')
Wang_2019_snRNA_no_GSM4119525 <- calc_feat_diff(Wang_2019_snRNA_no_GSM4119525, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
Wang_2019_snRNA_no_GSM4119525 <- calc_feat_diff(Wang_2019_snRNA_no_GSM4119525, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
saveRDS(Wang_2019_snRNA_no_GSM4119525, file = file.path(output_dir, 'Wang_2019_snRNA_no_GSM4119525.rds'))
# redo plots
analysis_routine(Wang_2019_snRNA_no_GSM4119525, html_output_dir, feats_plot, 'Wang_2019_snRNA_no_GSM4119525')

sessionInfo()
