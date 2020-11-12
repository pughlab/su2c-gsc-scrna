###############################################################################
## Address Reviewer Comment Regarding Regressing Out Cell Cycle Effects
## Produce PCA plots resulting when S.Score and G2M score regressed out directly
source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')
data_dir <- '~/projects/su2c_v2/data/raw/scRNA/Data//'
tumour_gsc_dir <- file.path(data_dir, 'SU2C_BTSC_TumourCells_combined_Oct_2019')
gsc_dir <- file.path(data_dir, 'GSCs_G800Removed_Nov2019')

genesets_dir <- '~/projects/su2c_v2/data/preprocessed/GeneSets/'
genesets_and_info <- readRDS(file.path(genesets_dir, 'genesets_and_info.rds'))
genesets <- genesets_and_info$gene_set_list
reports_dir <- '~/projects/su2c_v2/scripts/scRNA/reports/laura_rebuttal_July_2020_Cell_Cycle_Regress_Out'

genesets_regex <- 'Neftel|cahoy_(astrocyte|oligodendrocyte)|RNA.GSC.c[1-2]|glioma.stem.cell'
genesets_score <- genesets[grep(genesets_regex, names(genesets))]

# run for combined tumour + line cells
load(file.path(tumour_gsc_dir, 'BTSC_TumourCells_G800Removed_AUCell_Seurat.Rdata'))
BTSC_TumourCells <- subset_seurat(BTSC_TumourCells, subset_size = 20000)
BTSC_TumourCells <- AUCell_batch(BTSC_TumourCells, genesets = genesets_score, num_batches = 1)
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff')
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
BTSC_TumourCells@meta.data$Dev_IR_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
BTSC_TumourCells@meta.data$Stem_Astro_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)
sample_id_regex <- '(G|BT)[0-9]+-*[A-Z]*_[TL]'
cell_ids <- rownames(BTSC_TumourCells@meta.data)
BTSC_TumourCells@meta.data$SampleID <- regmatches(cell_ids, regexpr(sample_id_regex, cell_ids))
BTSC_TumourCells@meta.data$SampleID <- sub('-[A-Z]', '', BTSC_TumourCells@meta.data$SampleID)

feat_regex <- 'Neftel|orig.ident|cahoy_(astrocyte|oligodendrocyte)|zscore|CC.Difference|S.Score|G2M.Score'
feats_plot_common <- colnames(BTSC_TumourCells@meta.data)[grep(feat_regex, colnames(BTSC_TumourCells@meta.data))]
feats_plot_GSC_Tumour <- c('SampleID', feats_plot_common)


BTSC_TumourCells <- seurat_pipeline(BTSC_TumourCells, vars_regress = c('percent.mito', 'nUMI', 'S.Score', 'G2M.Score'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_G2M_S_regress.html',
                  params = list(inp_data = BTSC_TumourCells,
                                feats_plot = feats_plot_GSC_Tumour,
                                dr_use = c('pca', 'tsne')
                                ),
                  envir = new.env())
gc(full = TRUE)
BTSC_TumourCells <- seurat_pipeline(BTSC_TumourCells, vars_regress = c('percent.mito', 'nUMI', 'CC.Difference'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_CC_diff_regress.html',
                  params = list(inp_data = BTSC_TumourCells,
                                feats_plot = feats_plot_GSC_Tumour,
                                dr_use = c('pca', 'tsne')
                                ),
                  envir = new.env()
)
BTSC_TumourCells <- seurat_pipeline(BTSC_TumourCells, vars_regress = c('percent.mito', 'nUMI'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_no_CC_regress.html',
                  params = list(inp_data = BTSC_TumourCells,
                                feats_plot = feats_plot_GSC_Tumour,
                                dr_use = c('pca', 'tsne')
                  ),
                  envir = new.env()
)

rm(BTSC_TumourCells)

gc(full = TRUE)

# run for line cells alone
GSC_seurat <- readRDS(file.path(gsc_dir, "BTSC_G800L_removed_AUCell_SeuratObj.rds"))
GSC_seurat <- subset_seurat(GSC_seurat, subset_size = 20000)
GSC_seurat <- AUCell_batch(GSC_seurat, genesets = genesets_score, num_batches = 1)
GSC_seurat <- calc_feat_diff(GSC_seurat, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff')
GSC_seurat <- calc_feat_diff(GSC_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
GSC_seurat@meta.data$Dev_IR_Diff_zscore <- scale(GSC_seurat@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
GSC_seurat@meta.data$Stem_Astro_Diff_zscore <- scale(GSC_seurat@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)

feats_plot_GSC <- feats_plot_common

GSC_seurat <- seurat_pipeline(GSC_seurat, vars_regress = c('percent.mito', 'nUMI', 'S.Score', 'G2M.Score'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_PCA_G2M_S_regress.html',
                  params = list(inp_data = GSC_seurat,
                                feats_plot = feats_plot_GSC,
                                dr_use = c('pca', 'tsne')
                                ),
                  envir = new.env()
)
gc(full = TRUE)
GSC_seurat <- seurat_pipeline(GSC_seurat, vars_regress = c('percent.mito', 'nUMI', 'CC.Difference'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_PCA_CC_diff_regress.html',
                  params = list(inp_data = GSC_seurat,
                                feats_plot = feats_plot_GSC,
                                dr_use = c('pca', 'tsne')
                                ),
                  envir = new.env()
                  )
gc(full = TRUE)
GSC_seurat <- seurat_pipeline(GSC_seurat, vars_regress = c('percent.mito', 'nUMI'))
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_PCA_no_CC_regress.html',
                  params = list(inp_data = GSC_seurat,
                                feats_plot = feats_plot_GSC,
                                dr_use = c('pca', 'tsne')
                  ),
                  envir = new.env()
)
gc(full = TRUE)
sessionInfo()
