###############################################################################
## Address Reviewer Comment Regarding Regressing Out Cell Cycle Effects
## Produce PCA plots resulting when S.Score and G2M score regressed out directly
source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')
data_dir <- '~/projects/su2c_v2/data/raw/scRNA/Data//'
tumour_gsc_dir <- file.path(data_dir, 'SU2C_BTSC_TumourCells_combined_Oct_2019')
# gsc_dir <- file.path(data_dir, 'GSCs_G800Removed_Nov2019')

genesets_dir <- '~/projects/su2c_v2/data/preprocessed/GeneSets/'
genesets_and_info <- readRDS(file.path(genesets_dir, 'genesets_and_info.rds'))
genesets <- genesets_and_info$gene_set_list
reports_dir <- '~/projects/su2c_v2/scripts/scRNA/reports/laura_rebuttal_July_2020_PC_gene_selection'

genesets_regex <- 'Neftel|cahoy_(astrocyte|oligodendrocyte)|RNA.GSC.c[1-2]|glioma.stem.cell'
genesets_score <- genesets[grep(genesets_regex, names(genesets))]


# function to calculate number of samples per clusteer
samples_per_cluster <- function(seurat_obj) {
  meta_data <- seurat_obj@meta.data
  counts <- c()
  for (clust_id in unique(meta_data$cluster)) {
    counts <- c(counts, length(unique(meta_data$SampleID[meta_data$cluster == clust_id])))
  }
  return(counts)
}

# run for combined tumour + line cells
load(file.path(tumour_gsc_dir, 'BTSC_TumourCells_G800Removed_AUCell_Seurat.Rdata'))
BTSC_TumourCells <- subset_seurat(BTSC_TumourCells, subset_size = 20000)
sample_id_regex <- '(G|BT)[0-9]+-*[A-Z]*_[TL]'
cell_ids <- rownames(BTSC_TumourCells@meta.data)
BTSC_TumourCells@meta.data$SampleID <- regmatches(cell_ids, regexpr(sample_id_regex, cell_ids))
BTSC_TumourCells@meta.data$SampleID <- sub('-[A-Z]', '', BTSC_TumourCells@meta.data$SampleID)
BTSC_TumourCells <- AUCell_batch(BTSC_TumourCells, genesets = genesets_score, num_batches = 1)
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'RNA.GSC.c1_AUC', 'RNA.GSC.c2_AUC', 'Dev_IR_Diff')
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
BTSC_TumourCells@meta.data$Dev_IR_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Dev_IR_Diff, scale = TRUE, center = TRUE)
BTSC_TumourCells@meta.data$Stem_Astro_Diff_zscore <- scale(BTSC_TumourCells@meta.data$Stem_Astro_Diff, scale = TRUE, center = TRUE)

feat_regex <- 'Neftel|orig.ident|cahoy_(astrocyte|oligodendrocyte)|zscore|CC.Difference|S.Score|G2M.Score'
feats_plot_common <- colnames(BTSC_TumourCells@meta.data)[grep(feat_regex, colnames(BTSC_TumourCells@meta.data))]
feats_plot_GSC_Tumour <- c('SampleID', feats_plot_common)


BTSC_TumourCells <- seurat_pipeline(BTSC_TumourCells, vars_regress = c('percent.mito', 'nUMI', 'CC.Difference'))
BTSC_TumourCells <- RunTSNE(BTSC_TumourCells, reduction.use = 'pca', dims.use = 1:5)
BTSC_TumourCells@meta.data$res.0.8 <- NULL
BTSC_TumourCells <- FindClusters(BTSC_TumourCells, dims.use = 1:5, force.recalc = TRUE)
BTSC_TumourCells@meta.data$cluster <- BTSC_TumourCells@meta.data$res.0.8
gc(full = TRUE)
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_no_ribo.html',
                  params = list(inp_data = BTSC_TumourCells,
                                feats_plot = feats_plot_GSC_Tumour,
                                dr_use = c('pca', 'tsne')
                  ),
                  envir = new.env()
)
rmarkdown::render('seurat_barcharts_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_no_ribo_clusters_barcharts.html',
                  params = list(inp_data = BTSC_TumourCells,
                                x_axis_var = 'cluster', 
                                y_axis_var = 'SampleID'),
                  envir = new.env()
)

samples_per_cluster_no_ribo <- samples_per_cluster(BTSC_TumourCells)

# Do PCA on top 500 and top 2000 most variable genes in scaled data, excluding ribosomal genes
ribo.genes <- c(rownames(BTSC_TumourCells@scale.data)[grep("^RP[1:9]", rownames(BTSC_TumourCells@data))],
                rownames(BTSC_TumourCells@scale.data)[grep("^RP[L,S]", rownames(BTSC_TumourCells@data))])
valid.genes <- setdiff(rownames(BTSC_TumourCells@scale.data), ribo.genes)
top.n.genes <- names(sort(rowVars(BTSC_TumourCells@scale.data[valid.genes,]), decreasing = TRUE))
# do top 500 genes
BTSC_TumourCells <- RunPCA(BTSC_TumourCells, pc.genes = top.n.genes[1:500])
BTSC_TumourCells <- RunTSNE(BTSC_TumourCells, reduction.use = 'pca', dims.use = 1:5)
BTSC_TumourCells@meta.data$res.0.8 <- NULL
BTSC_TumourCells <- FindClusters(BTSC_TumourCells, dims.use = 1:5, force.recalc = TRUE)
BTSC_TumourCells@meta.data$cluster <- BTSC_TumourCells@meta.data$res.0.8
rmarkdown::render('seurat_dr_plot_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_no_ribo_top_500.html',
                  params = list(inp_data = BTSC_TumourCells,
                                feats_plot = feats_plot_GSC_Tumour,
                                dr_use = c('pca', 'tsne')
                  ),
                  envir = new.env()
)
rmarkdown::render('seurat_barcharts_generic.Rmd',
                  output_dir = reports_dir,
                  output_file = 'GSCs_Tumours_PCA_no_ribo_top_500_clusters_barcharts.html',
                  params = list(inp_data = BTSC_TumourCells,
                                x_axis_var = 'cluster', 
                                y_axis_var = 'SampleID'),
                  envir = new.env()
)
gc(full = TRUE)

samples_per_cluster_no_ribo_top500 <- samples_per_cluster(BTSC_TumourCells)

# compare number of samples per cluster obtained using different sets of genes for PCA
p_val <- wilcox.test(samples_per_cluster_no_ribo, samples_per_cluster_no_ribo_top500,
                     alternative = 'two.sided')$p.value
png(file = file.path(reports_dir, 'all_genes_vs_top_500_samples_per_cluster.png'))
boxplot_list <- list(all_genes = samples_per_cluster_no_ribo,
                     top_500 = samples_per_cluster_no_ribo_top500)
boxplot(boxplot_list,
        ylab = 'samples per cluster',
        main = paste('wilcoxon RS test p-value:', p_val))
dev.off()
# # do top 2000 genes
# BTSC_TumourCells <- RunPCA(BTSC_TumourCells, pc.genes = top.n.genes[1:2000])
# BTSC_TumourCells <- RunTSNE(BTSC_TumourCells, reduction.use = 'pca', dims.use = 1:5)
# BTSC_TumourCells@meta.data$res.0.8 <- NULL
# BTSC_TumourCells <- FindClusters(BTSC_TumourCells, dims.use = 1:5, force.recalc = TRUE)
# BTSC_TumourCells@meta.data$cluster <- BTSC_TumourCells@meta.data$res.0.8
# rmarkdown::render('seurat_dr_plot_generic.Rmd',
#                   output_dir = reports_dir,
#                   output_file = 'GSCs_Tumours_PCA_no_ribo_top_2000.html',
#                   params = list(inp_data = BTSC_TumourCells,
#                                 feats_plot = feats_plot_GSC_Tumour,
#                                 dr_use = c('pca', 'tsne')
#                   ),
#                   envir = new.env()
# )
# rmarkdown::render('seurat_barcharts_generic.Rmd',
#                   output_dir = reports_dir,
#                   output_file = 'GSCs_Tumours_PCA_no_ribo_top_2000_clusters_barcharts.html',
#                   params = list(inp_data = BTSC_TumourCells,
#                                 x_axis_var = 'cluster', 
#                                 y_axis_var = 'SampleID'),
#                   envir = new.env()
# )
# gc(full = TRUE)


sessionInfo()
