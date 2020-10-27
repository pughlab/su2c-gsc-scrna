###############################################################################
## Diffusion Map on Subset of GSCs

source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')

top_dir <- '~/projects/su2c_v2'
data_dir <- file.path(top_dir, 'data')
raw_data_dir <- file.path(data_dir, 'raw/scRNA/Data/GSCs_G800Removed_Nov2019')
raw_data_file <- 'BTSC_G800L_removed_AUCell_SeuratObj.rds'
gene_sets_dir <- file.path(data_dir, 'preprocessed/GeneSets')
output_dir <- file.path(data_dir, 'preprocessed/scRNA/GSCs_Diffusion_Map_Sep_2020')
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

gene_sets_file <- 'genesets_and_info.rds'
genesets_and_info <- readRDS(file = file.path(gene_sets_dir, gene_sets_file))
genesets <- genesets_and_info$gene_set_list
# rename RNA.GSC.c1, RNA.GSC.c2 to Developmental, Injury Response, repectively
genesets$Developmental <- genesets$RNA.GSC.c1
genesets$RNA.GSC.c1 <- NULL
genesets$Injury_Response <- genesets$RNA.GSC.c2
genesets$RNA.GSC.c2 <- NULL
# filter for specific genesets
geneset_regex <- 'Developmental|Injury_Response|Neftel|cahoy_(astrocyte|oligodendrocyte|neuron)|glioma.stem.cell|microglia'
genesets <- genesets[grep(geneset_regex, names(genesets), ignore.case = TRUE)]

seurat_obj <- readRDS(file = file.path(raw_data_dir, raw_data_file))

# remove AUC cols
auc_cols <- which(grep('AUC$', colnames(seurat_obj@meta.data)))
seurat_obj@meta.data <- seurat_obj@meta.data[,-auc_cols]