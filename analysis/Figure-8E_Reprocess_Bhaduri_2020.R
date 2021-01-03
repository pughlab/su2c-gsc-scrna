###############################################################################
### Preprocess Bhaduri et al. 2020 data
source('whitley_scRNA_helpers.R')
top_dir <- '~/projects/su2c_v2'
data_dir <- file.path(top_dir, 'data')
gene_sets_dir <- file.path(data_dir, 'preprocessed/GeneSets')
inp_dir <- file.path(data_dir, 'raw/scRNA/Bhaduri_2020')
preproc_data_dir <- file.path(data_dir, 'preprocessed/scRNA/Bhaduri_2020')
if (!dir.exists(preproc_data_dir)) {
  dir.create(preproc_data_dir)
}

###############################################################################
## Subroutines
seurat_subroutine <- function(raw_data, meta_data) {
  print('Creating Seurat Object')
  seurat_obj <- Seurat::CreateSeuratObject(raw.data = raw_data, project = "SeuratProject", min.cells = 20,
                                           min.genes = 200, is.expr = 0, do.center = FALSE, meta.data = meta_data)
  print('initial dims after filtering for min 20 cells per gene, min 200 genes per cell')
  dim(seurat_obj@raw.data)
  print('running Seurat Pipeline')
  print(Sys.time())
  seurat_obj <- seurat_pipeline(seurat_obj)
  seurat_obj <- Seurat::RunTSNE(seurat_obj)
  seurat_obj <- Seurat::RunDiffusion(seurat_obj)
  print(Sys.time())
  print('Clustering')
  seurat_obj <- Seurat::FindClusters(seurat_obj, print.output = FALSE)
  return(seurat_obj)
}


scoring_subroutine <- function(seurat_obj, genesets, preproc_data_dir, prefix) {
  print('Filtering Genesets')
  orig_length <- unlist(lapply(genesets, length))
  orig_gensets_condensed <- unlist(lapply(genesets, FUN = function(x) {paste(x, collapse = ';')}))
  genesets_filt <- lapply(genesets, FUN = function(x) {return(x[x %in% rownames(seurat_obj@raw.data)])})
  filt_length <- unlist(lapply(genesets_filt, length))
  filt_gensets_condensed <- unlist(lapply(genesets_filt, FUN = function(x) {paste(x, collapse = ';')}))
  kept <- filt_length >= 5
  genesets_summary <- data.frame(geneset_name = names(genesets), 
                                 original_length = orig_length, 
                                 original_genes = orig_gensets_condensed,
                                 filtered_length = filt_length, 
                                 filtered_genes = filt_gensets_condensed,
                                 kept = kept)
  genesets_summary_file <- paste0(prefix, '_genesets_summary.csv')
  write.csv(genesets_summary, file = file.path(preproc_data_dir, genesets_summary_file))
  print('Running Scoring')
  print(Sys.time())
  # Run Scoring
  seurat_obj <- AUCell_batch(inp_data = seurat_obj, genesets = genesets_filt, num_batches = 1)
  print(Sys.time())
  return(seurat_obj)
}
###############################################################################
## Load Data ##
print('Loading Data')
print(Sys.time())
# Meta Data
meta_data_file <- 'meta.tsv'
meta_table <- read.delim(file.path(inp_dir, meta_data_file), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# Expression Data
exp_data_file <- 'exprMatrix.tsv'
exp_data <- read.delim(file.path(inp_dir, exp_data_file), header = TRUE, row.names = 'gene', stringsAsFactors = FALSE)
exp_data <- as.matrix(exp_data)
# Gene Sets
gene_sets_file <- 'genesets_and_info.rds'
genesets_and_info <- readRDS(file = file.path(gene_sets_dir, gene_sets_file))
print(Sys.time())
## Run Checks ##
print('All cell ids match between metadata + table')
all(sort(colnames(exp_data)) == sort(rownames(meta_table)))
print('number of common ids')
common_ids <- intersect(colnames(exp_data), rownames(meta_table))
length(common_ids)
print('number of non-shared ids setdiff(union, intersect)')
diff_ids <- setdiff(union(colnames(exp_data), rownames(meta_table)), common_ids)
length(diff_ids)
head(diff_ids[diff_ids %in% colnames(exp_data)])
head(diff_ids[diff_ids %in% rownames(meta_table)])
print('it appears that column names in exp_data matrix were renamed with dots instead of dashes, substitute dots + dashes with underscore')
colnames(exp_data) <- gsub('\\.', '_', colnames(exp_data))
rownames(meta_table) <- gsub('-', '_', rownames(meta_table))
print('All cell ids match between metadata + table')
all(sort(colnames(exp_data)) == sort(rownames(meta_table)))

print('dimensions of data')
print(dim(exp_data))
print('dimensions of metadata')
print(dim(meta_table))

genesets <- genesets_and_info$gene_set_list
# rename RNA.GSC.c1, RNA.GSC.c2 to Developmental, Injury Response, repectively
genesets$Developmental <- genesets$RNA.GSC.c1
genesets$RNA.GSC.c1 <- NULL
genesets$Injury_Response <- genesets$RNA.GSC.c2
genesets$RNA.GSC.c2 <- NULL
# filter for specific genesets
geneset_regex <- 'Developmental|Injury_Response|Neftel|cahoy_(astrocyte|oligodendrocyte|neuron)|glioma.stem.cell|microglia'
genesets <- genesets[grep(geneset_regex, names(genesets), ignore.case = TRUE)]

# ## Load in Biomart mapping (hgnc to chr) for chromosome average quantification
# listEnsemblArchives()
# listMarts(host = "http://Dec2015.archive.ensembl.org")
# 
# ensembl_100 <- useMart(host="http://apr2020.archive.ensembl.org", 
#                       biomart='ENSEMBL_MART_ENSEMBL')
# 
# listDatasets(ensembl_100)[grep('hsapiens', listDatasets(ensembl_100)$dataset), ]
# ensembl_100 <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = ensembl_100)
# attr <- listAttributes(ensembl_100)
# BM.mapping <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'hgnc_symbol'),  mart = ensembl_100)

## Run Preprocessing on Tumour Cells only
meta_table$Tumor.Normal.Classification[meta_table$Tumor.Normal.Classification == ''] <- 'unclassified'
tumor_cells <- rownames(meta_table)[meta_table$Tumor.Normal.Classification == 'Tumor']
print('Using Tumor Cells only')
print(table(meta_table$Tumor.Normal.Classification))
bhaduri_2020_seurat <- seurat_subroutine(raw_data = exp_data[,tumor_cells], meta_data = meta_table[tumor_cells,])
print('Saving Results')
print(Sys.time())
output_fname <- 'bhaduri_2020_tumor_seurat.rds'
saveRDS(bhaduri_2020_seurat, file = file.path(preproc_data_dir, output_fname))
bhaduri_2020_seurat <- scoring_subroutine(bhaduri_2020_seurat, genesets, preproc_data_dir = preproc_data_dir, prefix = 'bhaduri_2020_tumor')
# add feature difference scores
bhaduri_2020_seurat <- calc_feat_diff(bhaduri_2020_seurat, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
bhaduri_2020_seurat <- calc_feat_diff(bhaduri_2020_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
print('Saving Results')
print(Sys.time())
output_fname2 <- 'bhaduri_2020_tumor_seurat_scored.rds'
saveRDS(bhaduri_2020_seurat, file = file.path(preproc_data_dir, output_fname2))
print('Finished')
print(Sys.time())

sessionInfo()