###############################################################################
### Preprocess Neftel et al. 2019 data
source('laura_SU2C_GBM_GSCs_combined_explore_helpers.R')
library(readxl)
library(biomaRt)
top_dir <- '~/projects/su2c_v2'
data_dir <- file.path(top_dir, 'data')
gene_sets_dir <- file.path(data_dir, 'preprocessed/GeneSets')
inp_dir <- file.path(data_dir, 'raw/scRNA/Neftel_2019')
preproc_data_dir <- file.path(data_dir, 'preprocessed/scRNA/Neftel_2019')
if (!dir.exists(preproc_data_dir)) {
  dir.create(preproc_data_dir)
}
###############################################################################
## Get mappings for chromosomes
listEnsemblArchives()
ensembl_hg19 <- useMart(host="http://grch37.ensembl.org", 
                           biomart='ENSEMBL_MART_ENSEMBL')

# listDatasets(ensembl_hg19)[grep('hsapiens', listDatasets(ensembl_hg19)$dataset), ]
ensembl_hg19 <- useDataset(dataset = 'hsapiens_gene_ensembl', 
                              mart = ensembl_hg19)
attr <- listAttributes(ensembl_hg19)
download_biomart <- TRUE
if (!download_biomart) {
  reload <- TRUE
  tryCatch({
    BM.mapping <- readRDS(file = file.path(preproc_data_dir, 'BM_mapping.rds'))
    reload <- FALSE
  })
  if (reload) {
    print('reloading ENSEMBL biomart mappings as none detected in preproc_data_dir')
  }
  download_biomart <- reload
}
if (download_biomart) {
  BM.mapping <- getBM(attributes =  c('hgnc_symbol',
                                      'chromosome_name'), 
                      mart = ensembl_hg19)
  saveRDS(BM.mapping, file = file.path(preproc_data_dir, 'BM_mapping.rds'))
}

###############################################################################
## Load Data ##
print('Loading Data')
print(Sys.time())
# Meta Data. Filter fo only 10X data
meta_data_file <- "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
meta_table <- as.data.frame(readxl::read_xlsx(file.path(inp_dir, meta_data_file), skip = 43))
meta_table <- meta_table[meta_table$`processed data file` == '10X_GBM_IDHwt_processed_TPM.tsv',]
rownames(meta_table) <- meta_table$`Sample name`
# Expression Data.
exp_data_file <- "GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv"
exp_data <- read.delim(file.path(inp_dir, exp_data_file), header = TRUE, row.names = 'GENE', stringsAsFactors = FALSE)
exp_data <- as.matrix(exp_data)
gc(full = TRUE)
# remove X out of colnames of exp_data, to match with meta_data
colnames(exp_data) <- sub('^X', '', colnames(exp_data))

# Gene Sets. Note that Gene sets are based on hgnc from hg38. hg19 was used for alignment here,
# but assumption is that vast majority of hgnc symbols, if kept between hg19 and hg38, refer to the same
# feature
gene_sets_file <- 'genesets_and_info.rds'
genesets_and_info <- readRDS(file = file.path(gene_sets_dir, gene_sets_file))
print(Sys.time())
## Run Checks ##
print('All cell ids match between metadata + table')
all(sort(colnames(exp_data)) == sort(rownames(meta_table)))

## Filter for Adult samples
print('cells, by age group')
print(table(meta_table$`adult/pediatric`))
# appears that no further action is needed w.r.t. subsetting data

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

print('Running Seurat Pipeline + Scoring')
print(table(meta_table$Tumor.Normal.Classification))
# filtering out low quality cells was already done by Neftel et al. Was not done for low expression genes
neftel_2019_seurat <- seurat_subroutine(raw_data = exp_data, meta_data = meta_table, min_genes = 20, min_cells = 0)
neftel_2019_seurat <- Seurat::RunTSNE(neftel_2019_seurat, reduction.use = 'pca', dims.use = 1:5)
print('Saving Results')
print(Sys.time())
output_fname <- 'neftel_2019_seurat.rds'
saveRDS(neftel_2019_seurat, file = file.path(preproc_data_dir, output_fname))
## Run Chromosome + Signature Scoring
# calculate z-scored avg chromosome expression
chr_avg_output_list <- calcAvgChrMat(neftel_2019_seurat, BM.mapping, chr.use = as.character(1:22))
chr_mat <- chr_avg_output_list$output.mat
rownames(chr_mat) <- paste0('chr.', rownames(chr_mat))
neftel_2019_seurat <- Seurat::AddMetaData(neftel_2019_seurat, as.data.frame(t(chr_mat)))
# signature scoring
neftel_2019_seurat <- scoring_subroutine(neftel_2019_seurat, genesets, preproc_data_dir = preproc_data_dir, prefix = 'neftel_2019')
# add feature difference scores
neftel_2019_seurat <- calc_feat_diff(neftel_2019_seurat, 'Developmental_AUC', 'Injury_Response_AUC', 'Dev_IR_Diff')
neftel_2019_seurat <- calc_feat_diff(neftel_2019_seurat, 'glioma.stem.cell_AUC', 'cahoy_astrocyte_AUC', 'Stem_Astro_Diff')
print('Saving Results')
print(Sys.time())
output_fname2 <- 'neftel_2019_seurat_scored.rds'
saveRDS(neftel_2019_seurat, file = file.path(preproc_data_dir, output_fname2))
print('Finished')
print(Sys.time())

sessionInfo()