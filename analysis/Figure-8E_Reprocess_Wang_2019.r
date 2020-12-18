###############################################################################
library(Seurat)
library(biomaRt)
source('whitley_scRNA_helpers.R')
top_dir <- '~/projects/su2c_v2'
raw_data_dir <- file.path(top_dir, 'data/raw/scRNA/Wang_2019_CancerDiscov')
preproc_data_dir <- file.path(top_dir, 'data/preprocessed/scRNA/Wang_2019_CancerDiscov')
gene_sets_dir <- file.path(top_dir, '/data/preprocessed/GeneSets')

if (!dir.exists(preproc_data_dir)) {
  dir.create(preproc_data_dir)
}

# Load Genesets
gene_sets_file <- 'genesets_and_info.rds'
genesets_and_info <- readRDS(file = file.path(gene_sets_dir, gene_sets_file))
genesets <- genesets_and_info$gene_set_list
# rename RNA.GSC.c1, RNA.GSC.c2 to Developmental, Injury Response, repectively
genesets$Developmental <- genesets$RNA.GSC.c1
genesets$RNA.GSC.c1 <- NULL
genesets$Injury_Response <- genesets$RNA.GSC.c2
genesets$RNA.GSC.c2 <- NULL

all_files <- dir(raw_data_dir, recursive = TRUE)
matrix_files <- all_files[grep('matrix.gene_vs_barcode.tsv$', all_files)]

meta_data <- data.frame(CellID = character(0), SampleID = character(0))

ensembl_current <- useMart(host="www.ensembl.org", 
                      biomart='ENSEMBL_MART_ENSEMBL')

# listDatasets(ensembl_current)[grep('hsapiens', listDatasets(ensembl_current)$dataset), ]
ensembl_current <- useDataset(dataset = 'hsapiens_gene_ensembl', 
                         mart = ensembl_current)
attr <- listAttributes(ensembl_current)
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
                      mart = ensembl_current)
  saveRDS(BM.mapping, file = file.path(preproc_data_dir, 'BM_mapping.rds'))
}

###############################################################################
## Define Routines

add_zeros <- function(x, new_genes) {
  x_names <- rownames(x)
  x <- rbind(x, matrix(0, nrow = length(new_genes), ncol = ncol(x)))
  rownames(x) <- c(x_names, new_genes)
  return(x)
}



Wang_preprocessing_routine <- function(input_dir, files_use, output_dir, genesets, BM.mapping, prefix) {
  print('combining data')
  print(Sys.time())
  first_loaded <- TRUE
  pb <- txtProgressBar(min = 0, max = length(files_use), style = 3)
  for (i in 1:length(files_use)) {
    f <- files_use[i]
    loaded_mat <- as.matrix(read.delim(file = file.path(input_dir, f), header = TRUE, row.names = 1))
    base_filename <- basename(f)
    SampleID <- regmatches(base_filename, regexpr('^GSM[0-9]*', base_filename))
    meta_data <- rbind(meta_data, data.frame(CellID = paste0(SampleID, colnames(loaded_mat)),
                                             SampleID = rep(SampleID, ncol(loaded_mat))))
    if (first_loaded) {
      first_loaded <- FALSE
      all_genes <- rownames(loaded_mat)
      combined_mat <- loaded_mat
    } else {
      all_genes <- union(rownames(loaded_mat), all_genes)
      diff_genes_1 <- setdiff(all_genes, rownames(loaded_mat))
      diff_genes_2 <- setdiff(all_genes, combined_mat)
      # set zeros for genes where no counts detected
      loaded_mat <- add_zeros(loaded_mat, diff_genes_1)
      combined_mat <- add_zeros(combined_mat, diff_genes_2)
      combined_mat <- cbind(combined_mat[all_genes,], loaded_mat[all_genes, ])
    }
    rm(loaded_mat)
    gc(full = TRUE)
    setTxtProgressBar(pb, i)
  }
  rownames(meta_data) <- colnames(combined_mat) <- meta_data$CellID
  rownames(combined_mat) <- all_genes
  # Run Seurat Pipeline
  seurat_obj <- seurat_subroutine(combined_mat, meta_data)
  rm(combined_mat)
  rm(meta_data)
  gc(full = TRUE)
  # Run TSNE
  seurat_obj <- Seurat::RunTSNE(seurat_obj)
  # Run Clustering
  seurat_obj <- Seurat::FindClusters(seurat_obj, force.recalc = TRUE, print.output = FALSE)
  seurat_obj@meta.data$cluster <- seurat_obj@meta.data$res.0.8
  # scoring
  seurat_obj <- scoring_subroutine(seurat_obj, genesets, preproc_data_dir, paste0(prefix, '_full'))
  # calculate z-scored avg chromosome expression
  chr_avg_output_list <- calcAvgChrMat(seurat_obj, BM.mapping, chr.use = as.character(1:22))
  chr_mat <- chr_avg_output_list$output.mat
  rownames(chr_mat) <- paste0('chr.', rownames(chr_mat))
  seurat_obj <- Seurat::AddMetaData(seurat_obj, as.data.frame(t(chr_mat)))
  # write.csv(chr_avg_output_list$chr.mapping, file = file.path(output_dir, paste0(prefix, '_full_chr_mapping.csv')))
  # write.csv(chr_avg_output_list$chr.summary, file = file.path(output_dir, paste0(prefix, '_full_chr_summary.csv')))
  # seurat_obj <- identify.glioma(seurat_obj)
  # # save results
  print('Saving Full Data')
  print(Sys.time())
  saveRDS(seurat_obj, file = file.path(output_dir, paste0(prefix, '_full_seurat.rds')))
  # ## Rerun Seurat pipeline on just glioma cells
  # # Run Seurat Pipeline
  # obj_raw_data <- seurat_obj@raw.data
  # obj_meta_data <- seurat_obj@meta.data
  # rm(seurat_obj)
  # gc(full = TRUE)
  # glioma_cells <- rownames(obj_meta_data)[obj_meta_data$is.glioma == 'glioma']
  # obj_raw_data <- obj_raw_data[,glioma_cells]
  # obj_meta_data <- obj_meta_data[glioma_cells,]
  # gc(full = TRUE)
  # seurat_obj_glioma <- seurat_subroutine(obj_raw_data, obj_meta_data)
  # # Run TSNE
  # seurat_obj_glioma <- Seurat::RunTSNE(seurat_obj_glioma)
  # # Run Clustering
  # seurat_obj_glioma <- Seurat::FindClusters(seurat_obj_glioma, force.recalc = TRUE, print.output = FALSE)
  # seurat_obj_glioma@meta.data$cluster <- seurat_obj_glioma@meta.data$res.0.8
  # # scoring
  # seurat_obj_glioma <- scoring_subroutine(seurat_obj_glioma, genesets, preproc_data_dir, paste0(prefix, '_glioma'))
  # print('Saving Filtered Data')
  # print(Sys.time())
  # saveRDS(seurat_obj_glioma, file = file.path(output_dir, paste0(prefix, '_glioma_seurat.rds')))
  # rm(seurat_obj_glioma)
  # gc(full = TRUE)
  # print('Finished')
  # print(Sys.time())
}


###############################################################################

## Do for snRNA-seq samples
snRNA_samples <- matrix_files[1:10]
# snRNA_samples <- matrix_files[1:2]
print('snRNA samples')
basename(snRNA_samples)
Wang_preprocessing_routine(input_dir = raw_data_dir, 
                           files_use = snRNA_samples, 
                           output_dir = preproc_data_dir, 
                           genesets = genesets, 
                           BM.mapping = BM.mapping, 
                           prefix = 'Wang_snRNA')
gc(full = TRUE)
## Do for scRNA-seq samples
scRNA_samples <- matrix_files[11:length(matrix_files)]
# scRNA_samples <- matrix_files[11:12]
print('snRNA samples')
basename(scRNA_samples)
Wang_preprocessing_routine(input_dir = raw_data_dir, 
                           files_use = scRNA_samples, 
                           output_dir = preproc_data_dir, 
                           genesets = genesets, 
                           BM.mapping = BM.mapping, 
                           prefix = 'Wang_scRNA')
