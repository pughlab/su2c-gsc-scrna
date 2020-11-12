## Add Laura Richards' AUCell scores to BTSC June 2019 data. All data was downloaded from UHN H4H cluster
## See README in data directory for details.
## We do this so that we're using a consistent set of scores for analyses for manuscript.

library(Seurat)
library(SummarizedExperiment)

## load scRNA data
data_dir <- '../../data/raw/scRNA/Data/GSCs_June2019/' 
data_file <- 'Global_SU2C_BTSCs_CCregressed_noRibo.Rdata'

load(file.path(data_dir, data_file))

## Add Laura's calculated AUC scores to metadata
scores_file <- "GSC_AllSamples_AUCscores.rds"
AUC_scores <- readRDS(file.path(data_dir, scores_file))
stopifnot(all(rownames(BTSC@meta.data) == rownames(AUC_scores)))
BTSC@meta.data <- cbind(BTSC@meta.data, AUC_scores)

## Save Data
preproc_data_dir <- '../../data/preprocessed/scRNA/GSCs_June_2019'
if (!dir.exists(preproc_data_dir)) {
  dir.create(preproc_data_dir)
}
saveRDS(BTSC, file = file.path(preproc_data_dir, 'laura_SU2C_GSCs_June_2019_scored.rds'))
Sys.time()
sessionInfo()