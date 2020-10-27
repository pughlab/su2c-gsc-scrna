###############################################################################
### Convert Laura's Combined Tumour + Line file to rds

load('../../data/raw/scRNA/Data/SU2C_BTSC_TumourCells_combined_Oct_2019/BTSC_TumourCells_AUCell_Seurat.Rdata')
preproc_dir <- '../../data/preprocessed/scRNA/GSCs_Tumour_combined'
if (!dir.exists(preproc_dir)) {
  dir.create(preproc_dir)
}
source('./laura_SU2C_GBM_GSCs_combined_explore_helpers.R')
# clean/add metadata
unique(BTSC_TumourCells@meta.data$SampleID)
# [1] "BT127_L"   "BT147_L"   "BT48_L"    "BT67_L"    "BT73_L"    "BT84_L"    "BT89_L"    "BT94_L"   
# [9] "G523_L"    "G549_L"    "G564_L"    "G566_L"    "G583_L"    "G620_L"    "G637_L"    "G729_L"   
# [17] "G797_L"    "G799_L"    "G800_L"    "G837_L"    "G851_L"    "G876_L"    "G885_L"    "G895_L"   
# [25] "G945-I_L"  "G945-J_L"  "G945-K_L"  "G946-J_L"  "G946-K_L"  "G1003-A_T" "G1003-B_T" "G1003-C_T"
# [33] "G1003-D_T" "G620_T"    "G910-A_T"  "G910-B_T"  "G910-C_T"  "G910-D_T"  "G910-E_T"  "G945-I_T" 
# [41] "G945-J_T"  "G945-K_T"  "G946-I_T"  "G946-J_T"  "G946-K_T"  "G967-A_T"  "G967-B_T"  "G967-C_T" 
# [49] "G967-D_T"  "G983-A_T"  "G983-B_T"  "G983-C_T" 
# condense sample ids for replicates of same sample
BTSC_TumourCells@meta.data$SampleID <- sub('-[A-Z]_', '_', BTSC_TumourCells@meta.data$SampleID)
# calculate C1/C2 diff
BTSC_TumourCells <- calc_feat_diff(BTSC_TumourCells, var1 = 'RNA.GSC.c1_AUC', var2 = 'RNA.GSC.c2_AUC',
                                  new_feat_name = 'C1_C2_diff')
# manually add culture_cond, sample type
Lab <- ifelse(grepl('^G', BTSC_TumourCells@meta.data$SampleID), 'Dirks', 'Weiss')
SampleType <- ifelse(grepl('_L$', BTSC_TumourCells@meta.data$SampleID), 'Line', 'Tumour')
BTSC_TumourCells@meta.data$SampleType <- SampleType
BTSC_TumourCells@meta.data$Lab <- Lab
culture_cond <- rep('Tumour', length(SampleType))
culture_cond[Lab == 'Dirks' & SampleType == 'Line'] <- 'Adherent'
culture_cond[Lab == 'Weiss' & SampleType == 'Line'] <- 'Sphere'
BTSC_TumourCells@meta.data$culture_cond <- culture_cond
BTSC_TumourCells@meta.data$CultureMethod <- BTSC_TumourCells@meta.data$culture_cond
# add seurat clusters
BTSC_TumourCells <- FindClusters(BTSC_TumourCells, 
                                 reduction.type = 'pca', 
                                 dims.use = 1:10, 
                                 resolution = 0.8)
BTSC_TumourCells@meta.data$seurat_clusters <- BTSC_TumourCells@meta.data$res.0.8
# add metadata to note if we're looking at 637 or 800
BTSC_TumourCells@meta.data$is_637_L_800_L <- rep('other', nrow(BTSC_TumourCells@meta.data))
BTSC_TumourCells@meta.data$is_637_L_800_L[BTSC_TumourCells@meta.data$SampleID == 'G800_L'] <- 'G800_L'
BTSC_TumourCells@meta.data$is_637_L_800_L[BTSC_TumourCells@meta.data$SampleID == 'G637_L'] <- 'G637_L'
# run tsne
BTSC_TumourCells <- RunTSNE(BTSC_TumourCells, dims.use = 1:10)
saveRDS(BTSC_TumourCells, file.path(preproc_dir, 'laura_BTSC_Tumour_combined_Oct_2019.rds'))
