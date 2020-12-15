
library("Seurat", 
        lib ="/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript/lib/" )
library(AUCell)

print("Load Data")
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/SeuratObj/BTSCs/Global_SU2C_BTSCs_CCregressed_noRibo.Rdata")

print("Reset Ident")
BTSC <- SetIdent(BTSC, 
                  ident.use = factor(BTSC@meta.data$SampleID)
                 )

print(table(BTSC@ident))

print("Remove G800_L")
raw.data <- as.matrix(BTSC@raw.data)
dim(raw.data)

remove.cells <- WhichCells(object = BTSC,
                           ident = c("G800_L")
                          )

raw.data <- raw.data[ ,!colnames(raw.data) %in% remove.cells]
print(dim(raw.data))

print("Create new seurat object")

BTSC <- CreateSeuratObject(raw.data = raw.data)


print("Update metadata")
BTSC@meta.data$orig.ident <- paste0(BTSC@meta.data$orig.ident, "_L")

print("Calculate percent mito")

mito.genes <- grep(pattern = "^MT-", 
                   x = rownames(x = BTSC@data), 
                   value = TRUE
                  )
percent.mito <- Matrix::colSums(BTSC@raw.data[mito.genes, ])/Matrix::colSums(BTSC@raw.data)


BTSC <- AddMetaData(object = BTSC, 
                                metadata = percent.mito, 
                                col.name = "percent.mito"
                               )

print("Normlaize data")
BTSC <- NormalizeData(object = BTSC, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 10000
                                 )

print("Find variable genes")
BTSC <- FindVariableGenes(object = BTSC, 
                                      mean.function = ExpMean, 
                                      dispersion.function = LogVMR, 
                                      x.low.cutoff = 0.0125, 
                                      x.high.cutoff = 3, 
                                      y.cutoff = 0.5
                                     )

print("Calculate CC Difference")
cc.genes <- readLines(con = "/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GSC_Tumour_PCA/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

BTSC <- CellCycleScoring(object = BTSC, 
                                     s.genes = s.genes, 
                                     g2m.genes = g2m.genes, 
                                     set.ident = FALSE
                                    )

BTSC@meta.data$CC.Difference <- BTSC@meta.data$S.Score - BTSC@meta.data$G2M.Score
head(BTSC@meta.data)

#################
print("Scale the data")
BTSC <- ScaleData(object = BTSC, 
                              vars.to.regress = c("nUMI", "percent.mito", "CC.Difference")
                 )

print("Run PCA (exclude ribo genes)")

ribo.genes <- c(rownames(BTSC@data)[grep("^RP[1:9]", rownames(BTSC@data))],
                rownames(BTSC@data)[grep("^RP[L,S]", rownames(BTSC@data))]
               ) #2006

BTSC <- RunPCA(BTSC, 
               pc.genes = rownames(BTSC@data)[!rownames(BTSC@data) %in% ribo.genes],
               pcs.compute = 100,
               genes.print = 10
              )


print("Score data with AUCell")
print("Load gene signatures")

##load gene sigs
load("~/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/GBM_stemnessProgram/AUCell_signatures.Rdata")
str(sigs)

print("Isolate expression matrix")
exprMatrix <- as.matrix(BTSC@data)
exprMatrix[1:10, 1:10]
print(dim(exprMatrix))

print("Calculate cell rankings")
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=FALSE)

print("Calculate AUC")
cells_AUC <- AUCell_calcAUC(sigs, cells_rankings)

print("Explore thresholds")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                             plotHist=FALSE, 
                                             assign=TRUE) 

AUC <- t(as.data.frame(cells_AUC@assays$data$AUC))
colnames(AUC) <- paste0(colnames(AUC), "_AUC")
AUC <- data.frame(AUC)
print(AUC[1:2, 1:2])

print("Add AUC scores to metadata")

BTSC <- AddMetaData(BTSC, metadata = AUC)
print(head(BTSC@meta.data))

print("Save data")

#PCA data
pc <- BTSC@dr$pca
saveRDS(pc, file = "BTSC_G800L_removed_PCA.rds")

#metadata
meta <- BTSC@meta.data
saveRDS(meta, file = "BTSC_G800L_removed_AUC_metadata.rds")

##save seruat object 
saveRDS(BTSC, file = "BTSC_G800L_removed_AUCell_SeuratObj.rds")

