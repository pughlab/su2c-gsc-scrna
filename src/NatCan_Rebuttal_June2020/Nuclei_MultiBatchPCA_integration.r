

### https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html#6_other_utilities
pca.out <- multiBatchPCA(A=sce1,
                         B=sce2,
                         d= 10, ## keep 10 PCs
                         get.all.gene = TRUE,
                         subset.row=chosen.hvgs,
                         BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE)
                       )
names(pca.out)
