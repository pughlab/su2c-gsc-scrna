###############################################################################
## Run preprocessing

library(knitr)
library(rmarkdown)
file.list <- c('owen_RNA_preprocess_coh1234_nonSU2C.Rmd', 'owen_RNA_preprocess_adult_GBM_coh1234_nonSU2C.Rmd', 'owen_RNA_preprocess_all_GBM_coh1234_nonSU2C.Rmd', 'owen_RNA_preprocess_adult_GBM_line_Dirks_coh1234_nonSU2C.Rmd')
for (fname in file.list) {
  rmarkdown::render(input = fname, output_dir = './reports', envir = new.env())
}
