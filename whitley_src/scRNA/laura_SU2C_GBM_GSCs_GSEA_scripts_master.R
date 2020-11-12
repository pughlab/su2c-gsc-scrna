html_output_dir <- './laura_SU2C_GBM_GSCs_GSEA_scripts'
rmarkdown::render('laura_SU2C_GSCs_Nov_2019_GSEA.Rmd', output_dir = html_output_dir, envir = new.env())
gc(full = TRUE)
rmarkdown::render('laura_SU2C_GBM_Nov_2019_GSEA.Rmd', output_dir = html_output_dir, envir = new.env())
gc(full = TRUE)
