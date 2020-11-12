#!/bin/bash
Rscript -e "rmarkdown::render('laura_SU2C_GBM_Nov_2019_explore.Rmd', output_dir = 'reports/laura_SU2C_GBM_Nov_2019')"
Rscript -e "rmarkdown::render('laura_SU2C_GBM_Nov_2019_GSEA.Rmd', output_dir = 'reports/laura_SU2C_GBM_Nov_2019')"
Rscript -e "rmarkdown::render('laura_SU2C_GSCs_Nov_2019_explore.Rmd', output_dir = 'reports/laura_SU2C_GSCs_Nov_2019')"
Rscript -e "rmarkdown::render('laura_SU2C_GSCs_Nov_2019_GSEA.Rmd', output_dir = 'reports/laura_SU2C_GSCs_Nov_2019')"
