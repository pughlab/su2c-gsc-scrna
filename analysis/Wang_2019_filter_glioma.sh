#!/bin/bash
Rscript -e "rmarkdown::render('Wang_2019_snRNA_determine_glioma.Rmd', output_dir = 'reports/Wang_2019_filter_glioma')"
Rscript -e "rmarkdown::render('Wang_2019_scRNA_determine_glioma.Rmd', output_dir = 'reports/Wang_2019_filter_glioma')"
