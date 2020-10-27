#!/bin/bash
Rscript -e "rmarkdown::render('Neftel_2019_scRNA_determine_glioma.Rmd', output_dir = 'reports/Neftel_2019')"
Rscript Neftel_2019_explore_filtered_data.R > Neftel_2019_explore_filtered_data_out.txt 2>Neftel_2019_explore_filtered_data_stderr.txt
