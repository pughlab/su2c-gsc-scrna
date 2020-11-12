#!/bin/bash
# Run analyses pertaining to Nuc Seq data (Tumours + Lines) from Pugh Lab
Rscript Nuc_Seq_Explore_master.R > Nuc_Seq_Explore_master_out.txt 2>Nuc_Seq_Explore_master_stderr.txt
# Explore results of integration methods. Must run corresponding scripts called in Nuc_Seq_preprocess.sh
# Note that scripts to run MNN + Liger must be uncommented in that script and below lines must be uncommented as well
# Rscript -e "rmarkdown::render('Nuc_Seq_MNN_evaluate.Rmd', output_dir = 'reports/Nuc_Seq_MNN/')" > Nuc_Seq_MNN_evaluate_out.txt 2>Nuc_Seq_MNN_evaluate_stderr.txt
# Rscript -e "rmarkdown::render('Nuc_Seq_Liger_Explore.Rmd', output_dir = 'reports/Nuc_Seq_Liger_Explore.Rmd')" > Nuc_Seq_Liger_explore_out.txt 2>Nuc_Seq_Liger_explore_stderr.txt
