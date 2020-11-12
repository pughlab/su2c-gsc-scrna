#!/bin/bash

# # clean preprocessed data directory
# find ../../data/preprocessed/scRNA -mindepth 1 -maxdepth 1 | xargs -l rm -rf
set -e
preproc_dir="../../data/preprocessed/scRNA/"
top_archive_preproc_dir="../../data/preprocessed/archived"
timestamp=$(date +"%Y_%b_%d_%H_%M_%S")
if [ ! -d $top_archive_preproc_dir ]; then
        mkdir $top_archive_preproc_dir
fi

if [ -d $preproc_dir ]; then
        preproc_dir_basename=$(basename $preproc_dir)
        archived_preproc_dir="${top_archive_preproc_dir}/${preproc_dir_basename}_archived_${timestamp}"
        echo "moving $preproc_dir contents to $archived_preproc_dir"
        if [ -d $archived_preproc_dir ]; then
                rm -rf $archived_preproc_dir
        fi
        mv $preproc_dir $archived_preproc_dir
fi

mkdir $preproc_dir

# # make a reports directory for most recent html reports.
# # archive previous results in a timestamped directory
# reports_dir="./reports_preprocessing"
# timestamp=$(date +"%Y_%b_%d")
# if [ -d $reports_dir ]; then
#         archived_reports_dir="${reports_dir}_archived_${timestamp}"
#         echo "moving $reports_dir contents to $archived_reports_dir"
#         if [ -d $archived_reports_dir ]; then
#                 rm -rf $archived_reports_dir
#         fi
#         mv $reports_dir $archived_reports_dir
# fi
# 
# mkdir $reports_dir

# add scores to June 2019 data
Rscript laura_SU2C_GSCs_June_2019_add_scores.R > laura_SU2C_GSCs_June_2019_add_scores_out.txt 2>laura_SU2C_GSCs_June_2019_add_scores_stderr.txt

# convert files to rds

Rscript laura_SU2C_BTSC_TumourCells_Oct2019_no_G800_L_to_rds.R > laura_SU2C_BTSC_TumourCells_Oct2019_no_G800_L_to_rds_out.txt 2>laura_SU2C_BTSC_TumourCells_Oct2019_no_G800_L_to_rds_stderr.txt

Rscript laura_SU2C_BTSC_TumourCells_Oct2019_to_rds.R > laura_SU2C_BTSC_TumourCells_Oct2019_to_rds_out.txt 2>laura_SU2C_BTSC_TumourCells_Oct2019_to_rds_stderr.txt

# # Nuc Seq Data (Pugh Lab provided combined set of data produced by Nuc Seq Technology + Live scRNA Data)
# ## Preprocess Nuc Seq Data alone
# bash Nuc_Seq_preprocess.sh > Nuc_Seq_preprocess_out.txt 2>Nuc_Seq_preprocess_stderr.txt

# Wang 2019 preprocessing
Rscript Wang_2019_preprocess.R > Wang_2019_preprocess_out.txt 2>Wang_2019_preprocecss_stderr.txt
# Neftel 2019 preprocessing
Rscript Neftel_2019_10X_preprocess.R > Neftel_2019_10X_preprocess_out.txt 2>Neftel_2019_10X_preprocess_stderr.txt
# Bhaduri 2020 preprocessing
Rscript Bhaduri_2020_preprocess.R > Bhaduri_2020_preprocess_out.txt 2>Bhaduri_2020_preprocess_stderr.txt
# RNA velocity preprocessing
bash run_scvelo_preprocess.sh > run_scvelo_preprocess_out.txt 2>run_scvelo_preprocess_stderr.txt
## Try Data Integration Schemes to integrate Nuc Seq/Live Data (comment out if you don't want to run)
# bash Nuc_Seq_Liger.sh
# bash Nuc_Seq_MNN.sh
# bash Nuc_Seq_ComBAT.sh
