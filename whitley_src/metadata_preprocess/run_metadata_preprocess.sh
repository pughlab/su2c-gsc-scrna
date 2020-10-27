#!/bin/bash
# run assembly of SU2C metadata from files sent by collaborators. For details, see script run and read OVERVIEW, see README in folders for any files from /data/raw loaded in
set -e
preproc_dir="../../data/preprocessed/metadata"
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
Rscript SU2C_all_metadata_preprocess.R > SU2C_all_metadata_preprocess_out.txt 2>SU2C_all_metadata_preprocess_stderr.txt
