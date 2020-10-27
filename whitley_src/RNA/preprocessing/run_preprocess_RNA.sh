#!/bin/bash
# remove previously generated preprocessed data if it exists.
set -e
timestamp=$(date +"%Y_%b_%d_%H_%M_%S")
preproc_dir="../../../data/preprocessed/RNA/"
top_archive_preproc_dir="../../../data/preprocessed/archived"
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

# make a reports directory for most recent html reports.
# archive previous results in a timestamped directory
reports_dir="./reports"
timestamp=$(date +"%Y_%b_%d")
if [ -d $reports_dir ]; then
        archived_reports_dir="${reports_dir}_archived_${timestamp}"
        echo "moving $reports_dir contents to $archived_reports_dir"
        if [ -d $archived_reports_dir ]; then
                rm -rf $archived_reports_dir
        fi
        mv $reports_dir $archived_reports_dir
fi

mkdir $reports_dir

Rscript run_preprocess_RNA.R > run_preprocess_RNA_out.txt 2>run_preprocess_RNA_stderr.txt
