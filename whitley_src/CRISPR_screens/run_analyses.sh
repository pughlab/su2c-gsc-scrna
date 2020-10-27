#!/bin/bash
# OVERVIEW:

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
# make a subdirectory of results for scRNA data if it doesn't exist already
# if it does exist, archive old results
results_dir="../../results/CRISPR_screens"
top_archive_results_dir="../../results/archived"

if [ ! -d $top_archive_results_dir ]; then
	mkdir $top_archive_results_dir
fi

if [ -d $results_dir ]; then
	results_dir_basename=$(basename $results_dir)
	archived_results_dir="${top_archive_results_dir}/${results_dir_basename}_archived_${timestamp}"
	echo "moving $results_dir contents to $archived_results_dir"
	if [ -d $archived_results_dir ]; then
		rm -rf $archived_results_dir
	fi
        mv $results_dir $archived_results_dir
fi

mkdir $results_dir

# run analyses. note that you can comment/uncomment things if you only want to run a subset. 
## GSEA analysis for TKOv3 library as used in Richards et al. 2020 manuscript
Rscript -e "rmarkdown::render('CRISPR_GSCs_avg_qBF_Dev_v_IR_GSEA.Rmd', output_dir = 'reports')" > CRISPR_GSCs_avg_qBF_Dev_v_IR_GSEA_out.txt 2 > CRISPR_GSCs_avg_qBF_Dev_v_IR_GSEA_stderr.txt
