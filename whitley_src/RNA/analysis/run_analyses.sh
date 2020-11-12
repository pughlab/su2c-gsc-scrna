#!/bin/bash
# kill on error
set -e
# make a reports directory for most recent html reports. 
# archive previous results in a timestamped directory 
reports_dir="./reports" 
timestamp=$(date +"%Y_%b_%d_%H_%M_%S") 
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
results_dir="../../../results/RNA"
top_archive_results_dir="../../../results/archived"

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

## RUN ANALYSES
# main analyses for Richards et al. 2019 manuscript
main_files=$(cat main_analysis_notebooks.txt)
# main_files=$(cat main_analysis_notebooks_last3.txt)
for p in $main_files; do
	parent_dir=$(dirname $p) 
	reports_dir="${parent_dir}/reports"
	Rscript -e "rmarkdown::render('${p}', output_dir = '${reports_dir}')"
done
# supplementary analyses
supplementary_files=$(cat supplementary_analysis_notebooks.txt)
 
for p in $supplementary_files; do
	parent_dir=$(dirname $p) 
 	reports_dir="${parent_dir}/reports"
 	Rscript -e "rmarkdown::render('${p}', output_dir = '${reports_dir}')"
done
# Rscript run_analyses.R > run_analyses_out.txt 2>run_analyses_stderr.txt
