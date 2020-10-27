#!/bin/bash
# OVERVIEW:
#	run analyses by calling 'master' shell scripts. can run individual sets of analyses by running individual shell scripts, but they assume that the reports/results directories made in this script have already been set up. Some of the analyses assume that all analyses for RNA-seq have already been run. Will have to have similar result generation/archiving workflow done for RNA-seq analyses, and have GSEA comparison analysis pull from an archived RNA-seq analysis so that we can pull from earlier results while continuously generating 'clean' reruns of results.
# 	Also note that for the 'GBM/GSC combined' analyses, due to memory limitson machine used for project, will need to run individual analyses in master R script 1 at a time as R for some reason is not returning memory even when gc(full = TRUE) is called. Do this by commenting/uncommenting individual analyses in laura_SU2C_GBM_GSCs_combined_explore_master.R, and run laura_SU2C_GBM_GSCs_combined_explore_master.h for all analyses.

# make a reports directory for most recent html reports.
# archive previous results in a timestamped directory
set -e
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
results_dir="../../results/scRNA"
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

# run analyses. note that you can comment/uncomment things if you only want to run a subset. each of these analyses should be able to be run on their own provided RNA-seq preprocessing + analyses run previously, and scRNA preprocessing done
bash laura_SU2C_GBM_GSCs_combined_explore_master.sh
bash laura_SU2C_GBM_GSCs_Nov_2019_explore_master.sh
bash laura_SU2C_GBM_GSCs_GSEA_scripts_master.sh
# supplementary analyses used to address reviewers' questions for Richards et al. 2020 manuscript
bash run_bhaduri_full.sh > run_bhaduri_full_out.txt 2>run_bhaduri_full_stderr.txt
bash Neftel_2019_analyses.sh > Neftel_2019_analyses_out.txt 2>Neftel_2019_analyses_stderr.txt
bash Wang_2019_analyses.sh > Wang_2019_analyses_out.txt 2>Wang_2019_analyses_stderr.txt
Rscript laura_PCA_gene_selection.R > laura_PCA_gene_selection_out.txt 2>laura_PCA_gene_selection_stderr.txt
bash run_scvelo.sh > run_scvelo_out.txt 2>run_scvelo_stderr.txt
bash laura_SU2C_GBM_GSCs_combined_explore_varimax_PC1_PC2_master.sh > laura_SU2C_GBM_GSCs_combined_explore_varimax_PC1_PC2_master_out.txt 2>laura_SU2C_GBM_GSCs_combined_explore_varimax_PC1_PC2_master_stderr.txt
Rscript laura_cell_cycle_regress_plots.R > laura_cell_cycle_regress_plots_out.txt 2>laura_cell_cycle_regress_plots_stderr.txt
