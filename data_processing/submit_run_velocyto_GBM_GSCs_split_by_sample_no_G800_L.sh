#!/bin/bash
#SBATCH --job-name=submit_velocyto
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=owen.whitley@mail.utoronto.ca
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o submit_run_velocyto_GBM_GSCs_split_by_sample_no_G800_L.out
# submit velocyto jobs for samples with merged files

# setup directories + options
top_dir='/cluster/projects/su2c_csc/scRNAseq_analysis/scRNA_seq_GBM_tissue/velocyto_tissue'
tumour_bams_dir="${top_dir}/../scRNA_seq_GBM_tissue_bams"
gsc_bams_dir="${top_dir}/../../scRNA_seq_GSC/laura_GSC_scRNA_bams"
barcodes_dir="${top_dir}/sample_barcodes_GBM_GSCs_split_by_sample_no_G800_L"
output_dir="${top_dir}/velocyto_GBM_GSCs_split_by_sample_no_G800_L"
if [ ! -d $output_dir ];
then
	mkdir $output_dir
fi
gtf_file="${top_dir}/gtf_file/gene_annotations.gtf"

# get tumour bams
tumour_bams=$(find $tumour_bams_dir | grep -E '.bam$')
# get GSC line bams
gsc_bams=$(find $gsc_bams_dir | grep -E '.bam$')
# setup bams
all_bams=$(echo $tumour_bams $gsc_bams)
# echo $all_bams
#
# merged_file_dirs=$(find ${top_dir} -type d | grep merged_files)

# get barcode files
barcode_files=$(find $barcodes_dir | grep tsv)
#echo "$barcode_files"
rerun=1

# function for performing basic calculations
calc() { awk "BEGIN{print $*}"; }

# for each directory containing merged files, get sample name + 
# appropriate files and submit a job to run velocyto



for bam_file in $all_bams
do
	echo "loop iteration: ${bam_file}"
	sample_name=$(echo $bam_file  | grep -oE '(G|BT)[0-9]*-*[A-Z]*_[LT]')
	new_dir="${output_dir}/${sample_name}_velocyto"
	
	if [ -d $new_dir ];
	then
	# if rerunning everything, remove the output directory
	# and make a new one
		if [ $rerun -eq 1 ];
		then 
			rm -rf $new_dir
			mkdir $new_dir
		fi
	else
	# if output directory doesn't exist, make a new one
		mkdir $new_dir
	fi

	# get barcodes_file
	barcode_file=$(echo "$barcode_files" | grep "${sample_name}_barcodes.tsv")
	num_barcode_files=$(echo "$barcode_files" | grep "${sample_name}_barcodes.tsv" | wc -l)
	if [ $num_barcode_files -ne 1 ];
	then
		echo "more or less than 1 barcode file for ${sample_name}, continuing loop"
		echo $barcode_file	
		echo $num_barcode_files
		continue
	fi
	bam_file_size=$(ls -alhr $bam_file | grep -oE '[0-9]*G' |  grep -oE '[0-9]*')
	# max vmem already set at 30g
	# walltime, assuming 6 hours for 20G and linear scaling of runtime
	walltime=$(calc "${bam_file_size}*6/20 + 12")
	walltime=$(echo $walltime | grep -oE '^[0-9]*')
	walltime="${walltime}:00:00"

	# submit job
	echo "submitting velocyto job for ${sample_name}"
	sbatch --export=barcodes_file=$barcode_file,bam_file=$bam_file,gtf_file=$gtf_file,new_dir=$new_dir,do_sort=true -o "${sample_name}_run_velocyto.out"  -t $walltime /cluster/home/hpc_owenwhitley/scripts/run_velocyto.sh		

done

exit 0
