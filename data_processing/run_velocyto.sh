#!/bin/bash
#PBS -N velocyto
#PBS -l nodes=1:ppn=1

#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=owen.whitley@mail.utoronto.ca
#SBATCH -N 1
#SBATCH --mem=30G
# run velocyto for 1 sample

set -e
module load velocyto/0.17.13
module load samtools

barcodes_file=$barcodes_file
bam_file=$bam_file
gtf_file=$gtf_file
new_dir=$new_dir
threads=1
do_sort=$do_sort



# velocyto run -b $barcodes_file -o $new_dir -@ $threads $CB_sorted_bam $gtf_file

if [[ $do_sort = true ]];
then
	echo "sorting bamfile"
	date
	# use to sort if bam files not pre-sorted
	sorted_fname=$(echo $bam_file | xargs basename | sed s/.bam$/_sorted.bam/g)
	sorted_fname="${new_dir}/${sorted_fname}"
	samtools sort -m  "25G" -@ $threads -o $sorted_fname $bam_file 	
	bam_file=$sorted_fname
	date
	echo "sorting complete"
fi

echo "velocyto: running"
date

velocyto run -o $new_dir -@ $threads -b $barcodes_file $bam_file $gtf_file

echo "velocyto: finished"
date

exit 0
