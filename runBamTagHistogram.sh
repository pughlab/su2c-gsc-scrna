### Run 'BamTagHistogram' on 10X Genomics Bam to obtain cell barcode readcounts
### Output used to create knee plot for Dropbead to determine cell count
### L.Richards
### NOTE: takes about 5-10min to run per library
### 

#!/bin/bash
#
#$ -cwd


### Example Usage:
### qsub runBamTagHistogram10X.sh SAMPLE_NAME
### qsub runBamTagHistogram10X.sh BT127_L

DROPSEQ_ROOT='/mnt/work1/users/pughlab/bin/Drop-seq_tools-1-2.12'
#BAM_FILE=$1
SAMPLE_NAME=$1

echo 'Running BamTagHistogram'
echo $(date)

TOOL='BAMTagHistogram'
$DROPSEQ_ROOT'/'$TOOL \
I=/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/data/outs/$SAMPLE_NAME/outs/possorted_genome_bam.bam \
O=$SAMPLE_NAME.cell_readcounts.txt.gz \
READ_QUALITY=0 \
TAG=CB

echo 'BamTagHistogram Finished'
echo $(date)

## END
