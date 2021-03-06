### Run Cell Ranger (v2.0.0) on 10X Genomics fastqs
### L.Richards
### $1 is the sample prefix
### $2 is the prefex on the fastq
### Example Usage on Mordor: qsub -q highmem.q runCellRanger_v2.0.0.sh BT127_L BT127

#!/bin/bash
#
#$ -cwd

mkdir /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/CellRanger_v2.0.0/$1
cd /mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/CellRanger_v2.0.0/$1

CELLRANGER_ROOT='/mnt/work1/users/pughlab/bin/CellRanger/v2.0/cellranger-2.0.0/cellranger'
TOOL='count'

$CELLRANGER_ROOT $TOOL --id=$1 \
--sample=$2 \
--transcriptome=/mnt/work1/users/pughlab/references/refdata-cellranger-GRCh38-1.2.0 \
--fastqs=/mnt/work1/users/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/data/fastqs/$1 \
--expect-cells=3000 \
--jobmode=sge \
--mempercore=6 \
--maxjobs=40
