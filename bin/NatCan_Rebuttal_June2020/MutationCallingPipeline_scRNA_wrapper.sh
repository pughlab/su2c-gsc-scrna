#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --mem=60G
#SBATCH -p himem
#SBATCH -c 30
#SBATCH -N 1
#SBATCH --account=pughlab

##############################################################
#                scRNA Mutation Calling Pipeline             #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Extract cell barcodes from seurat object
### 2) Run bamslice to generate 1 bam per cell barcode
### 3) sbatch HaplotypeCaller script on each bam
##############################################################

module load R/3.6.1
module load cellranger-dna/1.1.0
module load samtools/1.10
module load gatk/4.0.5.1

##############################################################
### EXAMPLE EXECUTION ON H4H
### sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/MutationCallingPipeline_scRNA_wrapper.sh
###
###
###
###
##############################################################

#### DEVELPOMENT ####
SEURAT_OBJ=/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/G523_L_res.0.2.RData
SAMPLE_ID=G523_L
BAM_FILE=/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/bams/G523_L.possorted_genome_bam.bam
REF=/cluster/projects/pughlab/references/cellranger_10x/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa

WORKING_DIR=$(pwd)



##############################################################
### 1) Extract and split barcodes from Seurat Object
###    'ExtractCellBarcodes_Seurat.r'
##############################################################
echo '#####################################'
echo '1) Extract Cell Barcodes'
date
echo '#####################################'
start=`date +%s`
mkdir $SAMPLE_ID
mkdir $SAMPLE_ID/cell_barcodes
cd ./$SAMPLE_ID/cell_barcodes/

Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/ExtractCellBarcodes_Seurat.r --seurat.obj $SEURAT_OBJ --outFilePrefix $SAMPLE_ID
wait

echo "Duration: $((($(date +%s)-$start)/60)) minutes"



##############################################################
### 2) Run CellRanger-DNA bamslice on barcode chunks
###    'ExtractCellBarcodes_Seurat.r'
##############################################################
echo '#####################################'
echo '2) CellRangerDNA bamslice'
date
echo '#####################################'
start=`date +%s`
cd $WORKING_DIR/$SAMPLE_ID
mkdir bamslice
cd bamslice

CB_CSV=$WORKING_DIR/$SAMPLE_ID/cell_barcodes/${SAMPLE_ID}_BamSlice_Config.csv

bash /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/CellRangerDNA_bamslice.sh $SAMPLE_ID $CB_CSV $BAM_FILE
wait

#### add code to clean up big uneeded files from bam BamSlice

echo "Duration: $((($(date +%s)-$start)/60)) minutes"



##############################################################
### 3) sbatch mutation calling pipeline
##############################################################
### This creates tons of jobs...1 job per cell/bam
### Try submitting a few jobs with low memory (per cell bams are pretty small)
echo '#####################################'
echo '3) Mutation Calling Pipeline (HaplotypeCaller)'
date
echo '#####################################'

start=`date +%s`
cd $WORKING_DIR/$SAMPLE_ID
mkdir HaplotypeCaller
cd HaplotypeCaller

ls $WORKING_DIR/${SAMPLE_ID}/bamslice/${SAMPLE_ID}/outs/subsets/*bam > cellbams.txt

### here we need to sbatch lots of miniscripts
for CELL_BAM in $(cat cellbams.txt); do
    #echo $CELL_BAM
    sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/HaplotypeCaller_MutationCalling_scRNA.sh $CELL_BAM $REF
done


##############################################################
echo '#####################################'
echo 'End of pipeline'
date
echo '#####################################'
