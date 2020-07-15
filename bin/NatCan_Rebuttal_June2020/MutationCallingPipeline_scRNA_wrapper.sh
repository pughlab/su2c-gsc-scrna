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
### Reference: https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/using/bamslice
###
### Example execution on H4H:
### sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/CellRangerDNA_bamslice.sh /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/bams/BT147_L.possorted_genome_bam.bam

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Extract cell barcodes from seurat object
### 2) Use bamslice to make 1 bam per cell barcode
### 3) Run Haplotyper on each bam
##############################################################

module load R/3.6.1
module load cellranger-dna/1.1.0

##############################################################
### EXAMPLE EXECUTION ON H4H
###
##############################################################

SEURAT_OBJ=/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/Broad_Portal/seuratObjs/G523_L_res.0.2.RData
SAMPLE_ID=G523_L
BAM_FILE=/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/bams/G523_L.possorted_genome_bam.bam

WORKING_DIR=$(pwd)

##############################################################
### 1) Extract and split barcodes from Seurat Object
###    'ExtractCellBarcodes_Seurat.r'
##############################################################

mkdir $SAMPLE_ID
mkdir $SAMPLE_ID/cell_barcodes
cd ./$SAMPLE_ID/cell_barcodes/

Rscript /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/ExtractCellBarcodes_Seurat.r --seurat.obj $SEURAT_OBJ --outFilePrefix $SAMPLE_ID
wait


##############################################################
### 2) Run CellRanger-DNA bamslice on barcode chunks
###    'ExtractCellBarcodes_Seurat.r'
##############################################################

cd $WORKING_DIR/$SAMPLE_ID
mkdir bamslice
cd bamslice

CB_CSV=$WORKING_DIR/$SAMPLE_ID/cell_barcodes/${SAMPLE_ID}_BamSlice_Config.csv

sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/CellRangerDNA_bamslice.sh $SAMPLE_ID $CB_CSV $BAM_FILE

##############################################################
### 3) Run mutation calling pipeline
##############################################################
