#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=150G
#SBATCH -p veryhimem
#SBATCH -c 30
#SBATCH -N 1
#SBATCH --account=pughlab

##############################################################
#             Split 10X Genomics Bam using bamslice          #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/using/bamslice
###
### Example execution on H4H:
### sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/CellRangerDNA_bamslice.sh /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/bams/BT147_L.possorted_genome_bam.bam

module load cellranger-dna/1.1.0

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Extract cell barcodes from Seurat Object
### 2)
### 3)
### 4)
### 5)
##############################################################

start=`date +%s`
echo ""
echo "********************"
echo "######"
date
echo "********************"

cellranger-dna bamslice --id=tumor_subsets \
                       --csv=tumor_subsets.csv \
                       --bam=/home/jdoe/runs/tumor_mixed/outs/possorted_bam.bam
