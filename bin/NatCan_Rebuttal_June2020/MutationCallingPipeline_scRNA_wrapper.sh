#!/bin/bash
#SBATCH -t 24:00:00
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
### 1) Run bamslice from cellranger-dna (v1.1.0)
### Splits up bams by cell barcode (ExtractCellBarcodes_Seurat.r)
##############################################################

##############################################################
### EXAMPLE EXECUTION ON H4H
###
###
##############################################################


##############################################################
### 1) Extract and split barcodes from Seurat Object
###    'ExtractCellBarcodes_Seurat.r'
##############################################################




##############################################################
### 2) Run CellRanger-DNA bamslice on barcode chunks
###    'ExtractCellBarcodes_Seurat.r'
##############################################################


###clean up directories

##############################################################
### 3) Extract and split barcodes from Seurat Object
##############################################################
### -in = output from bamslice
### sbatch lots of mini scripts
bamtools split -tag CB -in G523_L.bam
