#!/bin/bash

##############################################################
#             Mutation Calling pipeline in scRNA data        #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: Aiden from Gary Bader's Lab

##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Sort bam with samtools
### 2) Run mnnCorrect from Batchelor package
### 3) Recore corrected matrix with Dev and IR
### 4) Save results
##############################################################

module load samtools/1.10
module load gatk/4.0.5.1

start=`date +%s`
echo "********************"
echo "Set up variables"
date
echo "********************"
echo ""
### path to scRNA bam outputted by CellRanger
### ie. BT147_L.possorted_genome_bam.bam
#inp=$1
#bam=$1
bam=/cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/BT147_L.possorted_genome_bam.bam
echo $bam

### input sample name
### this will be appended to out files
inp=$(basename $bam .possorted_genome_bam.bam)
echo $inp

### path to reference genome
ref=/cluster/projects/pughlab/references/cellranger_10x/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
echo $ref

### set up directory schema
mkdir $inp
cd $inp


echo "********************"
echo "Sort bam with samtools"
date
echo "********************"
echo ""
samout="${inp}_temp_samsorted.bam"
samtools sort $bam -o $samout


echo "********************"
echo "GATK AddOrReplaceReadGroups (Picard)"
date
echo "********************"
echo ""
ARGout="${inp}__temp_sort_addgroup.bam"
gatk AddOrReplaceReadGroups --INPUT $samout  --OUTPUT $ARGout --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 --SORT_ORDER coordinate



echo "********************"
echo "GATK MarkDuplicates (Picard)"
date
echo "********************"
echo ""
MDout="${inp}_temp_sort_addgroup_rawdedupped.bam"
MDmetrix="${inp}_temp_sort_addgroup_rawdedupped_metrics.bam"
gatk MarkDuplicates --INPUT $ARGout --OUTPUT $MDout --METRICS_FILE $MDmetrix
rm $samout
rm $ARGout


echo "********************"
echo "GATK SplitNCigarReads"
date
echo "********************"
echo ""
SNCout="${inp}_temp_sort_addgroup_dedupped_splitN.bam"
gatk SplitNCigarReads -R $ref -I $MDout -O $SNCout
rm $MDmetrix
rm $MDout


echo "********************"
echo "Run HaplotypeCaller"
date
echo "********************"
echo ""
varfile="vars/${inp}_gatk_variants.vcf"
gatk HaplotypeCaller -R $ref -I $SNCout --dont-use-soft-clipped-bases  --standard-min-confidence-threshold-for-calling 20.0 -O $varfile
rm $SNCout
rm *_splitN.*


echo "********************"
echo "END"
echo "Duration: $((($(date +%s)-$start)/60)) minutes"
echo "********************"
echo ""
