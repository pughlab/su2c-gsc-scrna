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

##############################################################
### 1) Set up parameters and load modules
##############################################################

module load samtools/1.10
module load picard/2.10.9

inp=$1 ### scRNA bam outputted by CellRanger


##############################################################
# 1) Sort bam file with Samtools
##############################################################
echo "********************"
echo "Sort bam with samtools"
date
echo "********************"

samout="${inp}_temp_samsorted.bam"
samtools sort $inp -o $samout
ARGout="${inp}__temp_sort_addgroup.bam"


echo "********************"
echo "GATK AddOrReplaceReadGroups (picard)"
date
echo "********************"

gatk AddOrReplaceReadGroups --INPUT $samout  --OUTPUT $ARGout --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 --SORT_ORDER coordinate
MDout="${inp}_temp_sort_addgroup_rawdedupped.bam"
MDmetrix="${inp}_temp_sort_addgroup_rawdedupped_metrics.bam"

echo "********************"
echo "GATK MarkDuplicates "
date 
echo "********************"

gatk MarkDuplicates --INPUT $ARGout --OUTPUT $MDout --METRICS_FILE $MDmetrix
rm $samout
rm $ARGout
SNCout="${inp}_temp_sort_addgroup_dedupped_splitN.bam"
gatk SplitNCigarReads -R ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I $MDout -O $SNCout
rm $MDmetrix
rm $MDout
varfile="vars/${inp}_gatk_variants.vcf"
gatk HaplotypeCaller -R ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I $SNCout --dont-use-soft-clipped-bases  --standard-min-confidence-threshold-for-calling 20.0 -O $varfile
rm $SNCout
rm *_splitN.*
