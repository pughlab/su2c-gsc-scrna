#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=10G
#SBATCH -p all
#SBATCH -c 30
#SBATCH -N 1
#SBATCH --account=pughlab

##############################################################
#             Mutation Calling pipeline in scRNA data        #
#                         L.Richards                         #
#                         July 2020                          #
##############################################################
### Reference: Aiden from Gary Bader's Lab
### Example execution on H4H:
### sbatch /cluster/home/lrichard/github/SU2C_GSC_scRNA/bin/NatCan_Rebuttal_June2020/MutationCalling_scRNA.sh /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/bams/BT147_L.possorted_genome_bam.bam

module load gatk/4.0.5.1
module load samtools/1.10


##############################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Sort bam with samtools
### 2) GATK AddOrReplaceReadGroups, MarkDuplicates, SplitNCigarReads
### 3) Run HaplotypeCaller (outputs vcf)
##############################################################

start=`date +%s`
echo ""
echo "********************"
echo "Set up variables"
date
echo "********************"
### path to scRNA bam outputted by CellRanger
### ie. /cluster/projects/pughlab/projects/BTSCs_scRNAseq/Manuscript_G607removed/NatCan_Rebuttal/MutationCalling/BT147_L.possorted_genome_bam.bam
bam=$1
### input sample name; this will be appended to out files
inp=$(basename $bam .bam)
### path to reference genome
#ref=/cluster/projects/pughlab/references/cellranger_10x/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
ref=$2

echo $inp
echo $bam
echo $ref
### set up outs directory
#mkdir $inp
#cd $inp

echo ""
echo "********************"
echo "Sort bam with samtools"
date
echo "********************"
samout="${inp}_temp_samsorted.bam"
samtools sort $bam -o $samout


echo ""
echo "********************"
echo "GATK AddOrReplaceReadGroups (Picard)"
date
echo "********************"
ARGout="${inp}__temp_sort_addgroup.bam"
gatk AddOrReplaceReadGroups --INPUT $samout --OUTPUT $ARGout --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 --SORT_ORDER coordinate


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


echo ""
echo "********************"
echo "GATK SplitNCigarReads"
date
echo "********************"
SNCout="${inp}_temp_sort_addgroup_dedupped_splitN.bam"
gatk SplitNCigarReads -R $ref -I $MDout -O $SNCout
rm $MDmetrix
rm $MDout


echo ""
echo "********************"
echo "Run HaplotypeCaller"
date
echo "********************"
#varfile="vars/${inp}_gatk_variants.vcf"
varfile="${inp}_gatk_variants.vcf"
gatk HaplotypeCaller -R $ref -I $SNCout --dont-use-soft-clipped-bases  --standard-min-confidence-threshold-for-calling 20.0 -O $varfile
rm $SNCout
rm *_splitN.*


echo ""
echo "********************"
echo "END"
echo "Duration: $((($(date +%s)-$start)/60)) minutes"
echo "********************"
