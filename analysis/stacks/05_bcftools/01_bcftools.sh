#!/bin/bash 
#SBATCH --job-name=bcftools
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load bcftools/1.9
module load htslib/1.9

# input/output directories, supplemental files

INDIR=../results/trimmed_aligned

# bam list
find $INDIR -name "*bam" | sort >$INDIR/bams.list
BAMLIST=$INDIR/bams.list

OUTDIR=../results/bcftools
mkdir -p $OUTDIR

# ref genome
REFERENCE=../results/redundans_metaquast_filtered.nomt.masked.fasta

# set a variable for the sites targeted
TARGETS=../../../metadata/targets.bed


bcftools mpileup \
	-f $REFERENCE \
	-b $BAMLIST \
	-q 20 -Q 30 \
	-a "DP,AD"
	--targets-file $TARGETS | \
bcftools call -m -v -Oz -o $OUTDIR/chrysoperla.vcf.gz



