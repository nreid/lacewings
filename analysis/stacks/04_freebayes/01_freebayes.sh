#!/bin/bash 
#SBATCH --job-name=freebayes
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load freebayes/1.1.0
module load htslib/1.9


# input/output directories, supplemental files

INDIR=../results/trimmed_aligned

# bam list
find $INDIR -name "*bam" | sort >$INDIR/bams.list
BAMLIST=$INDIR/bams.list

OUTDIR=../results/freebayes
mkdir -p $OUTDIR

# ref genome
REFERENCE=../results/redundans_metaquast_filtered.nomt.masked.fasta

# set a variable for the sites targeted
TARGETS=../../../metadata/targets.bed

freebayes -f $REFERENCE -L $BAMLIST -t $TARGETS -m 30 -q 20 -k -V | \
bgzip -c >$OUTDIR/chrysoperla.vcf.gz

tabix -p vcf $OUTDIR/chrysoperla.vcf.gz
