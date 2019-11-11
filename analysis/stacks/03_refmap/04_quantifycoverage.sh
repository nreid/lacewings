#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general


module load samtools/1.9
module load htslib/1.9
module load R/3.6.0
module load bedtools/2.27.1
module load bamtools/2.5.1


# set/make directories
INDIR=../results/trimmed_aligned

# make an output directory
OUTDIR=../results/trimmed_aligned_stats
mkdir -p $OUTDIR


# make a list of bam files
find $INDIR -name "*bam" >$INDIR/bams.list

bamtools merge -list $INDIR/bams.list | \
bamtools filter -in - -mapQuality ">30" | \
samtools depth -d 300000 /dev/stdin | \
bgzip > $OUTDIR/depthperbase.txt.gz


