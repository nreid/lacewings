#!/bin/bash
#SBATCH --job-name=prepend
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[1-33]


echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# set input, output directories, variables, files
INDIR=../results/demultiplexed_fastqs
OUTDIR=../results/prepended_fastqs
mkdir -p $OUTDIR

BARCODE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../../metadata/henrysamples.txt | cut -f 1)
BQ=JJJJJ
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../../metadata/henrysamples.txt | cut -f 2)
FASTQ=${SAMPLE}.fq.gz

# prepend sequences using awk. 
	# if line # + 2 divisible by 4, add barcode. 
	# if line # divisible by 4, add base qualities
	# else, print line. 
	# gzip output, place in new directory
zcat $INDIR/$FASTQ | \
awk -v bc=$BARCODE -v bq=$BQ '{if((NR+2)%4 == 0) print bc $1; else if(NR%4 == 0) print bq $1; else print $1}' | \
gzip >$OUTDIR/$FASTQ

