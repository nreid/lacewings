#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20


hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9

# input, output files and directories
INDIR=../../data/trimmed_data

OUTDIR=../../data/aligned
mkdir -p $OUTDIR

# indexed reference genome
REFERENCE=../../../genome/redundans_metaquast_filtered.nomt.masked.fasta

FASTQS=($(ls -1 $INDIR/*.fastq.gz))

INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
OUTFILE=$(echo $INFILE | sed 's/fastq.gz/bam/')

SAM=$(echo $OUTFILE | sed 's/\..*//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

bwa mem -t 4 -R $RG $REFERENCE $INDIR/$INFILE | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

date
