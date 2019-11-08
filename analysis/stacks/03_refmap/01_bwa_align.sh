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
#SBATCH --array=[0-137]%20


hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9

# input, output files and directories
INDIR=../results/trimmed_data

mkdir -p ../results/trimmed_aligned
OUTDIR=../results/trimmed_aligned

# indexed reference genome
# used in script: /home/FCAM/ktaylor/rad_phy/scripts/bwa_array.sh
REFERENCE=/home/FCAM/ktaylor/rad_phy/reference_genomes/redundans_metaquast_filtered.nomt.masked.fasta

# note that 3 samples were dropped (no sequences)
FASTQS=($(ls -1 $INDIR/*.fastq.gz))

INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
OUTFILE=$(echo $INFILE | sed 's/fastq.gz/bam/')

SAM=$(echo $OUTFILE | sed 's/\..*//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)


bwa mem -t 4 -R $RG $REFERENCE $INDIR/$INFILE | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

date
