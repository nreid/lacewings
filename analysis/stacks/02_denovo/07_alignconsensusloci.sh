#!/bin/bash 
#SBATCH --job-name=alignconsensusloci
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9
module load bedtools/2.29.0

# indexed reference genome
# used in script: /home/FCAM/ktaylor/rad_phy/scripts/bwa_array.sh
REFERENCE=/home/FCAM/ktaylor/rad_phy/reference_genomes/redundans_metaquast_filtered.nomt.masked.fasta

FASTA=../results/denovo_ustacks/populations.loci.fa
BAM=../results/denovo_ustacks/populations.loci.bam
BED=../results/denovo_ustacks/populations.loci.bed

bwa mem -t 4 $REFERENCE $FASTA | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$BAM - >$BAM

bedtools bamtobed -i $BAM >$BED