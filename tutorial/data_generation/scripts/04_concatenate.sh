#!/bin/bash 
#SBATCH --job-name=concatenate
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

INDIR=../results/prepended_fastqs

FASTQS=($(ls -1 $INDIR/*fq.gz))

cat ${FASTQS[@]} >pool.fq.gz

