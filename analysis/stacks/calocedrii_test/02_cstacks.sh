#!/bin/bash 
#SBATCH --job-name=cstacks
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.41

# input, output files, directories
INDIR=../results/calocedrii_test

# doesn't need an output directory? 
# mkdir -p ../results/denovo_cstacks
# OUTDIR=../results/denovo_cstacks

# generate pop map if it doesn't exist
POPMAP=../../../metadata/popmap_calocedrii.txt

paste \
<(ls $INDIR | grep snps | sed 's/\..*//') \
<(ls $INDIR | grep snps | sed 's/\..*//' | cut -c 5-) | \
grep calocedrii \
>$POPMAP

cstacks -P $INDIR -M $POPMAP -p 20 --max-gaps 10 -n 15
