#!/bin/bash 
#SBATCH --job-name=populations
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


INDIR=../results/calocedrii_test

POPMAP=../../../metadata/popmap_calocedrii.txt

populations \
-P $INDIR \
-M $POPMAP \
--vcf \
--fasta-samples \
--fasta-loci \
--treemix \
-t 8