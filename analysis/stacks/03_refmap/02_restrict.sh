#!/bin/bash 
#SBATCH --job-name=restrict
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --qos=general
#SBATCH --partition=general


# software load
module load R/3.6.0

# reference genome
REFERENCE=/home/FCAM/ktaylor/rad_phy/reference_genomes/redundans_metaquast_filtered.nomt.masked.fasta

#output directory
OUTDIR=../../../metadata

# find restriction sites
cat $REFERENCE | Rscript restriction.R > $OUTDIR/restriction_sites.bed
