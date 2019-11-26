#!/bin/bash 
#SBATCH --job-name=refmap.pl
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.41

# input, output files, directories

INDIR=../results/trimmed_aligned
OUTDIR=../results/stacks_refmap
mkdir -p $OUTDIR

# refmap.pl -s option is broken. 
POPMAP=../../../metadata/popmap.txt
cat $POPMAP | sed 's/\s/.trim\t/' >${POPMAP}.trim
POPMAP=../../../metadata/popmap.txt.trim

ref_map.pl --samples $INDIR -s trim --popmap $POPMAP -o $OUTDIR -T 10