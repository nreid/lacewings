#!/bin/bash 
#SBATCH --job-name=ipyrad_denovo
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=noah.reid@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load ipyrad/0.9.19


# create population and barcodes files from existing files in case they don't exist:

BC=../../metadata/henrysamples.txt
BCP=../../metadata/henrysamples_ipyrad.txt
POP=../../metadata/popmap.ipyrad.txt

paste <(cut -f 2 $BC) <(cut -f 1 $BC) >$BCP
paste <(cut -f 2 $BC) <(cut -f 2 $BC | sed 's/^....//') >$POP
cut -f 2 ../../metadata/popmap.txt | sort | uniq | sed 's/$/:0 /' | tr "\n" " " | sed 's/^/# /' >>$POP

# create output directory

mkdir -p ../../results/ipyrad/denovo

# run ipyrad

ipyrad -p params-henryetal_denovo.txt -s 1234567 -c 32 -r 

