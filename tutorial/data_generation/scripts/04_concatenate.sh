#!/bin/bash 
#SBATCH --job-name=concatenate
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
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

# concatenate sequences
cat ${FASTQS[@]} >$INDIR/prepool.fq.gz

# shuffle sequences
zcat $INDIR/prepool.fq.gz | \
awk '{OFS="\t"; getline seq; \
                getline sep; \
                getline qual; \
                print $0,seq,sep,qual}' | \
shuf | \
awk '{OFS="\n"; print $1,$2,$3,$4}' | \
gzip >$INDIR/pool.fq.gz


