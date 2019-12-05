#!/bin/bash 
#SBATCH --job-name=subsample
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

module load seqkit/0.10.0

INDIR=../results/prepended_fastqs

FULLPOOL=prepool.fq.gz

PROP=0.1

seqkit sample -p $PROP $INDIR/$FULLPOOL | gzip >$INDIR/pool.fq.gz
