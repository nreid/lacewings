#!/bin/bash 
#SBATCH --job-name=discards
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

module load seqkit/0.10.0

INDIR=../results/demultiplexed_fastqs

OUTDIR=../results/prepended_fastqs
mkdir -p $OUTDIR

OUTFILE=discard_samples.fq.gz

FASTQS=($(ls -1 $INDIR/pool*discards))

PROP=0.25

seqkit sample -p $PROP ${FASTQS[0]} | gzip >$$OUTDIR/$OUTFILE
seqkit sample -p $PROP ${FASTQS[1]} | gzip >>$$OUTDIR/$OUTFILE
seqkit sample -p $PROP ${FASTQS[2]} | gzip >>$$OUTDIR/$OUTFILE

