#!/bin/bash
#SBATCH --job-name=process_radtags_full
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general


echo "host name : " `hostname`

module load stacks/2.41

#input/output directories, supplementary files
POOL=../../data_generation/results/prepended_fastqs/prepool.fq.gz

BARCODES=../../metadata/henrysamples.txt

# make demultiplexed directory if it doesn't exist
OUTDIR=../../results/demultiplexed_fastqs_FULL
mkdir -p $OUTDIR

process_radtags \
-f $POOL \
-o $OUTDIR \
-b $BARCODES \
-i gzfastq \
-y gzfastq \
-e sbfI \
-c \
-q \
-t 145 \
-s 20 \
--adapter_1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
--adapter_2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_mm 2



