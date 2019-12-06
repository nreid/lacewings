#!/bin/bash
#SBATCH --job-name=process_radtags
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
INDIR=../../results

POOL=$INDIR/pool.fq.gz
BARCODES=../../metadata/henrysamples.txt

# make demultiplexed directory if it doesn't exist
OUTDIR=$INDIR/demultiplexed_fastqs
mkdir -p $OUTDIR

echo demultiplexing file $FASTQ using barcode set $BC

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



