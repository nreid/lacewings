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
#SBATCH --array=[0-2]


echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load stacks/2.41

#input/output directories, supplementary files
INDIR=/home/FCAM/ktaylor/rad_analysis/data

POOLS=($(ls -1 $INDIR/*fastq.gz))
BARCODES=($(ls -1 ../../../metadata/*barcode*))

# make demultiplexed directory if it doesn't exist
mkdir -p ../results/demultiplexed_fastqs
OUTDIR=../results/demultiplexed_fastqs

FASTQ=$(echo ${POOLS[$SLURM_ARRAY_TASK_ID]})
BC=$(echo ${BARCODES[$SLURM_ARRAY_TASK_ID]})

echo demultiplexing file $FASTQ using barcode set $BC

process_radtags \
-f $FASTQ \
-o $OUTDIR \
-b $BC \
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



