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
TRIMDIR=../results/trimmed_data/

POOLS=($(ls -1 $TRIMDIR/*fastq.gz))
BARCODES=($(ls ../../../metadata/ | grep barcodes))

# make demultiplexed directory if it doesn't exist
mkdir -p ../results/demultiplexed_fastqs
OUTDIR=../results/demultiplexed_fastqs

FASTQ=$(echo ${POOLS[$SLURM_ARRAY_TASK_ID]})
BC=$(echo ${POOLS[$SLURM_ARRAY_TASK_ID]})


process_radtags \
-f $FASTQ \
-o $OUTDIR \
-b $BC \
-i gzfastq \
-y gzfastq \
-e sbfI \
-c \
-q \
-t 140 \
-s 20



