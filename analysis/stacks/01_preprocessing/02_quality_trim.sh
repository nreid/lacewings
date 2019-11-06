#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-140]%20

module load java
module load Trimmomatic/0.36

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

#input/output directories, supplementary files

ADAPT=../../../metadata/TruSeq3-SE.fa

INDIR=../results/demultiplexed_fastqs
FASTQS=($(ls -1 $INDIR/*fastq.gz))
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
OUTFILE=$(echo $INFILE | sed 's/fastq/trim.fastq/')

# make trimmed directory if it doesn't exist
mkdir -p ../results/trimmed_data
OUTDIR=../results/trimmed_data

java -jar $Trimmomatic SE \
-threads 10 \
$INDIR/$INFILE \
$OUTDIR/$OUTFILE \
ILLUMINACLIP:$ADAPT:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:45

