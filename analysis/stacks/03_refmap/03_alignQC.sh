#!/bin/bash
#SBATCH --job-name=alignQC
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general


module load subread/1.6.0
module load samtools/1.7
module load R/3.6.0
module load bedtools/2.27.1

# This script does some QC on the alignments. 

INDIR=../results/trimmed_aligned

# make a list of bam files
find $INDIR -name "*bam" >$INDIR/bams.list

# # make the bam indexes
find $INDIR -name "*bam" | xargs -P 12 -I {} samtools index {}

# make an output directory
OUTDIR=../results/trimmed_aligned_stats
mkdir -p $OUTDIR

# samtools bam statistics
for file in $(find $INDIR -name "*bam"); 
do samtools stats $file >${file}.stats
echo $file;
done

mv $INDIR/*stats $OUTDIR

# put the basic stats all in one file. 
grep ^SN $OUTDIR/Pool1_BC01.bam.stats | cut -f 2 > $OUTDIR/SN.txt
for file in $(find $OUTDIR -name "Pool*stats" | sort)
do paste $OUTDIR/SN.txt <(grep ^SN $file | cut -f 3) > $OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
done

# add a header
find $OUTDIR -name "Pool*stats" | sort | sed 's/.bam.*//' | sed 's/.*\///' | tr "\n" "\t" | sed 's/\t$/\n/'>$OUTDIR/SN2.txt
cat \n >>$OUTDIR/SN2.txt
cat $OUTDIR/SN.txt >>.$OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt

# run featureCounts to get a sense of the number of fragments on RAD sites
	# do left and right reads separately
# first generate a SAF file for featureCounts:
	# use dummy gene IDs

cat ../../../metadata/sbf1_left.bed ../../../metadata/sbf1_right.bed | sort -V >sbf1_flanks.bed

echo "GeneID	Chr	Start	End	Strand" >../../../metadata/sbf1_flanks.saf
cat ../../../metadata/sbf1_flanks.bed | \
awk '{OFS="\t"}{s=$2+1}{print NR,$1,s,$3,"+"}' >>../../../metadata/sbf1_flanks.saf

# run featurecounts
featureCounts \
-a ../../../metadata/sbf1_flanks.saf \
-o $OUTDIR/sbf1_flanks_counts.txt \
-Q 30 \
-F SAF \
--primary \
-p \
-T 12 \
$(cat ../results/aligned_ref/bams.list | tr "\n" " ")





