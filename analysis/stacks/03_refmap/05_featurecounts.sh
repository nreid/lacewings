#!/bin/bash
#SBATCH --job-name=featurecounts
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


INDIR=../results/trimmed_aligned
find $INDIR -name "*bam" >$INDIR/bams.list

OUTDIR=../results/trimmed_aligned_stats

META=../../../metadata

# run featureCounts to get a sense of the number of fragments on RAD sites
	# do left and right reads separately
# first generate a SAF file for featureCounts:
	# use dummy gene IDs

cat $META/sbf1_left.bed $META/sbf1_right.bed | sort -V >$META/sbf1_flanks.bed

echo "GeneID	Chr	Start	End	Strand" >$META/sbf1_flanks.saf
cat $META/sbf1_flanks.bed | \
awk '{OFS="\t"}{s=$2+1}{print NR,$1,s,$3,"+"}' >>$META/sbf1_flanks.saf

# run featurecounts
featureCounts \
-a $META/sbf1_flanks.saf \
-o $OUTDIR/sbf1_flanks_counts.txt \
-Q 30 \
-F SAF \
--primary \
-p \
-T 12 \
$(cat $INDIR/bams.list | tr "\n" " ")


# run featurecounts on potential rad sites identified from coverage:

# identify and summarize intervals
zcat ../results/trimmed_aligned_stats/depthperbase.txt.gz | \
awk '$3 > 10' | \
awk '{OFS="\t"}{s=$2-1}{print $1,s,$2,$3}' | \
bedtools merge -i stdin -c 4 -o mean,median,min,max,count | \
awk '$5 > 50' | \
awk '$8 > 50' \
>$META/possibleradsites.bed

# generate SAF file
echo "GeneID	Chr	Start	End	Strand" >$META/possibleradsites.saf
cat $META/possibleradsites.bed | \
awk '{OFS="\t"}{s=$2+1}{print NR,$1,s,$3,"+"}' >>$META/possibleradsites.saf

# run featurecounts
featureCounts \
-a $META/possibleradsites.saf \
-o $OUTDIR/possibleradsites_counts.txt \
-Q 30 \
-F SAF \
--primary \
-p \
-T 12 \
$(cat $INDIR/bams.list | tr "\n" " ")


# also get the intersection of known rad sites and the potential rad sites (with their coverage)

bedtools intersect -wa -u \
-a <(zcat ../results/trimmed_aligned_stats/depthperbase.txt.gz | \
awk '$3 > 10' | \
awk '{OFS="\t"}{s=$2-1}{print $1,s,$2,$3}' | \
bedtools merge -i stdin -c 4 -o mean,median,min,max,count | \
awk '$5 > 50' | \
awk '$8 > 50') \
-b $META/sbf1_flanks.bed \
>$META/radsite_intersect.bed


