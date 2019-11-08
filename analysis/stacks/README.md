# Analysis with STACKS

## Preprocessing steps
1. 01_preprocessing/01_process_radtags.sh
	- demultiplex sequences from three pools of 47 individuals. 
	- excludes reads based on quality and presence of adapter
2. 01_preprocessing/02_fastqc.sh
	- run fastqc on all demultiplexed samples
3. 01_preprocessing/03_multiqc.sh
	- aggregate fastqc reports
4. 01_preprocessing/04_quality_trim.sh
	- run trimmomatic to quality trim sequences
	- adapter still present in some sequences after process_radtags

Results in two sets of fastq files, one straight from process_radtags, one that has also been trimmed. 

The first set will be de novo. The second will be reference aligned. 

## Analysis

### stacks de novo

01. 02_de_novo

```
for file in *snps.tsv.gz;  do zcat $file | awk -v var=$file '{OFS="\t"} { if($7 !~ /-/){n+=1}} END {print var,n,NR,n/NR}'| cut -c 5- ; done | column | sort
```