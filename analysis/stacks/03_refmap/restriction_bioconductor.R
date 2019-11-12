library(tidyverse)
library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()


# This R script finds sbf1 cut sites, flanking sites and writes them to bed files

# read in reference genome 
seqs <- readDNAStringSet("../results/redundans_metaquast_filtered.nomt.masked.fasta")

# find sbf1 sites
mp <- vmatchPattern("CCTGCAGG",seqs)

# make it a GRanges object
mp2 <- as(mp, "GRanges")

# make it a data frame
mp3 <- as.data.frame(mp2)

# 0-index start
mp3[,2] <- mp3[,2] - 1

# left flank
left <- mp3
	left[,2] <- left[,2] - 9
	left[,3] <- left[,3] - 9
left[left[,2] < 0,2] <- 0

# right flank
right <- mp3
	right[,2] <- right[,2] + 9
	right[,3] <- right[,3] + 9

# find sbf1 1-off sites
mp_1 <- vmatchPattern("CCTGCAGG",seqs,max.mismatch=1)

# make it a GRanges object
mp2_1 <- as(mp_1, "GRanges")

# make it a data frame
mp3_1 <- as.data.frame(mp2_1)

# 0-index start
mp3_1[,2] <- mp3_1[,2] - 1


# write tables
write.table(mp3[,1:3],"../../../metadata/sbf1.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(mp3_1[,1:3],"../../../metadata/sbf1_off.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(left[,1:3],"../../../metadata/sbf1_left.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(right[,1:3],"../../../metadata/sbf1_right.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

# There are 24 possible off-by-1 cut sites. 

cuts <- getSeq(seqs,mp2_1)

as.character(cuts) %>% table() %>% sort()

