library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()


# This R script finds sbf1 cut sites, flanking sites and writes them to bed files

seqs <- readDNAStringSet("../results/redundans_metaquast_filtered.nomt.masked.fasta")

mp <- vmatchPattern("CCTGCAGG",seqs)

mp2 <- as(mp, "GRanges")

mp3 <- as.data.frame(mp2)

mp3[,2] <- mp3[,2] - 1

left <- mp3
	left[,2] <- left[,2] - 9
	left[,3] <- left[,3] - 9
left[left[,2] < 0,2] <- 0

right <- mp3
	right[,2] <- right[,2] + 9
	right[,3] <- right[,3] + 9

write.table(mp3[,1:3],"../../../metadata/sbf1.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(left[,1:3],"../../../metadata/sbf1_left.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(right[,1:3],"../../../metadata/sbf1_right.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)