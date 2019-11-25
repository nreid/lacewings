library(tidyverse)
library(GenomicRanges)
library(beeswarm)


# This script examines coverage among individuals for intervals
# identified as having substantial sequencing coverage
# It writes out a set of target regions for variant calling. 

# genomic intervals with coverage
pr <- read.table("../results/trimmed_aligned_stats/possibleradsites_counts.txt",stringsAsFactors=FALSE,header=TRUE)
colnames(pr)[-c(1:6)] <- str_extract(colnames(pr)[-c(1:6)],regex("[0-9]+[^.]+"))

# sbf1 cut sites
sbf1 <- read.table("../../../metadata/sbf1.bed",stringsAsFactors=FALSE)

# sbf1 + 1-off cut sites
sbf1off <- read.table("../../../metadata/sbf1_off.bed",stringsAsFactors=FALSE)

# SN stats from samtools

sn <- read.table("../results/trimmed_aligned_stats/SN.txt",stringsAsFactors=FALSE,header=TRUE,sep="\t",quote="")
sn <- t(sn)
rownames(sn) <- rownames(sn) %>% gsub(".([0-9][^\\.]*).trim","\\1",.)

plot(sn[,1],sn[,7]/sn[,1])
text(sn[,1],sn[,7]/sn[,1],labels=rownames(sn))

# generate Granges objects
prg <- GRanges(
	seqnames = pr[,2],
	ranges = IRanges(start = pr[,3], end = pr[,4])
	)

sbf1g <- GRanges(
	seqnames = sbf1[,1],
	ranges = IRanges(start = sbf1[,2]+1, end = sbf1[,3])
	)

sbf1offg <- GRanges(
	seqnames = sbf1off[,1],
	ranges = IRanges(start = sbf1off[,2]+1, end = sbf1off[,3])
	)


# get overlaps between coverage intervals and cut sites
# warnings are fine. 
olaps <- findOverlaps(prg,sbf1g) %>% as.matrix()
olapsoff <- findOverlaps(prg,sbf1offg) %>% as.matrix()

# species names
spec <- gsub("^....","",colnames(pr)[-c(1:6)])

# vector of intervals with cut sites
rr <- rep(FALSE,nrow(pr))
rr[olaps[,1]] <- TRUE

# and 1-off cut sites. 
rroff <- rep(FALSE,nrow(pr))
rroff[olapsoff[,1]] <- TRUE

# missing data and read counts across individuals
md <- rowSums(pr[,-c(1:6)] > 0)
co <- rowSums(pr[,-c(1:6)])

# subset the sites
keep <- md > 20 & co < 8e5 & co > 3000
pr2 <- pr[keep,]
rr2 <- rr[keep]
rroff2 <- rroff[keep]
md2 <- md[keep]

# scaled counts
pr2s <- apply(pr2[,-c(1:6)],MAR=2,FUN=function(x){x/sum(x)*1e6})

# vector to exclude low cov individuals
colMeans(pr2s[md2 > 100,]>0) %>% plot()
ind <- colMeans(pr2s[md2 > 100,]>0) > 0.5

# plot missing data as a function of total coverage
plot(log(rowSums(pr2s),10),rowSums(pr2s>0),col=rr2+rroff2+1)

# plot missing data as a function of upper -tile coverage
tile <- apply(pr2s[,ind],MAR=1,FUN=quantile,prob=0.88)
plot(log(tile,10),rowSums(pr2s>0),col=rr2+rroff2+1,pch=20,cex=0.4)

# beeswarm as a function of no cut site (1), 1-off (2) and cut site (3)
cutstatus <- rr2+rroff2+1
beeswarm(rowSums(pr2s[,ind]>0) ~ cutstatus,pch=20,cex=.5)
beeswarm(log(tile,10) ~ cutstatus,pch=20,cex=.5)

# write out target regions

targets <- cbind(scaf=pr2[,2],start=pr2[,3]-1,end=pr2[,4])
write.table(targets,"../../../metadata/targets.bed",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)




