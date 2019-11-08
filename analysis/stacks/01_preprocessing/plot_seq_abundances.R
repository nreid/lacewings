library(tidyverse)
library(ape)

# read in 100k sequences, their gc content
f <- pipe("bioawk -c fastx '{print $seq,gc($seq)}' 107_downK.fq.gz | head -n 200000")
# f <- pipe("bioawk -c fastx '{print $seq,gc($seq)}' 007_downesi.fq.gz | head -n 200000")
# f <- pipe("bioawk -c fastx '{print $seq,gc($seq)}' ../trimmed_data/007_downesi.trim.fastq.gz | head -n 200000")
# f <- pipe("bioawk -c fastx '{print $seq,gc($seq)}' 135_downesi.fq.gz | head -n 200000")
gc <- read.table(f,stringsAsFactors=FALSE)

# bin gc content. only do this for trimmed data. 
binsize <- 0.02
gc[,2] <- (gc[,2] %/% binsize) * binsize

# get unique sequences, their counts, gc content
dgs <- group_by(gc,V1) %>% summarize(.,gc=mean(V2),len=length(V2)) %>% data.frame()

# group by gc content, make a list of sequence abundances, order it
dgss <- group_by(dgs,gc) %>% summarize(.,len=list(sort(c(len),decr=TRUE)))
dgss <- arrange(dgss,gc)

# cumulative sum of sequence abundances within gc content categories
dgsc <- dgss
dgsc[["len"]] <- lapply(dgss[["len"]],cumsum)

# plot sequence abundance distribution by gc content
unnest(dgsc,cols=len) %>% plot(.,pch=20,cex=.2,ylab="sequence count",xlab="gc content")

# get most abundance sequences:

repseq <- dgs[,3] > 200
hab <- dgs[repseq,1]
# habm <- str_split(hab,"") %>% do.call(rbind,.)
habm <- c()

fill_na <- function(x){
				l <- length(x)
				if(l < 145){x <- c(x,rep(NA,145-l))}
				return(x)
	}

for(i in 1:length(hab)){

	habm <- rbind(
		habm,
		str_split(hab[i],"") %>% 
			unlist() %>% 
			fill_na()
		)


}


# toss most abundance sequences and those very similar to them

keep <- c()
for(i in 1:nrow(dgs)){
# for(i in 1:1000){
	tmp <- dgs[i,1] %>% str_split(.,"") %>% unlist() %>% fill_na()
	
	keep <- c(
		keep, 
		apply(habm,MAR=1,FUN=function(x){mean(x!=tmp,na.rm=TRUE) < 0.05}) %>% any()
		)

}

keep <- !keep

# repeat above analysis

# group by gc content, make a list of sequence abundances, order it
dgss <- group_by(dgs[keep,],gc) %>% summarize(.,len=list(sort(c(len),decr=TRUE)))
dgss <- arrange(dgss,gc)

# cumulative sum of sequence abundances within gc content categories
dgsc <- dgss
dgsc[["len"]] <- lapply(dgss[["len"]],cumsum)

# plot sequence abundance distribution by gc content
unnest(dgsc,cols=len) %>% plot(.,pch=20,cex=.2)
