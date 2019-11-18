library(tidyverse)
library(ape)
library(phytools)

f <- pipe('gzcat ../results/bcftools/chrysoperla.vcf.gz | head -n 50000 | grep "CHROM"')
h <- scan(f,what="character")
vcf <- read.table("../results/bcftools/chrysoperla.vcf.gz",stringsAsFactors=FALSE) 
colnames(vcf) <- h

bial <- nchar(vcf[,5])==1
hq <- vcf[,6] > 30

vcf2 <- vcf[bial & hq,]

hz <- apply(vcf2[,-(1:9)],MAR=2,FUN=function(x){c(sum(grepl("0/1",x)),sum(!grepl("\\./\\.",x)))}) %>% t()

hz <- cbind(hz,hz=hz[,1]/hz[,2])

gts <- apply(vcf2[,-c(1:9)],MAR=2,FUN=function(x){x <- gsub(":.*","",x); tr <- c("0/0"=0,"0/1"=0.5,"1/1"=1,"./."=NA); x <- tr[x]})

dfunc <- function(p1,p2){

	mean(p1 * (1-p2) + (1-p1) * p2,na.rm=TRUE)

}

dmat <- matrix(nrow=138,ncol=138)
rownames(dmat) <- h[-c(1:9)]
colnames(dmat) <- h[-c(1:9)]
for(i in 1:138){

	for(j in 1:138){

		dmat[i,j] <- dfunc(gts[,i],gts[,j])

	}

	print(i)
}

# roughly rescale divergence by sites genotyped
rescale <- nrow(gts)/600000
dmat <- dmat * rescale


spec <- gsub("^....","",colnames(gts))
uspec <- unique(spec)

keep <- colMeans(!is.na(gts)) > 0.2

within <- c()
for(i in 1:length(uspec)){
	tmp <- spec[keep]==uspec[i]
	within <- c(within,as.vector(dmat[keep,keep][tmp,tmp]))

}

btwn <- c()
for(i in 1:(length(uspec)-1)){

	tmp1 <- spec[keep]==uspec[i]
	for(j in (i+1):length(uspec)){
		tmp2 <- spec[keep]==uspec[j]
		print(length(as.vector(dmat[keep,keep][tmp1,tmp2])))
		btwn <- c(btwn, as.vector(dmat[keep,keep][tmp1,tmp2]))
	}
}

hist(btwn,breaks=seq(0,0.03,0.001),freq=FALSE,col=rgb(1,0,0,.5))
within %>% hist(.,freq=FALSE,add=TRUE,breaks=seq(0,0.03,0.001),col=rgb(0,1,0,.5))
diag(dmat[keep,keep]) %>% hist(.,freq=FALSE,add=TRUE,breaks=seq(0,0.03,0.001),col=rgb(0,0,1,0.5))
	

tr <- as.dist(dmat[keep,keep]) %>% nj() %>% midpoint.root() 
plot(tr,cex=.5,type="phylogram")
add.scale.bar()
plot(tr,show.tip.label=FALSE)
tiplabels(gsub("^....","",tr$tip.label),col=as.numeric(as.factor(spec[keep])),cex=0.5,frame="none")

# allele frequencies by species
freqs <- (group_by(data.frame(spec=spec[keep],t(gts[,keep])),spec) %>% summarize_all(funs(mean(., na.rm = TRUE))))
counts <- (group_by(data.frame(spec=spec[keep],t(gts[,keep])),spec) %>% summarize_all(funs(sum(!is.na(.), na.rm = TRUE))))

fd <- data.frame(freqs[-1]) %>% t()
colnames(fd) <- freqs[[1]]
cd <- data.frame(counts[-1]) %>% t()
colnames(cd) <- counts[[1]]


rowSums( fd > 0 & fd < 1,na.rm=TRUE) %>% table()

data.frame(spec=colnames(fd),var=colMeans( fd > 0 & fd < 1,na.rm=TRUE),na=colMeans(is.na(fd))) %>% arrange(.,var)

sharedmat <- matrix(nrow=length(uspec),ncol=length(uspec))
colnames(sharedmat) <- colnames(fd)
rownames(sharedmat) <- colnames(fd)

for(i in 1:length(uspec)){

	for(j in 1:length(uspec)){

		sharedmat[i,j] <- mean(fd[,i] > 0 & fd[,i] < 1 & fd[,j] > 0 & fd[,i] < 1,na.rm=TRUE)

	}

}

ftr <- dist(sharedmat) %>% nj()


# look at heterozygote deficits

# genotype frequencies
species <- "calocedrii"
nsam <- sum(spec==species)
gtf <- (apply(gts[cd[,species]==nsam,spec==species],MAR=1,FUN=function(x){factor(x,levels=c(0,0.5,1)) %>% table()}) %>% t())

sim <- apply(
	cbind(as.vector(table(fd[cd[,species]==7,species])),alt=(0:(nsam*2))/(nsam*2)),
	MAR=1,
	FUN=function(x){
		rmultinom(n=x[1],prob=c((1-x[2])^2,2*(1-x[2])*x[2],x[2]^2),size=nsam)
		}) %>% 
	do.call(cbind,.) %>% t()

plot(jitter(fd[cd[,species]==7,species]^2,factor=5),jitter(gtf[,3]/nsam),pch=20,col=rgb(0,0,0,.3))
points(jitter(((sim[,2]+2*sim[,3])/(nsam*2))^2,factor=5),jitter(sim[,3])/nsam,pch=20,col=rgb(1,0,0,.3))



# look at allele depths
gtmat <- head(vcf2[,-c(1:9)],n=5000) %>% unlist() %>% str_split(.,":|,") %>% do.call(rbind,.) %>% data.frame(.,stringsAsFactors=FALSE)

sg <- as.numeric(gtmat[,5]) > 15

plot(as.numeric(gtmat[sg,6]),as.numeric(gtmat[sg,7]),col=as.factor(gtmat[sg,1]),pch=20,cex=.2)

hist(as.numeric(gtmat[sg,5]),breaks=500)

# simulate allele depths

gt_dep_counts <- table(gtmat[sg,1],as.numeric(gtmat[sg,5])) %>% as.matrix() %>% t()

sim_ad <- c()

for(i in 1:nrow(gt_dep_counts)){
	tmp1 <- rbinom(gt_dep_counts[i,1],prob=0.99,size=15+i)
	tmp2 <- rbinom(gt_dep_counts[i,2],prob=0.50,size=15+i)
	tmp3 <- rbinom(gt_dep_counts[i,3],prob=0.01,size=15+i)

	tmp <- rbind(
			data.frame(gen=rep("0/0",length(tmp1)),ref=tmp1,alt=(15+i-tmp1)),
			data.frame(gen=rep("0/1",length(tmp2)),ref=tmp2,alt=(15+i-tmp2)),
			data.frame(gen=rep("1/1",length(tmp3)),ref=tmp3,alt=(15+i-tmp3))
		)



	sim_ad <- rbind(
		sim_ad,
		tmp
		)
}

# as expected, the empirical allele depth plot is nuts. 
plot(sim_ad[,2],sim_ad[,3],col=as.factor(sim_ad[,1]),pch=20,cex=.2)
