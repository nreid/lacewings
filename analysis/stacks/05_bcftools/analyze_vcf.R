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

# missing data filter threshold
nm <- rowSums(is.na(gts)) < 55

# total site depth:
str_extract(vcf2[,8],regex("(?<=DP=)[0-9]+")) %>% as.numeric() %>% log(.,10) %>% hist(.,breaks=500)
dep <- str_extract(vcf2[,8],regex("(?<=DP=)[0-9]+")) %>% as.numeric()
# depth threshold
dp <- dep > 50000 & dep < 400000

gts <- gts[nm & dp,]
vcf2 <- vcf2[nm & dp,]

str_extract(vcf2[,8],regex("(?<=DP=)[0-9]+")) %>% as.numeric() %>% log(.,10) %>% hist(.,breaks=500)


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
# rescale <- nrow(gts)/600000
# dmat <- dmat * rescale


spec <- gsub("^....","",colnames(gts))
uspec <- unique(spec)

# samples to retain (vector are phylogenetically misplaced samples, possible mis-IDs)
keep <- colMeans(!is.na(gts)) > 0.2
keep <- keep & !(colnames(gts) %in% c("079_heidarii","120_rufilabris","129_comanche","171_adamsi","096_carnea","134_plorabunda","042_lucasina","010_downesi","007_downesi"))


# tree and MDS plot
tr <- as.dist(dmat[keep,keep]) %>% nj() %>% midpoint.root() 
plot(tr,cex=.5,type="phylogram",tip.color=as.numeric(as.factor(gsub("^....","",tr$tip.label))))
plot(tr,cex=.5,type="phylogram")
add.scale.bar()
plot(tr,show.tip.label=FALSE)
tiplabels(gsub("^....","",tr$tip.label),col=as.numeric(as.factor(gsub("^....","",tr$tip.label))),cex=0.5,frame="none")

mds <- as.dist(dmat[keep,keep]) %>% cmdscale(.,k=6)
i <- 1 ; j <- 2
plot(mds[,i],mds[,j],col=as.factor(gsub("^....","",rownames(mds))))
# text(mds[,i],mds[,j],labels=rownames(mds),col=as.numeric(as.factor(gsub("^....","",rownames(mds)))),cex=0.5)
text(mds[,i],mds[,j],labels=rownames(mds),cex=0.5)


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

hist(btwn,breaks=seq(0,0.15,0.001),freq=FALSE,col=rgb(1,0,0,.5))
within %>% hist(.,freq=FALSE,add=TRUE,breaks=seq(0,0.15,0.001),col=rgb(0,1,0,.5))
diag(dmat[keep,keep]) %>% hist(.,freq=FALSE,add=TRUE,breaks=seq(0,0.15,0.001),col=rgb(0,0,1,0.5))
	

# allele frequencies by species

freqs <- c()
counts <- c()
for(i in 1:length(uspec)){
	if(sum(spec[keep]==uspec[i]) > 1){
	freqs <- cbind(freqs, rowMeans(gts[,keep][,spec[keep]==uspec[i]],na.rm=TRUE))
	counts <- cbind(counts,rowSums(!is.na(gts[,keep][,spec[keep]==uspec[i]]))*2)
	}
	if(sum(spec[keep]==uspec[i]) == 1){
	freqs <- cbind(freqs, gts[,keep][,spec[keep]==uspec[i]])
	counts <- cbind(counts,2*(!is.na(gts[,keep][,spec[keep]==uspec[i]])))
	}

}

colnames(freqs) <- uspec
colnames(counts) <- uspec

pmat <- matrix(nrow=length(uspec),ncol=length(uspec))
rownames(pmat) <- uspec
colnames(pmat) <- uspec

for(i in 1:length(uspec)){

	for(j in 1:length(uspec)){
		if(i==j){next()}
		pmat[i,j] <- mean((freqs[,i] * (1 - freqs[,j])) + (freqs[,j] * (1 - freqs[,i])),na.rm=TRUE)

	}
}

pmat2 <- matrix(nrow=length(uspec),ncol=length(uspec))
rownames(pmat2) <- uspec
colnames(pmat2) <- uspec

for(i in 1:length(uspec)){

	for(j in 1:length(uspec)){
		if(i==j){next()}
		pmat2[i,j] <- mean(dmat[keep,keep][spec[keep]==uspec[i],spec[keep]==uspec[j]])

	}
}

par(mfrow=c(1,2))
ptr <- pmat %>% as.dist() %>% nj() %>% midpoint.root()
plot(ptr,"unrooted")
ptr2 <- pmat2 %>% as.dist() %>% nj() %>% midpoint.root()
plot(ptr2,"unrooted")


levelplot(pmat[hclust(as.dist(pmat))$order,hclust(as.dist(pmat))$order])
levelplot(pmat2[hclust(as.dist(pmat2))$order,hclust(as.dist(pmat2))$order])
levelplot(dmat[keep,keep][hclust(as.dist(dmat[keep,keep]))$order,hclust(as.dist(dmat[keep,keep]))$order])

# look at heterozygote deficits

# genotype frequencies
species <- "calocedrii"
nsam <- sum(spec==species)
gtf <- (apply(gts[counts[,species]==(nsam*2),keep & spec==species],MAR=1,FUN=function(x){factor(x,levels=c(0,0.5,1)) %>% table()}) %>% t())

sim <- apply(
	cbind(as.vector(table(freqs[counts[,species]==(nsam*2),species])),alt=(0:(nsam*2))/(nsam*2)),
	MAR=1,
	FUN=function(x){
		rmultinom(n=x[1],prob=c((1-x[2])^2,2*(1-x[2])*x[2],x[2]^2),size=nsam)
		}) %>% 
	do.call(cbind,.) %>% t()

plot(jitter(freqs[counts[,species]==(nsam*2),species]^2,factor=5),jitter(gtf[,3]/nsam),pch=20,col=rgb(0,0,0,.3))
points(jitter(((sim[,2]+2*sim[,3])/(nsam*2))^2,factor=5),jitter(sim[,3])/nsam,pch=20,col=rgb(1,0,0,.3))

par(mfrow=c(1,2))
plot(jitter(freqs[counts[,species]==(nsam*2),species]^2,factor=5),jitter(gtf[,3]/nsam),pch=20,col=rgb(0,0,0,.3))
abline(0,1)
plot(jitter(((sim[,2]+2*sim[,3])/(nsam*2))^2,factor=5),jitter(sim[,3])/nsam,pch=20,col=rgb(1,0,0,.3))
abline(0,1)

# look at allele depths
# pick a species
species <- "calocedrii"

# sites that vary within that species
sites <- which((rowMeans(gts[,spec==species],na.rm=TRUE) > 0) & (rowMeans(gts[,spec==species],na.rm=TRUE) < 1))

gtmat <- head(vcf2[sites,-c(1:9)],n=20000)[,spec==species] %>% unlist() %>% str_split(.,":|,") %>% do.call(rbind,.) %>% data.frame(.,stringsAsFactors=FALSE)
for(i in 2:7){gtmat[,i] <- as.numeric(gtmat[,i])}

sg <- as.numeric(gtmat[,5]) > 15 #& as.numeric(gtmat[,6]) > 0 & as.numeric(gtmat[,7]) > 0

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

# calculate Fst

fst <- function(sp1,sp2,fr=freqs,co=counts){

	pi1 <- mean((2 * (freqs[,sp1] * (1-freqs[,sp1])) * counts[,sp1] / (counts[,sp1]-1)),na.rm=TRUE)
	pi2 <- mean((2 * (freqs[,sp2] * (1-freqs[,sp2])) * counts[,sp2] / (counts[,sp2]-1)),na.rm=TRUE)
	pibtwn <- mean((freqs[,sp1] * (1 - freqs[,sp2])) + ((1-freqs[,sp1]) * freqs[,sp2]),na.rm=TRUE)

	ff <- 1 - (pi1 + pi2) / 2 / pibtwn

	data.frame(stats=c(fst=ff,pi1=pi1,pi2=pi2,pibtwn=pibtwn))

}


fst("johnsoni","adamsi")

