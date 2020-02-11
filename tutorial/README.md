### Introduction

Restriction-associated DNA sequencing (RAD-seq) approaches are a collection of methods used to assay genetic variation in groups of individuals. These approaches are a form of **reduced representation sequencing**. This means that some relatively small fraction of sites in the genome, usually widely dispersed, are targeted for sequencing. This is accomplished by digesting the genome using 1 or more restriction enzymes, and sequencing DNA adjacent to the cut site(s). Many types of questions (see discussion in Lowry et al. 2016 and Catchen et al. 2017) do not necessarily require the kind of genomically dense data produced by whole genome sequencing (WGS), so RAD-seq provides two benefits: 

	1. It is cost effective, allowing a much larger number of individuals to be genotyped than an equivalent WGS study. 
	2. It is computationally effective. A data file containing genotypes for several thousand sites is much more portable and easy to analyze than one potentially containing tens of millions of sites. 

This tutorial will cover two popular software packages for analyzing RAD data. For each it will demonstrate two approaches for generating genotypes: _de novo_ assembly and reference mapping. 


### Sources of data:

In this tutorial we will use data from two sources. 

The first source of data is **Henry et al. 2019** (HTJ):

    Henry, Charles S., Katherine L. Taylor, and James B. Johnson. "A new lacewing species of the Chrysoperla carnea species-group from central Asia associated with conifers (Neuroptera: Chrysopidae)." Journal of Natural History 53.21-22 (2019): 1277-1300.

These data are a subset of data collected for a larger project studying the phylogeny of green lacewings (_Chrysoperla_). The paper used these data to support the description of a new lacewing species in the _Chrysoperla carnea_ species group from central Asia, _Chrysoperla duelli_. Lacewings in this group are morphologically identical, but have distinctive courtship songs. _C. duelli's_ song, however, is highly similar to a song by another species in the group that lives in North America, _C. downesi_. 

The original data consist of 3 RAD libraries, each containing 47 barcoded samples. The authors used SbfI to digest genomic DNA, then randomly sheared and sequenced a single end on the Illumina platform. 

33 of these samples are used in HTJ. Here I synthesize a multiplexed pool of the 33 HTJ samples from the original larger pools to demonstrate a RAD analysis from start to finish (see directory "data_generation"). 

Scripts for exploratory analysis of the entire dataset are located in "../analysis"

The second source of data is an unpublished dataset by Wegrzyn and colleagues (WEA). It's a genetic association study of resistance of the green ash _Fraxinus pennsylvanica_ to an invasive beetle, the emerald ash borer (_Agrilus planipennis_). 

The total dataset consists of 6 libraries of 10-16 barcoded samples. WEA use ddRAD (Peterson et al. 2012). In ddRAD, 2 enzymes are used to digest genomic DNA and only fragments cut by both enzymes are sequenced. WEA additionally use 


### References

Catchen, Julian M., Paul A. Hohenlohe, Louis Bernatchez, W. Chris Funk, Kimberly R. Andrews, and Fred W. Allendorf. 2017. “Unbroken: RADseq Remains a Powerful Tool for Understanding the Genetics of Adaptation in Natural Populations.” Molecular Ecology Resources.

Lowry, David B., Sean Hoban, Joanna L. Kelley, Katie E. Lotterhos, Laura K. Reed, Michael F. Antolin, and Andrew Storfer. 2016. “Breaking RAD: An Evaluation of the Utility of Restriction Site-Associated DNA Sequencing for Genome Scans of Adaptation.” Molecular Ecology Resources, November. https://doi.org/10.1111/1755-0998.12635.

Peterson, Brant K., Jesse N. Weber, Emily H. Kay, Heidi S. Fisher, and Hopi E. Hoekstra. 2012. “Double Digest RADseq: An Inexpensive Method for de Novo SNP Discovery and Genotyping in Model and Non-Model Species.” PloS One 7 (5): e37135.



