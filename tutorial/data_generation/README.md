This document explains how the tutorial data were generated. 

1. demultiplex original RAD data. do as little filtering as possible. write discards to a file. 
2. prepend a barcode to demultiplexed sample fastqs. samples and barcodes are from tutorial/metadata/henrysamples.txt per Henry et al 2019. these are not the original barcodes. 
3. concatenate sample fastqs. 
4. sample sequences from discard file(s) and concatenate with pooled sample fastq. 
5. shuffle pooled fastq. 