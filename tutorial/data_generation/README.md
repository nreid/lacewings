This document explains how the tutorial data were generated. 

1. demultiplex original RAD data. do as little filtering as possible. write discards to a file. 
2. prepend a barcode to demultiplexed sample fastqs. samples and barcodes are from tutorial/metadata/henrysamples.txt per Henry et al 2019. these are not the original barcodes. 
3. sample sequences from discard file(s)
4. concatenate sample and discard fastqs. 
5. shuffle pooled fastq. 

This pool can then be analyzed in the tutorial as if it were a genuine pool of barcoded RAD sequence. 