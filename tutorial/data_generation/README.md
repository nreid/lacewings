This document explains how the tutorial data were generated. 

UConn EEB graduate student Katie Taylor generated a large set of RAD data, consisting of 141 samples from about 30 species of green lacewing (*Chrysoperla*). For the tutorial, we will use about 30 samples from 5 species that were published as part of the description of a new species. I generated the tutorial data subset as follows:

1. demultiplex original RAD data. do as little filtering as possible. write discards to a file. 
2. prepend a barcode to demultiplexed sample fastqs. samples and barcodes are from tutorial/metadata/henrysamples.txt per Henry et al 2019. these are not the original barcodes. 
3. sample sequences from discard file(s)
4. concatenate sample and discard fastqs. 
5. the final file is still large, with extremely high coverage, so subsample to 10% of original reads. 

This pool can then be analyzed in the tutorial as if it were a genuine pool of barcoded RAD sequence. 
