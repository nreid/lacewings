#!/bin/bash

jid1=$(sbatch --parsable 01_denovo-ustacks.sh)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_denovo-cstacks.sh)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 03_denovo-sstacks.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 04_denovo-tsv2bam.sh)
jid5=$(sbatch --parsable --dependency=afterok:$jid4 05_denovo-gstacks.sh)
jid6=$(sbatch --parsable --dependency=afterok:$jid5 06_denovo-populations.sh)
jid7=$(sbatch --parsable --dependency=afterok:$jid6 R01_refmappl.sh)
