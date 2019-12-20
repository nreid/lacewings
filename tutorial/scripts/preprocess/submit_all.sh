#!/bin/bash

# submit all preprocessing steps as dependencies. 

jid1=$(sbatch 01_processradtags.sh)
jid2=$(sbatch --dependency=afterok:$jid1 02_fastqc.sh)
jid3=$(sbatch --dependency=afterok:$jid2 03_multiqc.sh)
jid4=$(sbatch --dependency=afterok:$jid3 04_qualitytrim.sh)
jid5=$(sbatch --dependency=afterok:$jid4 05_trimfastqc.sh)
jid6=$(sbatch --dependency=afterok:$jid5 06_trimmultiqc.sh)
