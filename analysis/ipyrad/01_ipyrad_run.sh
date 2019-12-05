#!/bin/bash 
#SBATCH --job-name=ipyrad_full
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load ipyrad/0.9.19


# create population file:

cat ../../metadata/popmap.txt >../../metadata/popmap.ipyrad.txt
POP=../../metadata/popmap.ipyrad.txt
cut -f 2 ../../metadata/popmap.txt | sort | uniq | sed 's/$/:0 /' | tr "\n" " " | sed 's/^/# /' >>$POP

ipyrad -p params-fulldata.txt -s 1 -c 32

echo step 1 done

ipyrad -p params-fulldata.txt -s 2 -c 32

echo step 2 done

ipyrad -p params-fulldata.txt -s 3 -c 32

echo step 3 done

ipyrad -p params-fulldata.txt -s 4 -c 32

echo step 4 done

ipyrad -p params-fulldata.txt -s 5 -c 32

echo step 5 done

ipyrad -p params-fulldata.txt -s 6 -c 32

echo step 6 done

ipyrad -p params-fulldata.txt -s 7 -c 32

echo step 7 done

