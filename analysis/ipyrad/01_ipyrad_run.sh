#!/bin/bash 
#SBATCH --job-name=ipyrad_full
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load ipyrad/0.7.22


# create population file:

cat ../../metadata/popmap.txt >../../metadata/popmap.ipyrad.txt
POP=../../metadata/popmap.ipyrad.txt
cut -f 2 ../../metadata/popmap.txt | sort | uniq | sed 's/$/:0 /' | tr "\n" " " | sed 's/^/# /' >>$POP

ipyrad -p params-fulldata.txt -s 1

echo step 1 done

ipyrad -p params-fulldata.txt -s 2

echo step 2 done

ipyrad -p params-fulldata.txt -s 3

echo step 3 done

ipyrad -p params-fulldata.txt -s 4

echo step 4 done

ipyrad -p params-fulldata.txt -s 5

echo step 5 done

ipyrad -p params-fulldata.txt -s 6

echo step 6 done

ipyrad -p params-fulldata.txt -s 7

echo step 7 done

