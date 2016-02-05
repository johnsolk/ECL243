#!/bin/bash -l
#SBATCH -J ECL243_readcounts
#SBATCH -p med

echo "This is the number of sequences that start with 'GGTCA':"
zcat /home/ljcohen/ECL243/bkn009_R1.fq.gz | grep '^GGTCA' | wc -l
