#!/bin/bash -l
#SBATCH -J ECL243_challenge
#SBATCH -p med

for i in {A,T,G,C}{A,T,G,C}{A,T,G,C}; 
do 
echo ${i}":"
zcat /home/ljcohen/ECL243/bkn009_R1.fq.gz | grep ^${i} | wc -l; 
done
