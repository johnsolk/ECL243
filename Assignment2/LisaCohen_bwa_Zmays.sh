#!/bin/bash -l
#SBATCH -J ECL243_bwa
#SBATCH -p med

module load bwa/0.7.9a
module load samtools/1.2

# Zmays = domesticated maize
# Zmex = Zea mays spe. mexicana, wild teosinte from mountians of Mexico
# Tdip = Tripsacum dactyloides, "gamma grass"

bwa mem zea Zmays_R1.fq.gz Zmays_R2.fq.gz | samtools view -S -b - | samtools flagstat - > Zmays_flagstat.txt
