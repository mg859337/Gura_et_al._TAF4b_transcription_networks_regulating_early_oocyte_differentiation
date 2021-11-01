#!/bin/bash

# Request Resources
#SBATCH --time=42:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=108G

# Specify Job Name
#SBATCH -J getcounts

#Specify Output and Input Files

#SBATCH -o outputcounts.out
#SBATCH -e errorcounts.err

declare -a names=()

for thing in *.bam; do

echo "I see $thing"
names+=(${thing})

done

echo ${names[@]}

module load subread

featureCounts -p -t exon -g gene_id -a ~/data/genomes/Mouse/GRCm38/Mus_musculus.GRCm38.96.gtf -o counts.txt ${names[@]}
