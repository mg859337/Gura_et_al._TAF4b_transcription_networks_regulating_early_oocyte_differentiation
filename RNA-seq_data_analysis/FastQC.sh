#!/bin/bash

# Request Resources
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job Name
#SBATCH -J FastQC_RNAseq

#Specify Output and Input Files

#SBATCH -o outputfqc.out
#SBATCH -e errorfqc.err

module load fastqc

for item in *.fastq.gz; do
	echo I see ${item}
	fastqc ${item}
done
