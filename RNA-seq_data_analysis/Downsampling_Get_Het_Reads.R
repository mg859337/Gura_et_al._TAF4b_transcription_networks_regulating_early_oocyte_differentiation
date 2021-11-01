# A script to understand how many reads are in each E14.5 Het sample

#!/bin/bash

# Request Resources
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G

# Specify Job name
#SBATCH -J Get_reads

# Specify Job Output and Input Files
#SBATCH -o outputfilereads.out
#SBATCH -e errorfilereads.err

module load samtools

# check exact numbers of reads in each library 

for aligned in F-HET*.bam; do
	echo flagstat before downsample $aligned
	samtools flagstat $aligned > ${aligned}_predownsample.flagstat.txt
done
