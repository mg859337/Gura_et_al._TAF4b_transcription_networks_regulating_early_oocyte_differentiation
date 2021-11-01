#!/bin/bash

# This script downloads scRNA-seq data from Zhao et al., 2020: GSE130212

# Request resources: 
#SBATCH --time=48:00:00 
#SBATCH --cpus-per-task=32 
#SBATCH --mem=120G

# Specify a job name: 
#SBATCH -J Dwnld_Zhao

# Specify output and error files: 
#SBATCH -o outputfiledwnld.out
#SBATCH -e errorfiledwnld.err

module load sratoolkit

array=("SRR8946219" "SRR8946221" "SRR8946224" "SRR8946222" "SRR8946223" "SRR8946220")

for element in "${array[@]}"; do
	echo I see "${element}"
    fastq-dump --gzip --split-files ${element}
done
