#!/bin/bash

# A script that downsamples the Hets based on the randomly generated numbers

# Request Resources
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G

# Specify Job name
#SBATCH -J Down

# Specify Job Output and Input Files
#SBATCH -o outputfiledown.out
#SBATCH -e errorfiledown.err

# downsample files and get flagstat

module load samtools

samtools view -s 0.29 -b F-HET-1.bam > F-HET-1_downsampled.bam
samtools flagstat F-HET-1_downsampled.bam > F-HET-1_downsampled.flagstat.txt

samtools view -s 0.60 -b F-HET-2.bam > F-HET-2_downsampled.bam
samtools flagstat F-HET-2_downsampled.bam > F-HET-2_downsampled.flagstat.txt

samtools view -s 0.59 -b F-HET-3.bam > F-HET-3_downsampled.bam
samtools flagstat F-HET-3_downsampled.bam > F-HET-3_downsampled.flagstat.txt

samtools view -s 0.54 -b F-HET-4.bam > F-HET-4_downsampled.bam
samtools flagstat F-HET-4_downsampled.bam > F-HET-4_downsampled.flagstat.txt

