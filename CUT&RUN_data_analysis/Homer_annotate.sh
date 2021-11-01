#!/bin/bash

# Request Resources
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=246G

# Specify Job Name
#SBATCH -J Homer_CUTRUN_Annotate

#Specify Output and Input Files
#SBATCH -o outputhomer_annotate.out
#SBATCH -e errorhomer_annotate.err

module load homer
module load samtools

# Get general annotation of peaks
annotatePeaks.pl F-TAF4b_1_peaks.txt ~/data/genomes/Mouse/GRCm39/genome.fa -gtf ~/data/genomes/Mouse/GRCm39/Mus_musculus.GRCm39.103.gtf > F-TAF4b_1_annotated_peaks.txt
annotatePeaks.pl F-TAF4b_2_peaks.txt ~/data/genomes/Mouse/GRCm39/genome.fa -gtf ~/data/genomes/Mouse/GRCm39/Mus_musculus.GRCm39.103.gtf > F-TAF4b_2_annotated_peaks.txt

annotatePeaks.pl F-H3K4_1_regions.txt ~/data/genomes/Mouse/GRCm39/genome.fa -gtf ~/data/genomes/Mouse/GRCm39/Mus_musculus.GRCm39.103.gtf > F-H3K4_1_annotated_regions.txt
annotatePeaks.pl F-H3K4_2_regions.txt ~/data/genomes/Mouse/GRCm39/genome.fa -gtf ~/data/genomes/Mouse/GRCm39/Mus_musculus.GRCm39.103.gtf > F-H3K4_2_annotated_regions.txt