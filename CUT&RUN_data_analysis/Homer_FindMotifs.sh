#!/bin/bash

# Request Resources
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=246G

# Specify Job Name
#SBATCH -J Homer_CUTRUN_Motifs

#Specify Output and Input Files
#SBATCH -o outputhomermotifs.out
#SBATCH -e errorhomermotifs.err

module load homer
module load samtools

# find Taf4b motif

findMotifsGenome.pl F-TAF4b_1_peaks.txt ~/data/genomes/Mouse/GRCm39 F-TAF4b_1_motifs/ -size 200 -mask
findMotifsGenome.pl F-TAF4b_2_peaks.txt ~/data/genomes/Mouse/GRCm39 F-TAF4b_2_motifs/ -size 200 -mask
