#!/bin/bash

# Request Resources
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=246G

# Specify Job Name
#SBATCH -J Homer_CUTRUN_Peaks

#Specify Output and Input Files
#SBATCH -o outputhomerpeaks.out
#SBATCH -e errorhomerpeaks.err

module load homer
module load samtools

# call peaks


findPeaks ../F-TAF4b_1_1.120/ -style factor -P 0.0015 -fdr 0.01 -F 2.0 -L 2 -LP 0.0015 -i ../F-IgG_1_1.120/ -o ./F-TAF4b_1_peaks.txt
findPeaks ../F-TAF4b_2_1.120/ -style factor -P 0.0015 -fdr 0.01 -F 2.0 -L 2 -LP 0.0015 -i ../F-IgG_2_1.120/ -o ./F-TAF4b_2_peaks.txt

findPeaks ../F-H3K4_1_150.500/ -style histone -P 0.0015 -fdr 0.01 -F 2.0 -L 2 -LP 0.0015 -i ../F-IgG_1_150.500/ -o ./F-H3K4_1_peaks.txt
findPeaks ../F-H3K4_2_150.500/ -style histone -P 0.0015 -fdr 0.01 -F 2.0 -L 2 -LP 0.0015 -i ../F-IgG_2_150.500/ -o ./F-H3K4_2_peaks.txt
