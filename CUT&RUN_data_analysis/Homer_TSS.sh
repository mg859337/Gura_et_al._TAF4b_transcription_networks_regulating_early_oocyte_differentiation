#!/bin/bash

# Request Resources
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=246G

# Specify Job Name
#SBATCH -J Homer_CUTRUN_TSS

#Specify Output and Input Files
#SBATCH -o outputhomer_TSS.out
#SBATCH -e errorhomer_TSS.err

module load samtools
module load python/3.5.2

cd ~/Apps/
source Homer.venv/bin/activate
PATH=$PATH:/users/mgura/Apps/Homer.venv/bin/
cd Homer.venv/

annotatePeaks.pl tss GRCm39 -size 4000 -hist 20 -d ~/scratch/F-TAF4b_1_1.120/ > ~/scratch/F-TAF4b_1_1.120/F-TAF4b_1_TSS.txt
annotatePeaks.pl tss GRCm39 -size 4000 -hist 20 -d ~/scratch/F-TAF4b_2_1.120/ > ~/scratch/F-TAF4b_2_1.120/F-TAF4b_2_TSS.txt

annotatePeaks.pl tss GRCm39 -size 4000 -hist 20 -d ~/scratch/F-H3K4_1_150.500/ > ~/scratch/F-H3K4_1_150.500/F-H3K4_1_TSS.txt
annotatePeaks.pl tss GRCm39 -size 4000 -hist 20 -d ~/scratch/F-H3K4_2_150.500/ > ~/scratch/F-H3K4_2_150.500/F-H3K4_2_TSS.txt