#!/bin/bash

# Request Resources
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job Name
#SBATCH -J Homer_CUTRUN_Tags

#Specify Output and Input Files
#SBATCH -o outputhomertags.out
#SBATCH -e errorhomertags.err

module load homer
module load samtools

# make the tag directories for homer

makeTagDirectory F-TAF4b_1.120/ F-TAF4b_1.1_120.bam
makeTagDirectory F-TAF4b_2.120/ F-TAF4b_2.1_120.bam

makeTagDirectory F-H3K4_1_150.500/ F-H3K4_1.150_500.bam
makeTagDirectory F-H3K4_2_150.500/ F-H3K4_2.150_500.bam

makeTagDirectory F-IgG_1_1.120/ F-IgG_1.1_120.bam
makeTagDirectory F-IgG_2_1.120/ F-IgG_2.1_120.bam
makeTagDirectory F-IgG_1_150.500/ F-IgG_1.150_500.bam
makeTagDirectory F-IgG_2_150.500/ F-IgG_2.150_500.bam
