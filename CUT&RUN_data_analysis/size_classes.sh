#!/bin/bash

# Request Resources
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job Name
#SBATCH -J sizes_CUTRUN

#Specify Output and Input Files
#SBATCH -o outputsizes.out
#SBATCH -e errorsizes.err

module load samtools

for bam in *_filtsort.bam; do
OUTPUT=$(basename ${bam} | cut -f 1 -d "_")
samtools view -h -o ${OUTPUT}_filtsort.sam ${bam}

# get tf class sizes
cat <(head -n 68 ${OUTPUT}_filtsort.sam) <(awk ' $9 <= 120 && $9 >= 1 || $9 >= -120 && $9 <= -1 ' ${OUTPUT}_filtsort.sam) > ${OUTPUT}.1_120.sam
samtools view -h -t GRCm39_chrom.sizes -b -o ${OUTPUT}.1_120.bam ${OUTPUT}.1_120.sam

# get histone class sizes
cat <(head -n 68 ${OUTPUT}_filtsort.sam) <(awk ' $9 <= 500 && $9 >= 150 || $9 >= -500 && $9 <= -150 ' ${OUTPUT}_filtsort.sam) > ${OUTPUT}.150_500.sam
samtools view -h -t GRCm39_chrom.sizes -b -o ${OUTPUT}.150_500.bam ${OUTPUT}.150_500.sam

done
