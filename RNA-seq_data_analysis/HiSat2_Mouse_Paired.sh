#!/bin/bash

# The following specifies computing cluster parameters

# Request Resources
#SBATCH --time=44:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job name
#SBATCH -J HISAT2_RNAseq

# Specify Job Output and Input Files
#SBATCH -o outputfilealign.out
#SBATCH -e errorfilealign.err

module load hisat2
module load samtools

DBDIR="/users/mgura/data/genomes/Mouse/GRCm38"
GENOME="genome"

for fastq in *R1_001_val_1.fq.gz; do
	fastq2=$(echo $fastq | sed 's/R1_001_val_1.fq.gz/R2_001_val_2.fq.gz/g'); 
	echo "${fastq} ${fastq2}";
	OUTPUT=$(basename ${fastq} |cut -f 1 -d "_"); 
	echo "Output = $OUTPUT";

	hisat2 \
		--dta \
		-x ${DBDIR}/${GENOME} \
		-1 ${fastq} \
		-2 ${fastq2} \
		-S ${OUTPUT}.sam &> ${OUTPUT}.log

	samtools sort ${OUTPUT}.sam-o ${OUTPUT}.bam ${OUTPUT}.sam
  rm ${OUTPUT}.sam

done
