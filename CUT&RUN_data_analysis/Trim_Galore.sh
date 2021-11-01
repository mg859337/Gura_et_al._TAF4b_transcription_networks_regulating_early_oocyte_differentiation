#!/bin/bash
# Request Resources
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job name
#SBATCH -J Trim_CR

# Specify Output and Input Files
#SBATCH -o outputtrim.out
#SBATCH -e errortrim.err

module load cutadapt
module load fastqc
module load trimgalore

for read1 in *R1_001.fastq.gz; do
# OUTPUT2=$(basename${read1} | cut -f 1 -d ".");
read2=$(echo $read1 | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g');
echo "I see $read1"
echo "I see $read2"
trim_galore --paired -q 10 --fastqc ${read1} ${read2}
done
