#!/bin/bash

# Request Resources
#SBATCH --time=44:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=146G

# Specify Job name
#SBATCH -J STRINGTIE_TPM

# Specify Job Output and Input Files
#SBATCH -o outputfilestringtie.out
#SBATCH -e errorfilestringtie.err

module load stringtie


for aligned in *.bam; do
OUTPUT=$(basename ${aligned} |cut -f 1 -d ".");
echo "Output = $OUTPUT";

stringtie ${aligned} -o ~/scratch/gtfs/${OUTPUT}.gtf -G /users/mgura/data/genomes/Mouse/GRCm38/Mus_musculus.GRCm38.96.gtf -l ${OUTPUT} -A ~/scratch/tabs/${OUTPUT}.tab -e
done
