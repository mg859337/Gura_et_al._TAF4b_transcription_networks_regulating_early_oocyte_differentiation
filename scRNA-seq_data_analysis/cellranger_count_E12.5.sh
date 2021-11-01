#!/bin/bash

# Request resources: 
#SBATCH --time=40:00:00 
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G

# Specify a job name: 
#SBATCH -J cellranger_count

# Specify output and error files: 
#SBATCH -o outputfilecount_E12.out
#SBATCH -e errorfilecount_E12.err

module load python/3.5.2
cd ..
source cellranger.venv/bin/activate

cd cellranger.venv/
cd cellranger-5.0.0/

export PATH=/users/mgura/scratch/cellranger.venv/cellranger-5.0.0:$PATH

cellranger count 	--id=E12_germ_cells \
					--transcriptome=../refdata-gex-mm10-2020-A \
					--fastqs=/users/mgura/scratch/Zhao \
					--sample=E12
