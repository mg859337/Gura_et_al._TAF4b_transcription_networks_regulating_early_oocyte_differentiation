#!/bin/bash

# Request resources: 
#SBATCH --time=40:00:00 
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G

# Specify a job name: 
#SBATCH -J cellranger_aggr_Zhao

# Specify output and error files: 
#SBATCH -o outputfileaggr.out
#SBATCH -e errorfileaggr.err

module load python/3.5.2
cd ..
source cellranger.venv/bin/activate

cd cellranger.venv/
cd cellranger-5.0.0/

export PATH=/users/mgura/scratch/cellranger.venv/cellranger-5.0.0:$PATH

cellranger aggr --id=Zhao_aggr \
					--csv=/users/mgura/scratch/Zhao/Zhao_libraries.csv \
					--normalize=mapped 
