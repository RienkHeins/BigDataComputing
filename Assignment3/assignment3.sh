#!/bin/bash
#SBATCH --job-name=assignment3
#SBATCH --output=out.txt
#SBATCH --error=err.txt
#SBATCH --time 10:00:00
#SBATCH --cpus-per-task=17
#SBATCH --nodes=1
#SBATCH --partition=assemblix

# load conda environment and create variables for data files
source /commons/conda/conda_load.sh
export INDEX=/data/dataprocessing/MinIONData/all_bacteria.fna
export DATA=/data/dataprocessing/MinIONData/all.fq

# run minimap
for n in {1..16} ; do /usr/bin/time -o timings.txt --append -f "${n}\t%e" minimap2  -a -t ${n} $INDEX $DATA > /dev/null ; done
