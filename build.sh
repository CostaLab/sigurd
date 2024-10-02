#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --job-name=test
#SBATCH --time=10-10:00:00
#SBATCH --output=/data/mg000001/logs/output.%J.%x.txt

module unload R
module load scRNA

Rscript /data/mg000001/sigurd/build.R
