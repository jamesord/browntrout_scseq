#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="mikado_prep_210423.sh.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G

module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate mikado

mikado prepare -p 20 --json-conf configuration.yaml --exclude-redundant
