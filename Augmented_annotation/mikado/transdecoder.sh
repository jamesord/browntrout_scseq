#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="transdecoer_210423.sh.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --dependency=afterok:53304729

module load UHTS/Assembler/TransDecoder/5.3.0

TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta --cpu 20