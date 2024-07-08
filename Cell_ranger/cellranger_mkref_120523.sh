#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --job-name="cellranger_mkref_120523.sh"
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=02:00:00

module load UHTS/SingleCell/cellranger/6.0.1

cellranger mkref --genome=strutta_cellranger_120523 --fasta=Salmo_trutta.fSalTru1.1.dna.toplevel.fa --genes=Salmo_trutta.fSalTru1.1.104.filtered.mikado-augmented.gtf --memgb=80