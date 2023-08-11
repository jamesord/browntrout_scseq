#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="mikado_pick_210423.sh.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --dependency=afterok:53304760

module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate mikado

mikado pick -p 20 --configuration configuration.yaml --subloci-out mikado.subloci.gff3 \
--genome /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/Salmo_trutta.fSalTru1.1.dna.toplevel.fa \
--reference-update