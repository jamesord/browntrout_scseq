#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="mikado_serialise_210423.sh.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --dependency=afterok:53304741

module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate mikado

mikado serialise --json-conf configuration.yaml --orfs mikado_prepared.fasta.transdecoder.gff3 \
--junctions /storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/mikado/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed