#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="extract_splice_junctions_140423.sh"
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate portcullis

#cd /storage/homefs/jo20n766/trout_kidney_transcriptome/mapping/
#ls /storage/homefs/jo20n766/trout_kidney_transcriptome/mapping/*.bam > bamlist

#samtools merge -o merged.bam -@ 10 -b bamlist

#rm bamlist
#mv merged.bam /storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/mikado/
#cd /storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/mikado/

#samtools index merged.bam

portcullis full \
/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/Salmo_trutta.fSalTru1.1.dna.toplevel.fa \
merged.bam