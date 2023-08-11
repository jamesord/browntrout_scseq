#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="rnabloom_trout_SE_311022.sh"
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=20:00:00
#SBATCH --dependency=afterok:40829384

PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/samtools-1.15.1/bin/

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate rna-bloom

ref=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/ENSEMBL_cdna/Salmo_trutta.fSalTru1.1.cdna.all.fa

for sample in `cat ../../trimming/sample_list_se`;do
samtools fastq -@ 9 ../../mapping/${sample}_aligned_sorted.se.bam > ${sample}_aligned_sorted.single.fastq
done

cat *single.fastq > single.merged.fastq
rm *single.fastq

rnabloom -sef single.merged.fastq \
-t 10 -outdir rnabloom_trout_SE_311022 -ref $ref -norr

rm single.merged.fastq