#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="hisat2_trout_PE_311022.sh"
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2G

# load modules
module load vital-it/7
module load UHTS/Aligner/hisat/2.2.1
module load UHTS/Analysis/sambamba/0.7.1
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/samtools-1.15.1/bin/

ref_path=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/

for sample in `cat ../trimming/sample_list_pe`; do

hisat2 -x $ref_path/Salmo_trutta.fSalTru1.1.dna.toplevel \
-p 8 --known-splicesite-infile $ref_path/Salmo_trutta.fSalTru1.1.104.filtered.splices.txt \
-1 ../trimming/${sample}.trimmed.filtered.1.fastq.gz \
-2 ../trimming/${sample}.trimmed.filtered.2.fastq.gz \
-U ../trimming/${sample}.trimmed.filtered.s.fastq.gz \
| samtools sort -@ 7 -o ${sample}_sorted.pe.bam
sambamba view -F "not unmapped" -t 15 -f bam -o ${sample}_aligned_sorted.pe.bam ${sample}_sorted.pe.bam

rm ${sample}_sorted.pe.bam

done