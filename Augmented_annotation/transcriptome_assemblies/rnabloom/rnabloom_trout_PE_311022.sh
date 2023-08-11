#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="rnabloom_trout_PE_311022.sh"
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=60:00:00

module load UHTS/Analysis/sambamba/0.7.1
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/samtools-1.15.1/bin/

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate rna-bloom

ref=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/ENSEMBL_cdna/Salmo_trutta.fSalTru1.1.cdna.all.fa

for sample in `cat ../../trimming/sample_list_pe`;do

# first extract paored reads which aligned concordantly, discarding those pairs which did not
sambamba view -F "proper_pair" -t 10 -f bam -o ${sample}_aligned_sorted.pp.bam ../../mapping/${sample}_aligned_sorted.pe.bam
sambamba sort -n -t 10 -o ${sample}_aligned_namesort.pp.bam ${sample}_aligned_sorted.pp.bam
samtools fastq -@ 9 -1 ${sample}_aligned_sorted.1.fastq -2 ${sample}_aligned_sorted.2.fastq ${sample}_aligned_namesort.pp.bam
rm ${sample}_aligned_sorted.pp.bam ${sample}_aligned_namesort.pp.bam

# however, we give the benefit of the doubt to aligned reads which were previously orphaned during quality trimming...
sambamba view -F "not paired" -t 10 -f bam -o ${sample}_aligned_sorted.se.bam ../../mapping/${sample}_aligned_sorted.pe.bam
samtools fastq -@ 9 ${sample}_aligned_sorted.se.bam > ${sample}_aligned_sorted.s.fastq
rm ${sample}_aligned_sorted.se.bam

done

cat *1.fastq > 1.merged.fastq
cat *2.fastq > 2.merged.fastq
cat *s.fastq > s.merged.fastq
rm *.1.fastq *.2.fastq *.s.fastq *.bai

rnabloom -left 1.merged.fastq -right 2.merged.fastq -revcomp-right -sef s.merged.fastq \
-t 10 -outdir rnabloom_trout_PE_311022 -ref $ref -norr

rm 1.merged.fastq 2.merged.fastq s.merged.fastq