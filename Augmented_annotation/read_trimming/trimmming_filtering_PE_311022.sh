#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="trimmming_filtering_PE_311022.sh"
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G

# load modules
module load vital-it/7
module load UHTS/Analysis/BBMap/38.91
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Aligner/bowtie2/2.3.4.1
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/pigz-2.6
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/mira_4.9.6_linux-gnu_x86_64_static/bin
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/samtools-1.15.1/bin/

#mkdir fastqc

# raw fastq files obtained from ENA using wget
for sample in `cat sample_list_pe`;do

gunzip ${sample}_1.fastq.gz
gunzip ${sample}_2.fastq.gz

# quality / adapter trim with bbduk
bbduk.sh in=${sample}_1.fastq in2=${sample}_2.fastq out=${sample}_1.trimmed.fastq out2=${sample}_2.trimmed.fastq outs=${sample}_s.trimmed.fastq \
ref=adapters ktrim=r k=25 mink=11 tpe=t tbo=t hdist=1 trimq=10 qtrim=rl threads=10

# mirabait
mirabait -I -t 10 -j rrna -p ${sample}_1.trimmed.fastq ${sample}_2.trimmed.fastq ${sample}_s.trimmed.fastq

# filter for T. bryosalmonae transcripts
bowtie2 -x /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/t_bryo/Tbryo_AhmadFaberKumar \
-1 rRNA_miss_${sample}_1.trimmed.fastq \
-2 rRNA_miss_${sample}_2.trimmed.fastq \
-U rRNA_miss_${sample}_s.trimmed.fastq \
--un-conc-gz ${sample}.trimmed.filtered.fastq.gz \
--un-gz ${sample}.trimmed.filtered.s.fastq.gz \
-p 5 --no-unal --very-sensitive | samtools sort -@ 5 -o ${sample}_tbryo.bam
rm ${sample}_tbryo.bam

# renaming to retain proper extension
mv ${sample}.trimmed.filtered.fastq.1.gz ${sample}.trimmed.filtered.1.fastq.gz
mv ${sample}.trimmed.filtered.fastq.2.gz ${sample}.trimmed.filtered.2.fastq.gz

# run fastqc on the trimmed reads
fastqc -t 10 -o fastqc ${sample}.trimmed.filtered.1.fastq.gz
fastqc -t 10 -o fastqc ${sample}.trimmed.filtered.2.fastq.gz

# cleanup
pigz -p 10 ${sample}_1.fastq
pigz -p 10 ${sample}_2.fastq
rm ${sample}_*.trimmed.fastq
rm rRNA_*_${sample}_*.trimmed.fastq

done