#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="trimmming_filtering_SE_311022.sh"
#SBATCH --time=07:00:00
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

mkdir fastqc

# raw fastq files obtained from ENA using wget
for sample in `cat sample_list_se`;do

gunzip ${sample}.fastq.gz

# quality / adapter trim with bbduk
bbduk.sh in=${sample}.fastq out=${sample}.trimmed.fastq ref=adapters ktrim=r k=25 mink=11 tpe=t tbo=t hdist=1 trimq=10 qtrim=rl threads=10

# mirabait
mirabait -I -t 10 -j rrna ${sample}.trimmed.fastq

# filter for T. bryosalmonae transcripts
bowtie2 -x /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/t_bryo/Tbryo_AhmadFaberKumar -U rRNA_miss_${sample}.trimmed.fastq \
--un-gz ${sample}.trimmed.filtered.fastq.gz -p 5 --no-unal --very-sensitive | samtools sort -@ 5 -o ${sample}_tbryo.bam
rm ${sample}_tbryo.bam

# run fastqc on the trimmed reads
fastqc -t 10 -o fastqc ${sample}.trimmed.filtered.fastq.gz

pigz -p 10 ${sample}.fastq
rm ${sample}.trimmed.fastq
rm rRNA_miss_${sample}.trimmed.fastq
rm rRNA_match_${sample}.trimmed.fastq

done