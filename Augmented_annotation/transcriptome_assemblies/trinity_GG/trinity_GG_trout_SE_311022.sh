#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="trinity_GG_trout_SE_311022.sh"
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G

# optimising performance
# https://trinityrnaseq.github.io/performance/mem.html#:~:text=By%20default%2C%20Trinity%20runs%20two,the%20number%20of%20CPU's%20available.

# load modules
module load vital-it/7
module load Java/11.0.2
module load UHTS/Analysis/sambamba/0.7.1
module load UHTS/Analysis/jellyfish/2.3.0
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/samtools-1.15.1/bin/
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/salmon-1.9.0_linux_x86_64/bin/
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/bowtie2-2.4.5-linux-x86_64/
PATH=$PATH:/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/trinityrnaseq-v2.14.0

# set path to BAM files
bam_path=/storage/homefs/jo20n766/trout_kidney_transcriptome/mapping
# set output path on the workspace as we need more allowance for temporary files
out_path=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/temp

# make list of BAMs containing single read alignments and merge those BAMs
ls $bam_path/*.se.bam > bam.list.se
samtools merge -@ 15 -o merged.se.bam -b bam.list.se

# run trinity with full cleanup enabled
Trinity --genome_guided_bam merged.se.bam \
--genome_guided_max_intron 10000 --CPU 16 --max_memory 60G --output $out_path/Trinity_trout_SE_311022 --full_cleanup --bflyCalculateCPU

# remove merged bam
rm merged.se.bam bam.list.se

# recover output files from the workspace
mv $out_path/Trinity_trout_SE* ./