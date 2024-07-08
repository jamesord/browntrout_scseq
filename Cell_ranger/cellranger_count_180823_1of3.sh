#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --job-name="cellranger_count_180823_1of3.sh"
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7G
#SBATCH --time=10:00:00

module load UHTS/SingleCell/cellranger/6.0.1

# sample list was generated as follows:
# ls *.fastq.gz | cut -c-2 | sort -u  > sample.list

for i in `cat /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/scRNAseq_121021/raw/sample.list | head -n 3`; do

# make a temporary folder and copy fastq files of one sample
mkdir ${i}_fastqs
cp /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/scRNAseq_121021/raw/${i}_*.fastq.gz ${i}_fastqs

# rename the files
cd ${i}_fastqs
for j in L1 L2 L3 L4;do
mv ${i}_${j}_R1_001_*.fastq.gz ${i}_S1_${j}_R1_001.fastq.gz
mv ${i}_${j}_R2_001_*.fastq.gz ${i}_S1_${j}_R2_001.fastq.gz
done
cd ../

cellranger count --fastqs ${i}_fastqs \
--id ${i}_cellranger \
--transcriptome /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/strutta_cellranger_120523/ \
--nosecondary --no-bam --include-introns

# remove temporary fastqs folder
rm -r ${i}_fastqs

done