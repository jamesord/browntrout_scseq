#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="stringtie_trout_070423.sh"
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=150M

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate stringtie

ref_gtf=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/Salmo_trutta.fSalTru1.1.104.gtf

for bam in `ls ../../mapping/*.bam | xargs -n 1 basename`;do
sample=`echo ${bam} | cut -f1 -d'_'`
stringtie -o ${sample}_stringtie.gtf -G $ref_gtf -p 5 -l ${sample}_STRG -s 1000000000 ../../mapping/${bam}
done

ls *.gtf > gtf_list.txt
stringtie -o stringtie_trout_merged_070423.gtf --merge -F 3 -p 5 gtf_list.txt
stringtie -o stringtie_trout_refmerged_070423.gtf --merge -F 3 -p 5 gtf_list.txt -G $ref_gtf