#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="mikado_compare_210423.sh.sh"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --dependency=afterok:53304783

module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate mikado

mikado compare -r Salmo_trutta.fSalTru1.1.104.gtf --index -x 20

mikado compare -r Salmo_trutta.fSalTru1.1.104.gtf -p mikado.loci.gff3 -o compare -l compare.log -x 20