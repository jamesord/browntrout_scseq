#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="blastp_120523.sh"
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=250M

module load Blast/ncbi-blast/latest

blastp -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
-num_threads 10 -query  mikado.loci.novloc_noveltrans.aa -db /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/databases/Trinotate/uniprot_sprot.pep -out mikado.loci.novloc_noveltrans.swissprot.blastp.outfmt6