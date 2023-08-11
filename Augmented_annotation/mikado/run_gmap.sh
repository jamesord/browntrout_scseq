#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="gmap_140423.sh"
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00

module load UHTS/Analysis/gmap/2019.03.15

ref_dir=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/
ref_fasta=Salmo_trutta.fSalTru1.1.dna.toplevel.fa
ref_name=Salmo_trutta.fSalTru1.1.dna.toplevel
okayset=/storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/tr2aacds/Aus_Est_pHetero7/okayset/all_transcripts_AusEst.okay.mrna

#cd $ref_dir

#gmap_build -D $ref_dir -d $ref_name $ref_fasta

#cd /storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/mikado/

#mkdir gmap
cd gmap

gmap -f 3 --nthreads 10 -D $ref_dir -d $ref_name $okayset > all_transcripts_AusEst.okay.gmap.gff3