#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=james.ord@vetsuisse.unibe.ch
#SBATCH --job-name="tr2aacds_trout_pHetero7_141122.sh"
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G

ml Blast/ncbi-blast/2.10.1+
ml UHTS/Analysis/cd-hit/4.6.8
ml SequenceAnalysis/SequenceAlignment/exonerate/2.4.0
ml UHTS/Analysis/SeqKit/0.13.2
export evigene=/storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/software/evigene/

mkdir Aus_Est_pHetero7

# concatenate assemblies
cat ../rnabloom/rnabloom_trout_SE_311022/rnabloom.transcripts.fa \
../trinity/Trinity_trout_SE_311022.Trinity-GG.fasta \
../rnabloom/rnabloom_trout_PE_311022/rnabloom.transcripts.fa \
../trinity/Trinity_trout_PE_311022.Trinity-GG.fasta | awk '{print $1}' | awk '/^>/{$0=$0"_"(++i)}1' > Aus_Est_pHetero7/all_transcripts_AusEst.fasta

cd Aus_Est_pHetero7
$evigene/scripts/prot/tr2aacds.pl -species=Salmo_trutta -NCPU=16 -MAXMEM=160000 -cdnaseq all_transcripts_AusEst.fasta -pHeterozygosity 7

grep 'main\|noclass' okayset/all_transcripts_AusEst.pubids | cut -f1 -d$'\t' > okayset/primary_IDs.txt
seqkit grep -f okayset/primary_IDs.txt -j 16 okayset/all_transcripts_AusEst.okay.mrna > okayset/all_transcripts_AusEst.okay.primary.mrna

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate busco
busco -m transcriptome --cpu 16 -i okayset/all_transcripts_AusEst.okay.primary.mrna -o busco_primaries -l actinopterygii_odb10 --download_path $WORKSPACE/databases/busco_downloads/