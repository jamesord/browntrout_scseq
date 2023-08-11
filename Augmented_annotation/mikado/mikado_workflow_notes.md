# MIKADO workflow for S. trutta transcripts
James Ord | notes compiled 11/08/23

The following workflow used Mikado to pick candidate transcripts for augmenting the ENSEMBL S. trutta annotation.
Used Mikado v2.3.4 installed via conda. It follows loosely the tutorial workflow: https://mikado.readthedocs.io/en/stable/Tutorial/

Transcriptome assemblies were obtained from trout kidney RNAseq reads using Trinity (genome-guied), RNAbloom and Stringtie. High-quality splice junctions were obtained using portcullis, using the a merged BAM file of all RNAseq alignments.
The Trinity and RNAbloom contigs were combined and collapsed (redundant transcripts removed) using tr2aacds. The tr2aacds-collapsed transcriptome was aligned to the genome (fSalTru1.1) using GMAP.
Subsequently, three transcript assemblies were inputted to the Mikado workflow: tr2aacds, stringtie, and the ENSEMBL transcripts (Salmo_trutta.fSalTru1.1.104.gtf) (see list.txt).

The following were then run:
1) mikado_configure.sh
2) mikado_prep.sh
3) transdecoder.sh (run transdecoder on mikado_prepared.fasta)
4) mikado_serialise.sh
5) mikado_pick.sh
6) mikado_compare.sh

The resulting transcripts were outputted in the file 'mikado.loci.gff3'. These were subsequently filtered further to obtain only non-reference transcripts which met specific quality criteria (coding transcripts, complete ORFs, multiple exons).
The output files compare.tmap and mikado.loci.metrics.tsv were used for this filtering.
Further details of the filtering and postprocessing are in the folder 'postprocessing'.
