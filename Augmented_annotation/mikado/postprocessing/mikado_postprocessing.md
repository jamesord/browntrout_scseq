# Postprocessing of MIKADO transcripts
The goal here was to extract only high quality, 'novel' (i.e. non-reference) transcripts from the Mikado output and convert the corresponding annotations into GTF files which can be simply appended to existing ENSEMBL gene annotations.
The workflow is as follows:

1) Obtain and filter transcript and gene IDs of (1) novel gene transcripts and (2) novel transcripts of reference genes
2) Filter the Mikado gff3 for the above transcripts
3) Convert filtered gff3 files to gtf and replace Mikado gene IDs with ENSEMBL IDs (in case of reference genes)
4) Functional annotation of novel genes

The two gtf files this generates ('mikado.loci.refloc_noveltrans.ENS.gtf' and 'mikado.loci.novloc_noveltrans.gtf') were appended to the S. trutta ENSEMBL gtf to obtain the 'augmented' annotation.
## 1) GET IDs OF NOVEL TRANSCRIPTS (OF BOTH NOVEL GENES AND REFERENCE GENES), AND FILTER ACCORDING TO MIKADO TAGS
First get a list of all transcripts overlapping reference ones. This includes the reference transcripts themselves. It excludes fusion transcripts (tag "f"). The output includes mikado transcript and gene ID and the ENSEMBL IDs (gene and transcript) of overlapping transcripts.
```
awk '$10>0' ../compare.tmap | grep -v "f" | cut -f4,5,1,2 > ref_overlap_transcripts.txt
```
Get a list of all transcripts that do not overlap with reference transcripts (tag "u"), and their mikado gene IDs
```
awk '$3=="u"' ../compare.tmap | cut -f4,5 > non_ref_overlap_transcripts.txt
```
Get a list of all non-reference transcripts, further filtered to include only multi-exon (tag 27), coding (tag 20) and complete ORF (tag 36) transcripts  
```
cat ../mikado.loci.metrics.tsv | grep -v "re_ENS" | awk '$27>1' | awk '$20>0' | awk '$36=="True"' | cut -f1 > nonref_transcripts_filtered.txt
```
The high quality non-reference transcripts in 'nonref_transcripts_filtered.txt' are used to filter the first two files ('ref_overlap_transcripts.txt' and 'non_ref_overlap_transcripts.txt') using the R script 'filter_transcripts1.R'.
```
module load R
Rscript filter_transcripts_1.R
```
The above R script produces the files: refloc_noveltrans.txt, novloc_noveltrans.txt, noveltrans_ENSgenes.txt (the latter used later for conversion of mikado IDs to ENSEMBL IDs).
### Counting the number of genes / transcripts
```
# number of transcripts derived from novel genes:
wc -l novloc_noveltrans.txt | cut -f1 -d" " # 1383
# number of novel genes:
cut -f2 novloc_noveltrans.txt | uniq | wc -l # 1200

# number of novel transcripts derived from reference genes:
wc -l refloc_noveltrans.txt | cut -f1 -d" " # 13712
# number of reference genes with novel transcripts:
cut -f2 refloc_noveltrans.txt | uniq | wc -l # 8812
```
## 2) FILTER THE MIKADO GFF3 ONLY FOR THE TANSCRIPTS OF INTEREST
Use the mikado grep utility to extract only the novel transcripts from the mikado.loci.gff3 file.
```
module load Anaconda3
eval "$(conda shell.bash hook)" 
conda activate mikado
mikado util grep novloc_noveltrans.txt ../mikado.loci.gff3 mikado.loci.novloc_noveltrans.gff3
mikado util grep refloc_noveltrans.txt ../mikado.loci.gff3 mikado.loci.refloc_noveltrans.gff3
conda deactivate
```
## 3) OBTAIN FINAL GTF FILES OF NOVEL TRANSCRIPTS
Use agat_convert_sp_gff2gtf.pl from the agat toolkit to convert the novel transcript gff3 files to gtf
```
conda activate agat
# FOR NOVEL GENES:
agat_convert_sp_gff2gtf.pl --gff mikado.loci.novloc_noveltrans.gff3 -o mikado.loci.novloc_noveltrans0.gtf
cat mikado.loci.novloc_noveltrans.gtf | sed 's/Name.*//g' | sed 's/Parent.*//g' > mikado.loci.novloc_noveltrans_cut.gtf
# remove old / rename new
rm mikado.loci.novloc_noveltrans.gtf
mv mikado.loci.novloc_noveltrans_cut.gtf mikado.loci.novloc_noveltrans.gtf

# FOR NOVEL TRANSCRIPTS OF REFERENCE GENES:
agat_convert_sp_gff2gtf.pl --gff mikado.loci.refloc_noveltrans.gff3 -o mikado.loci.refloc_noveltrans0.gtf
cat mikado.loci.refloc_noveltrans.gtf | sed 's/Name.*//g' | sed 's/Parent.*//g' | grep -v -w "gene" > mikado.loci.refloc_noveltrans_cut.gtf 
# note removal of the 'gene' rows as this is intended to be appended to the ENSEMBL gtf which already contains the gene-level entries
# remove old / rename new
rm mikado.loci.refloc_noveltrans.gtf
mv mikado.loci.refloc_noveltrans_cut.gtf mikado.loci.refloc_noveltrans.gtf
```
Then, for novel transcript of reference genes, replace MIKADO gene IDs with ENSEMBL gene IDs. This enables them to be integrated with the existing ENSEMBL annotations. See the R script for more details.
```
Rscript insert_ENSEMBL_gids.R
```
This outputs mikado.loci.refloc_noveltrans.ENS.gtf which, along with mikado.loci.novloc_noveltrans.gtf can now be APPENDED to the ENSEMBL gtf to get the AUGMENTED gtf.
## 4) FUNCTIONAL ANNOTATION
Transcripts assigned to novel genes can do with functional annotation (gene name corresponding with best BLAST hit, gene descriptions, GO terms).

First fet amino acid sequences for novel gene transcripts:
```
agat_sp_extract_sequences.pl -aa -g mikado.loci.novloc_noveltrans.gtf \
-f /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/Salmo_trutta.fSalTru1.1.dna.toplevel.fa \
-o mikado.loci.novloc_noveltrans.aa
```
Run blast to get the top hit for each aa -> use to assign gene names
```
module load Blast/ncbi-blast/latest
blastp -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
-num_threads 10 -query  mikado.loci.novloc_noveltrans.aa \
-db /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/databases/Trinotate/uniprot_sprot.pep \
-out mikado.loci.novloc_noveltrans.swissprot.blastp.outfmt6
```
Get only one entry per hit per gene; there may still be multiple hits per gene. Deduplicated list of top gene hits obtained as follows:
```
cut -f1 mikado.loci.novloc_noveltrans.swissprot.blastp.outfmt6 | rev | cut -f2- -d"."| rev > blast_hit_genes
cut -f2 mikado.loci.novloc_noveltrans.swissprot.blastp.outfmt6 > topphits
paste blast_hit_genes topphits | uniq > blastp_topgenehits_dedup.txt
```
The output file blastp_topgenehits_dedup.txt was then checked MANUALLY as one gene may still have more than one best hit.
The manually edited version (alternate hits assigned to an extra column) is blastphits_genelevel.txt.

FURTHER ANNOTATION (GO TERMS AND GENE DESCRIPTIONS) WAS PERFORMED USING PANZZER2 WEB INTERFACE
http://ekhidna2.biocenter.helsinki.fi/sanspanz/
The file mikado.loci.novloc_noveltrans.aa was submitted as a batch query
