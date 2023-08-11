mikado configure \
--list list.txt \
--reference /storage/workspaces/vetsuisse_fiwi_epev/fiwi_epev001/datasets/dataset_collection/trout/transcriptome/Salmo_trutta.fSalTru1.1.dna.toplevel.fa \
--mode permissive \
--scoring mammalian.yaml  \
--copy-scoring mammalian.yaml \
--junctions /storage/homefs/jo20n766/trout_kidney_transcriptome/assembly/mikado/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed \
-bt uniprot_sprot_plants.fasta \
configuration.yaml