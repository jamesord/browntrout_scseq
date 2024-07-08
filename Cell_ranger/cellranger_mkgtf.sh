#!/bin/bash

module load UHTS/SingleCell/cellranger/6.0.1

# we filter to keep protein coding genes, lncRNA, and Immunoglobulin / TCR gene biotypes that are present in the annotation
# It is necessary to perform at least some filtering (e.g. small RNAs) as reads assigned to multiple features will not be counted downstream.
cellranger mkgtf Salmo_trutta.fSalTru1.1.104.gtf Salmo_trutta.fSalTru1.1.104.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lncRNA \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_J_gene