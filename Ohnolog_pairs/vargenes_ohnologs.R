library(dplyr)
library(biomaRt)

# print more info about variable genes
vargenes<-read.table("Seurat_cluster_ID/variable_gene_descriptions_161023.txt",sep="\t",
                     header=T, quote = "",stringsAsFactors = FALSE)

# set up ensembl databases
ensembl_trout<-useMart("ensembl",dataset="strutta_gene_ensembl")
eluc_hom<-as.data.frame(getBM(attributes = c('ensembl_gene_id','elucius_homolog_ensembl_gene','elucius_homolog_orthology_confidence',"chromosome_name"),
                              filters = 'ensembl_gene_id',
                              values = vargenes$gene_id,
                              mart = ensembl_trout))
# subset for high confidence homology
eluc_hom<-subset(eluc_hom,elucius_homolog_orthology_confidence==1)
# count number of E.lucius homologs for each trout gene
eluc_hom <- eluc_hom %>% add_count(ensembl_gene_id,name="n_trout")
# and number of trout homologs for each E.lucius gene
eluc_hom <- eluc_hom %>% add_count(elucius_homolog_ensembl_gene,name="n_eluc")

# notice one possible "mistake" in that the homologs do not seem to be properly paired across species (look at the numbers)
# this is presumably just something that can happen and we need to account for it at the end.
subset(eluc_hom,elucius_homolog_ensembl_gene=="ENSELUG00000020684")
# how to account for this and get only properly assigned 2:1 pairings?

# subset 2:1 pairings (one copy in trout and 2 in pike)
ohnolog_pairs<-subset(eluc_hom,n_trout==1&n_eluc==2)[c(1,2,4)]
ohnolog_pairs<-ohnolog_pairs[order(ohnolog_pairs$elucius_homolog_ensembl_gene),]
# keeping only trout genes with 2x E.lucius IDs should address the abovementioned inconsistency
ohnolog_pairs<-ohnolog_pairs%>%add_count(elucius_homolog_ensembl_gene);ohnolog_pairs<-subset(ohnolog_pairs,n==2)[c(1:3)]

# get only pairs on different trout chromosomes (count the number of trout chromosome on which a pike homolog occurs)
ohnolog_pairs<-ohnolog_pairs %>% add_count(elucius_homolog_ensembl_gene,chromosome_name)
ohnolog_pairs<-subset(ohnolog_pairs,n==1)[c(1,2)]

write.table(ohnolog_pairs,file="variable_genes_ohnolog_pairs.txt",row.names = F,col.names = T,quote=F,sep="\t")

