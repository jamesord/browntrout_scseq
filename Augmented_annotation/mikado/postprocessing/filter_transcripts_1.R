library(data.table)
nonreftr<-read.table("nonref_transcripts_filtered.txt",header=F)
ref_overlap<-read.table("ref_overlap_transcripts.txt",header=F)
novloctr<-read.table("non_ref_overlap_transcripts.txt",header=F)

refgene_noveltrans<-subset(ref_overlap,V3%in%nonreftr$V1)
write.table(refgene_noveltrans[c(3,4)],file="refloc_noveltrans.txt",row.names=F,col.names=F,quote=F,sep="\t")

novloctr<-subset(novloctr,V1%in%nonreftr$V1)
write.table(novloctr,file="novloc_noveltrans.txt",row.names=F,col.names=F,quote=F,sep="\t")

# for the conversion of novel transcript gene IDs to ENSEMBL
gids<-refgene_noveltrans[c(2,4)];gids<-gids[!duplicated(gids),]
write.table(gids,file="noveltrans_ENSgenes.txt",col.names=F,row.names=F,quote=F,sep="\t")
