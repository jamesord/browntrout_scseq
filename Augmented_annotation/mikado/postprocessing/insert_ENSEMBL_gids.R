library(data.table)
library(stringr)

gids <- read.table("noveltrans_ENSgenes.txt",header=F)
gtf <- read.table("mikado.loci.refloc_noveltrans.gtf",sep="\t",quote="",header=F)

gids$V1 <- paste('gene_id "',gids$V1,'";',sep="")
gids$V2 <- paste('gene_id "',gids$V2,'";',sep="")

genes<-NULL
for(i in seq(1:nrow(gids))){
gene<-subset(gtf,V9 %like% gids[i,2])
gene$V9<-str_replace_all(gene$V9,gids[i,2],gids[i,1])
genes<-rbind(genes,gene)
}

write.table(genes,col.names=F,row.names=F,file="mikado.loci.refloc_noveltrans.ENS.gtf",quote=F,sep="\t")