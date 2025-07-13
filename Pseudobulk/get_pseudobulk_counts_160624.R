library(dplyr)
library("Seurat")
library(data.table)
library(Matrix)

# Define pseudobulk function
get_pseudobulks<-function(seurat,out.folder){
dir.create(out.folder)
clust_counts<-data.frame(cluster=Idents(seurat),individual=seurat$orig.ident,row.names = NULL)
clust_counts<-clust_counts[!duplicated(clust_counts),] %>% group_by(cluster) %>% summarise(individuals=n())
pseudobulk_clusters<-(droplevels(subset(clust_counts,individuals==nlevels(as.factor(seurat$orig.ident)))$cluster))
cell_index<-data.frame(cluster=Idents(seurat),individual=seurat$orig.ident)
cell_index<-subset(cell_index,cluster %in% pseudobulk_clusters)
for (i in pseudobulk_clusters){
  y<-data.frame(row.names=row.names(seurat@assays$RNA@counts))
  for (j in levels(as.factor(seurat$orig.ident))){
    cellIDs<-rownames(subset(cell_index,individual==j&cluster==i))
    selection<-seurat@assays$RNA@counts[,cellIDs] %>% rowSums(na.rm = FALSE, dims = 1) %>% as.data.frame
    colnames(selection)<-j
    y<-cbind(y,selection)
  }
  write.table(y,file=paste(out.folder,"/pseudobulk_counts_cluster",i,".txt",sep=""),sep="\t",quote=FALSE)
}
}

# load test data
load("../trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.rds")

immune.combined.sct <- subset(immune.combined.sct, idents = c("1","5","12","21"))

# get pseudobulks
get_pseudobulks(immune.combined.sct,"pseudobulk_counts_160624")
