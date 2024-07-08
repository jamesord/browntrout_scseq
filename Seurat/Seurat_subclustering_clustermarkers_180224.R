# SEURAT subclustering and subcluster markers
# James Ord 04/02/24
# UPDATED 18/02/24: - get markers for all clusters/subclusters using RNA as default assay (as recommended by Seurat devs)
#                   - perform DE within cell types (neutrophils, monocytes, t-cells, b-cells) 

# see subclustering resolutions script and resulting subcluster UMAPs for rationale of specific resolutions / subclustering / merging decisions.

library(dplyr)
library("Seurat")
library("patchwork")
library(data.table)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(fossil)
library("forcats")
library(biomaRt)
library(clustree)
library(sctransform)
library(DESeq2)
library("glmGamPoi")
library(future)

# set maximum shared memory (100GB)
options(future.globals.maxSize = 100 * 1024 ^ 3)

# load data
load("trout_SCTprocessed_PC30k30_var.ft.filt_161023.rds")

# create subcluster variable
immune.combined.sct$subcluster <- as.character(Idents(immune.combined.sct))

# CLUSTER 8: SPLIT INTO 3 SUBCLUSTERS
immune.combined.sct.c8 <- subset(immune.combined.sct, idents = 8)
immune.combined.sct.c8 <- FindNeighbors(immune.combined.sct.c8, dims = 1:30)
immune.combined.sct.c8 <- FindClusters(immune.combined.sct.c8,res=0.05)

# merge subcluster 8.3 into 8.2
new_names <- c("0", "1", "2", "2")
names(new_names) <- levels(immune.combined.sct.c8)
immune.combined.sct.c8 <- RenameIdents(object = immune.combined.sct.c8, new_names)

# Change the information of cells containing sub-cluster information
immune.combined.sct$subcluster[rownames(immune.combined.sct.c8@meta.data)] <- paste("8.",Idents(immune.combined.sct.c8),sep="")

# CLUSTER 10: SPLIT INTO 3 SUBCLUSTERS
immune.combined.sct.c10 <- subset(immune.combined.sct, idents = 10)
immune.combined.sct.c10 <- FindNeighbors(immune.combined.sct.c10, dims = 1:30)
immune.combined.sct.c10 <- FindClusters(immune.combined.sct.c10,res=0.05)

# merge subcluster 10.3 into 10.1
new_names <- c("0", "1", "2", "1")
names(new_names) <- levels(immune.combined.sct.c10)
immune.combined.sct.c10 <- RenameIdents(object = immune.combined.sct.c10, new_names)

# Change the information of cells containing sub-cluster information
immune.combined.sct$subcluster[rownames(immune.combined.sct.c10@meta.data)] <- paste("10.",Idents(immune.combined.sct.c10),sep="")

# CLUSTER 19: SPLIT INTO TWO SUBCLUSTERS
immune.combined.sct.c19 <- subset(immune.combined.sct, idents = 19)
immune.combined.sct.c19 <- FindNeighbors(immune.combined.sct.c19, dims = 1:30)
immune.combined.sct.c19 <- FindClusters(immune.combined.sct.c19,res=0.02)
# Change the information of cells containing sub-cluster information
immune.combined.sct$subcluster[rownames(immune.combined.sct.c19@meta.data)] <- paste("19.",Idents(immune.combined.sct.c19),sep="")

# reset default ident to subcluster
immune.combined.sct<-SetIdent(immune.combined.sct, value = immune.combined.sct@meta.data$subcluster)

write.table(immune.combined.sct$subcluster,file="subcluster_IDs_040224.txt",row.names=T,sep="\t",quote=F,col.names=F)

###################################
## COUNT CELLS IN CLUSTERS, ECT. ##
###################################

# count up the numbers of cells in clusters, individuals etc.
clust_counts<-data.frame(cluster=Idents(immune.combined.sct),
                         individual=immune.combined.sct$orig.ident,
                         group=immune.combined.sct$group,
                         batch=immune.combined.sct$batch) %>% 
  add_count(cluster,name="cluster_count") %>% add_count(individual,name="individual_count") %>%
  add_count(individual,cluster,name="cluster_individual_count")
clust_counts<-clust_counts[!duplicated(clust_counts),]

# write table of cell counts in individuals / clusters (filtered)
write.table(clust_counts,file="cluster_individual_counts_subclusters_040224.txt",sep="\t",row.names=FALSE,quote=FALSE)

##########################
## GETTING MARKER GENES ##
########################## 

# run DE only on the variable features

immune.combined.sct1<-subset(x = immune.combined.sct, features = VariableFeatures(object = immune.combined.sct))
DefaultAssay(immune.combined.sct1) <- "RNA"
immune.combined.sct1 <- NormalizeData(immune.combined.sct1)

# get the markers # important to set assay as RNA, see e.g., https://github.com/satijalab/seurat/discussions/4000
markers <- FindAllMarkers(immune.combined.sct1, only.pos = TRUE, min.pct = 0.25, assay="RNA") # use 25% in case of clusters containing multiple cell types or subtypes
colnames(markers)[7]<-"gene_id"
# write out
write.table(markers,file="cluster_markers_final_180224.txt",sep="\t",row.names=FALSE,quote=FALSE)

# GET DE GENES / MARKERS AMONGST CELL TYPES --> omit 'only.pos = TRUE'
# Neutrophils
immune.combined.sct.N <- subset(immune.combined.sct1, idents = c(1,2,3,6,20,27))
markersN <- FindAllMarkers(immune.combined.sct.N, min.pct = 0.25, assay="RNA")
colnames(markersN)[7]<-"gene_id"
write.table(markersN,file="DE_Neutrophils_180224.txt",sep="\t",row.names=FALSE,quote=FALSE)
# Monocytes
immune.combined.sct.M <- subset(immune.combined.sct1, idents = c(5,"10.0","10.1","10.2","19.0","19.1"))
markersM <- FindAllMarkers(immune.combined.sct.N, min.pct = 0.25, assay="RNA")
colnames(markersM)[7]<-"gene_id"
write.table(markersM,file="DE_Monocytes_180224.txt",sep="\t",row.names=FALSE,quote=FALSE)
# T-cells
immune.combined.sct.T <- subset(immune.combined.sct1, idents = c("8.0","8.1","8.2",11,21,26,28))
markersT <- FindAllMarkers(immune.combined.sct.T, min.pct = 0.25, assay="RNA")
colnames(markersT)[7]<-"gene_id"
write.table(markersT,file="DE_Tcells_180224.txt",sep="\t",row.names=FALSE,quote=FALSE)
# B-cells
immune.combined.sct.B <- subset(immune.combined.sct1, idents = c(0,7,12,13,15,22,24))
markersB <- FindAllMarkers(immune.combined.sct.B, min.pct = 0.25, assay="RNA")
colnames(markersB)[7]<-"gene_id"
write.table(markersB,file="DE_Bcells_180224.txt",sep="\t",row.names=FALSE,quote=FALSE)

#####################################
## WRITE OUT UPDATED SEURAT OBJECT ##
#####################################

DefaultAssay(immune.combined.sct) <- "integrated"

# get a subsample of the (var features-subset) seurat object (10%) for optional visualisation / testing on a local machine
immune.combined.sct.downsampled <- immune.combined.sct[, sample(colnames(immune.combined.sct), size = ncol(immune.combined.sct)/10, replace=F)]
immune.combined.sct.downsampled<-subset(x = immune.combined.sct.downsampled, features = VariableFeatures(object = immune.combined.sct.downsampled))
save(immune.combined.sct.downsampled, file = "trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.downsamp.rds")

# write out full updated object
save(immune.combined.sct, file = "trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.rds")