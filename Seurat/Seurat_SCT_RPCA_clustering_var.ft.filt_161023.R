# SEURAT SCT TRANSFORM WITH RPCA AND MULTIPLE CLUSTERING RESOLUTIONS
# Using 30 PCs and K anchor of 30
# James Ord 16/10/23

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
library("glmGamPoi")
library(future)

# set maximum shared memory (300GB)
options(future.globals.maxSize = 300 * 1024 ^ 3)
options(bitmapType='cairo')

############################################################
## 1: Fetching data, generating and merging Seurat object ##
############################################################

rm(list=ls())

dirs<-list.dirs("input_data")[2:10]
batch.data <- Read10X(data.dir = dirs, gene.column = 1)
gcdata <- CreateSeuratObject(counts = batch.data, project = "trout", min.cells = 3, min.features = 200)
# add metadata
levels(gcdata$orig.ident)<-c("A4","A6","A8","B4","B6","B8","C4","C6","C8")
# add group identities
gcdata$batch <- ifelse(gcdata$orig.ident%like%"A","Batch_A",
                                   ifelse(gcdata$orig.ident%like%"B","Batch_B","Batch_C"))
gcdata$group<-ifelse(gcdata$orig.ident%like%"4","wld",
                                   ifelse(gcdata$orig.ident%like%"6","mix","frm"))

################################################################
## 2. Normalisation / variance stabilisation with SCTransform ##
################################################################

# read in features to filter from variable features (cell cycle regulators and hemoglobins)
ccyc<-read.table("cell_cycle.txt")
hemos<-read.table("hemoglobins.txt")
to_remove<-rbind(ccyc,hemos)

# store mitochondrial percentage in object meta data
gcdata <- PercentageFeatureSet(gcdata, pattern = "^MT-", col.name = "percent.mt")
# NOTE (10/10/2024): this was not the correct way to store the % MT in this data; resulting variable was empty

# the rest is following the example here:
#https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1

gcdata.list <- SplitObject(gcdata, split.by = "batch")
gcdata.list <- lapply(X = gcdata.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = c("percent.mt")) # vars.to.regress not in example
# NOTE (10/10/2024): as percent.mt was empty, regressing on it would have done NOTHING. However this should not have adversely affected anything else.

features <- SelectIntegrationFeatures(object.list = gcdata.list)
gcdata.list <- PrepSCTIntegration(object.list = gcdata.list, anchor.features = features) #NOTE: subsets scale.data slot to only include anchor features

gcdata.list <- lapply(X = gcdata.list, FUN = RunPCA, features = features[! features %in% to_remove$V1])

# 30 PCs, k=30

immune.anchors <- FindIntegrationAnchors(object.list = gcdata.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 30)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30, features.to.integrate = features)
immune.combined.sct <- RunPCA(immune.combined.sct,features=VariableFeatures(immune.combined.sct)[! VariableFeatures(immune.combined.sct) %in% to_remove$V1])
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- RunTSNE(immune.combined.sct, reduction = "pca", dims = 1:30)

# INITIAL PLOTS
elbow_plot<-ElbowPlot(immune.combined.sct)
PCA_plot<-DimPlot(immune.combined.sct, reduction = "pca", dims = c(1, 2), group.by = "batch", raster=T)
UMAP_batches <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "batch", raster=T)+
  ggtitle("30 PCs, k anchor 30")
TSNE_batches <- DimPlot(immune.combined.sct, reduction = "tsne", group.by = "batch", raster=T)+
  ggtitle("30 PCs, k anchor 30")

immune.combined.sct0<-immune.combined.sct

################
## CLUSTERING ##
################

immune.combined.sct <- FindNeighbors(immune.combined.sct, dims = 1:30)

# FINDING THE CLUSTERING RESOLUTION
# Perform clustering at a series of resolutions from 0 to 1
for (i in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)){
  immune.combined.sct <- FindClusters(immune.combined.sct, resolution = i, algorithm = 1, random.seed = 100)}

# Generate clustering tree
treeplot<-clustree(immune.combined.sct, prefix = "integrated_snn_res.")

cluster_sc3s<-as.data.frame(treeplot$data %>% group_by(integrated_snn_res.) %>% summarise(mean_sc3=mean(sc3_stability)))
most_stable_res<-subset(cluster_sc3s,mean_sc3==max(cluster_sc3s$mean_sc3))[,1]
print(paste("Most stable resolution determined to be:",most_stable_res))

# Make a plot of SC3 index values
sc3_plot<-ggplot(cluster_sc3s,aes(x=integrated_snn_res.,y=mean_sc3,group = 1))+
  theme_bw()+
  geom_line()+
  geom_point()+
  labs(x="Clustering resolution",y="SC3 stability index")

# UMAP plots at various resolutions  
umap0.1<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.1",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.1")
umap0.2<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.2",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.2")
umap0.3<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.3",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.3")
umap0.4<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.4",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.4")
umap0.5<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.5",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.5")
umap0.6<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.6",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.6")
umap0.7<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.7",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.7")
umap0.8<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.8",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.8")
umap0.9<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.9",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.0.9")
umap1<-DimPlot(immune.combined.sct, reduction = "umap", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.1",raster=T) + NoLegend()+
  ggtitle("UMAP plot: integrated_snn_res.1")

# TSNE plots at various resolutions  
tsne0.1<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.1",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.1")
tsne0.2<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.2",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.2")
tsne0.3<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.3",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.3")
tsne0.4<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.4",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.4")
tsne0.5<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.5",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.5")
tsne0.6<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.6",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.6")
tsne0.7<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.7",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.7")
tsne0.8<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.8",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.8")
tsne0.9<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.0.9",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.0.9")
tsne1<-DimPlot(immune.combined.sct, reduction = "tsne", label = TRUE,label.box=TRUE,group.by = "integrated_snn_res.1",raster=T) + NoLegend()+
  ggtitle("tsne plot: integrated_snn_res.1")

pdf("DimPlots_PC30k30_var.ft.filt_161023.pdf")
print(PCA_plot)
print(elbow_plot)
print(UMAP_batches)
print(TSNE_batches)
print(treeplot)
print(sc3_plot)
print(umap0.1)
print(tsne0.1)
print(umap0.2)
print(tsne0.2)
print(umap0.3)
print(tsne0.3)
print(umap0.4)
print(tsne0.4)
print(umap0.5)
print(tsne0.5)
print(umap0.6)
print(tsne0.6)
print(umap0.7)
print(tsne0.7)
print(umap0.8)
print(tsne0.8)
print(umap0.9)
print(tsne0.9)
print(umap1)
print(tsne1)
dev.off()

# save the object; clustering can be performed again after the 'best' resolution is identified
save(immune.combined.sct0, file = "trout_SCTprocessed_PC30k30_var.ft.filt_161023.rds")
