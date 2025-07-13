# FIG 3 PROLIFERATING CELL MARKERS
# UPDATED 18/07/24

library(dplyr)
library(biomaRt)
library(ggplot2)
library(data.table)
library(Seurat)
library(tidyr)
library(gridExtra)
library(forcats)
library("ggforce")
library(grid)
library("ComplexUpset")
library(ggtext)
library(ggplotify)
library("ggrastr")
`%not_in%`<-Negate(`%in%`)

# Load in dataset
load("C:/Users/James/Documents/R/trout_local/trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.downsamp.rds")
DefaultAssay(immune.combined.sct.downsampled)<-"integrated"

levels(Idents(immune.combined.sct.downsampled))<-gsub("8.0","8",levels(Idents(immune.combined.sct.downsampled)))
levels(Idents(immune.combined.sct.downsampled))<-gsub("10.0","10",levels(Idents(immune.combined.sct.downsampled)))
levels(Idents(immune.combined.sct.downsampled))<-gsub("19.0","19",levels(Idents(immune.combined.sct.downsampled)))

# read in cluster labels
clabs<-read.table("Seurat_cluster_ID/current/cluster_classifications_120724.txt",sep="\t",header=T)
colnames(clabs)[3]<-"ident"
level_order<-levels(immune.combined.sct.downsampled)
clabs <- clabs[match(level_order, clabs$cluster), ]
clabs$Group<-ifelse(clabs$Group=="Monocytes / Macrophages","Monocytes /\nMacrophages",clabs$Group) # for plotting purposes
new.cluster.ids <- as.character(clabs$ident)
names(new.cluster.ids) <- levels(immune.combined.sct.downsampled)
immune.combined.sct.downsampled <- RenameIdents(immune.combined.sct.downsampled, new.cluster.ids)
# read in marker data
fshmk_small<-read.delim("cell_marker_compilation/fish_compilation/all_potential_markers_long_290124.txt",sep="\t",header=T)
fshmk_small<-subset(fshmk_small, gene_id%in%VariableFeatures(immune.combined.sct.downsampled))

# # Fig 3A Prolif markers
nlevels(as.factor(subset(fshmk_small,primary_ct=="Proliferating_cells")$gene_id)) # 3 markers present
# Add module score
immune.combined.sct.downsampled<-AddModuleScore(object = immune.combined.sct.downsampled,
                                                 features = list(subset(fshmk_small,primary_ct=="Proliferating_cells")$gene_id),
                                                 ctrl = 10, name = 'Prolif')
# make initial violin plot and extract data (by merging with cluster labeling / metadata)
Pplot<-VlnPlot(immune.combined.sct.downsampled,features="Prolif1")
Pplotdat<-merge(Pplot$data,clabs,by="ident")
# reorder some factor levels etc. to ensure correct plotting
Pplotdat$Group <-factor(Pplotdat$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))
Pplotdat$ident<-factor(Pplotdat$ident,levels=levels(as.factor(clabs$ident)))
clabs<-clabs[order(clabs$ident),]

# generate plot
Pplot<-ggplot(Pplotdat,aes(x=ident,y=Prolif1,fill=ident))+
  theme_classic()+
  rasterise(geom_point(position=position_jitter(0.25),size=0.05,alpha=0.2))+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  geom_hline(yintercept=0,lty="dashed")+NoLegend()+
  labs(x="",y="Module expression score",title="A",subtitle="Proliferating cell markers (mki67, pcna x2)")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(size=30,face = "plain"))+
  facet_grid(.~Group,scales="free_x",space="free")

# UMAP plots of proliferation markers

# PROLIFERATION
# There are two homologs of PCNA: ENSSTUG00000049031 and ENSSTUG00000034708
# ONLY NEED TO SHOW ONE --> note that other homolog showed similar expression
p1<-FeaturePlot(immune.combined.sct.downsampled,features=c("ENSSTUG00000034708"), min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="pcna*",x="",y="",subtitle="ENSSTUG00000034708 | N3,B4,T6,MP\nM1 (<60% / L2FC<1)")+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(plot.title = element_text(face = "italic"))
# mki67
p2<-FeaturePlot(immune.combined.sct.downsampled,features=c("ENSSTUG00000004378"), min.cutoff = "q10", max.cutoff = "q90",raster=T) +
  labs(title="mki67*",x="",y="",subtitle="ENSSTUG00000004378 | N3,B4,T6\nM1,MP (<60%)")+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(plot.title = element_text(face = "italic"))

# see also H2AZ1
p3<-FeaturePlot(immune.combined.sct.downsampled,features=c("ENSSTUG00000000660"), min.cutoff = "q10", max.cutoff = "q90",raster=T) +
  labs(x="",y="",title="H2AZ",subtitle="ENSSTUG00000000660 | N3,M1,B4,T6,MP")+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(plot.title = element_text(face = "italic"))

# There are two homologs of PCLAF: ENSSTUG00000001118 and ENSSTUG00000046947
# ONLY NEED TO SHOW ONE --> note that other homolog showed similar expression
p4<-FeaturePlot(immune.combined.sct.downsampled,features=c("ENSSTUG00000046947"), min.cutoff = "q10", max.cutoff = "q90",raster=T) +
  labs(title="pclaf",x="",y="",subtitle="ENSSTUG00000046947 | N3,B4,T6\nM1,MP (<60%)")+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(plot.title = element_text(face = "italic"))

F3B<-as.ggplot(grid.arrange(p1,p2,p3,p4,ncol=4,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))

# Generate composite figure
lay <- rbind(c(1,1,1),
             c(2,2,2))

grobz <- lapply(list(Pplot,F3B), ggplotGrob)

pdf("Figure3_UPDATED_170724.pdf", width = 16, height = 8)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()