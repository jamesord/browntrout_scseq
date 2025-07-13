# FIGURES 1 AND 2
library(dplyr)
library(biomaRt)
library(ggplot2)
library(data.table)
library(Seurat)
library(tidyr)
library(gridExtra)
library(forcats)

# UPDATED 17/07/24

load("../trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.rds")
DefaultAssay(immune.combined.sct)<-"integrated"

levels(Idents(immune.combined.sct))<-gsub("8.0","8",levels(Idents(immune.combined.sct)))
levels(Idents(immune.combined.sct))<-gsub("10.0","10",levels(Idents(immune.combined.sct)))
levels(Idents(immune.combined.sct))<-gsub("19.0","19",levels(Idents(immune.combined.sct)))

fshmk_small<-read.table("all_potential_markers_long_290124.txt",sep="\t",header=T)
fshmk_small<-subset(fshmk_small, gene_id%in%VariableFeatures(immune.combined.sct))

clabs<-read.table("cluster_classifications_120724.txt",sep="\t",header=T)
colnames(clabs)[2]<-"ident"
clabs$Group<-ifelse(clabs$Group=="Monocytes / Macrophages","Monocytes /\nMacrophages",clabs$Group) # for plotting purposes

level_order<-levels(immune.combined.sct)
clabs <- clabs[match(level_order, clabs$cluster), ]
new.cluster.ids <- as.character(clabs$ident)
names(new.cluster.ids) <- levels(immune.combined.sct)
immune.combined.sct1 <- RenameIdents(immune.combined.sct, new.cluster.ids)

# Fig 1B UMAP
UMAP1<-UMAPPlot(immune.combined.sct1,cols=clabs$plot_colour,raster=T)+NoLegend()+ggtitle("B")+
  theme(plot.title = element_text(size=30,face = "plain"))
# In order to display the label colours in the correct order, the UMAP1 dataframe itself needs to be sorted according to the identity (cluster name) factor level
UMAP1$data$ident<-factor(UMAP1$data$ident,levels=levels(Idents(immune.combined.sct1)))
UMAP1$data<-UMAP1$data[order(UMAP1$data$ident),]
UMAP2<-LabelClusters(UMAP1, id = "ident", size = 4, repel = T,  box.padding = 1.5, box=T,color=clabs$label_text_colour, segment.color="black",max.overlaps=15)

# Fig. 1C cell counts in clusters
clustcounts<-read.table("cluster_individual_counts_subclusters_040224.txt",header=T,sep="\t")[c(1,5)] %>% distinct()
clustcounts<-merge(clustcounts,clabs,by="cluster")
clustcounts$Shorthand <- factor(clustcounts$Shorthand,levels=clustcounts$Shorthand)
clustcounts$ident <- factor(clustcounts$ident,levels=clustcounts$ident)
clustcounts$Group = factor(clustcounts$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))

cplot<-ggplot(clustcounts,aes(x=Shorthand,y=cluster_count,fill=ident))+
  theme_classic()+
  theme(panel.grid.major.y=element_line(color="darkgrey"))+
  geom_bar(stat="identity")+
  labs(x="",y="Cell count",title="C")+
  scale_fill_manual(values=clustcounts$plot_colour)+
  theme(plot.title = element_text(hjust = 0))+
  scale_y_continuous(n.breaks=10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(size=30,face = "plain"))+
  facet_grid(.~Group,scales="free_x",space="free")+NoLegend()

# blank space for A (exp diagram)
blankplot<-ggplot()+theme_minimal()+labs(title="A")+theme(plot.title = element_text(size=30,face = "plain"))

lay <- rbind(c(1,2),
             c(3,2))

grobz <- lapply(list(blankplot,UMAP2,cplot), ggplotGrob)

pdf("Figure1_FULL_170724.pdf", width = 20, height = 10)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()

# FIG 2

# Neutrophil markers
nlevels(as.factor(subset(fshmk_small,primary_ct=="Neutrophils")$gene_id)) # 9 markers present
immune.combined.sct1<-AddModuleScore(object = immune.combined.sct1,
                                                 features = list(subset(fshmk_small,primary_ct=="Neutrophils")$gene_id),
                                                 ctrl = 10, name = 'Neutrophils')
Nplot<-VlnPlot(immune.combined.sct1,features="Neutrophils1")
Nplotdat<-merge(Nplot$data,clabs,by="ident")
colnames(Nplotdat)[c(2)]<-"expression"
Nplotdat$marker_group<-"Neutrophil markers\n(N=9)"

# Macrophage markers
nlevels(as.factor(subset(fshmk_small,primary_ct=="Macrophages"|primary_ct=="Monocytes")$gene_id)) # 7 markers present
immune.combined.sct1<-AddModuleScore(object = immune.combined.sct1,
                                                 features = list(subset(fshmk_small,primary_ct=="Macrophages"|primary_ct=="Monocytes")$gene_id),
                                                 ctrl = 10, name = 'Macrophages')
Mplot<-VlnPlot(immune.combined.sct1,features="Macrophages1")
Mplotdat<-merge(Mplot$data,clabs,by="ident")
colnames(Mplotdat)[c(2)]<-"expression"
Mplotdat$marker_group<-"Monocyte / Macrophage\nmarkers (N=7)"

#  T-cell markers
nlevels(as.factor(subset(fshmk_small,primary_ct=="T-cells")$gene_id)) # 22 markers present
immune.combined.sct1<-AddModuleScore(object = immune.combined.sct1,
                                                 features = list(subset(fshmk_small,primary_ct=="T-cells")$gene_id),
                                                 ctrl = 10, name = 'T_cells')
Tplot<-VlnPlot(immune.combined.sct1,features="T_cells1")
Tplotdat<-merge(Tplot$data,clabs,by="ident")
colnames(Tplotdat)[c(2)]<-"expression"
Tplotdat$marker_group<-"T-cell markers\n(N=22)"

# B-cell markers
nlevels(as.factor(subset(fshmk_small,primary_ct=="B-cells")$gene_id)) # 12 markers present
immune.combined.sct1<-AddModuleScore(object = immune.combined.sct1,
                                                 features = list(subset(fshmk_small,primary_ct=="B-cells")$gene_id),
                                                 ctrl = 10, name = 'B_cells')
Bplot<-VlnPlot(immune.combined.sct1,features="B_cells1")
Bplotdat<-merge(Bplot$data,clabs,by="ident")
colnames(Bplotdat)[c(2)]<-"expression"
Bplotdat$marker_group<-"B-cell markers\n(N=12)"

allplotdat<-rbind(Nplotdat,Mplotdat,Tplotdat,Bplotdat)
allplotdat$Group = factor(allplotdat$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))
allplotdat$marker_group = factor(allplotdat$marker_group, levels=c("Neutrophil markers\n(N=9)","Monocyte / Macrophage\nmarkers (N=7)","T-cell markers\n(N=22)","B-cell markers\n(N=12)"))

allplot<-ggplot(allplotdat,aes(x=Shorthand,y=expression,fill=ident))+
  theme_classic()+
  geom_point(position=position_jitter(0.25),size=0.05,alpha=0.2)+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  geom_hline(yintercept=0,lty="dashed")+NoLegend()+
  labs(x="",y="Module expression score")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(marker_group~Group,scales="free_x",space="free")+
  geom_hline(aes(yintercept=-Inf)) + 
  coord_cartesian(clip="off")

pdf("Figure2_FULL_170724.pdf", width = 10, height = 8)
allplot
dev.off()

# output a table of the prior markers just of the main cell types
fshmk_smaller<-subset(fshmk_small,primary_ct=="T-cells"|
                        primary_ct=="B-cells"|
                        primary_ct=="Neutrophils"|
                        primary_ct=="Macrophages"|
                        primary_ct=="Monocytes"|
                        primary_ct=="T-cells")
# add ensembl gene names to this as well?
ensembl_trutta<-useEnsembl(biomart="ensembl",dataset="strutta_gene_ensembl",mirror="useast")
ensgnames<-as.data.frame(getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                                 filters = 'ensembl_gene_id',
                                 values = fshmk_smaller$gene_id,
                                 mart = ensembl_trutta))
colnames(fshmk_smaller)[3]<-"ensembl_gene_id"
fshmk_smaller<-merge(fshmk_smaller,ensgnames,by="ensembl_gene_id")
# finalise formatting and write out table
fshmk_smaller$primary_ct<-factor(fshmk_smaller$primary_ct,levels=c("Neutrophils","Monocytes","Macrophages","T-cells","B-cells"))
fshmk_smaller<-fshmk_smaller[order(fshmk_smaller$primary_ct),]
write.table(fshmk_smaller[c(2,3,5,1)],file="immune_markers_in_varft.txt",sep="\t",col.names = T,row.names = F,quote=F)
