# SUPP FIGURE: MT EXPRESSION, N GENES, TOTAL GENE EXPRESSION
library(dplyr)
library(ggplot2)
library(data.table)
library(Seurat)
library(tidyr)
library(gridExtra)
library(forcats)

load("../trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.rds")
DefaultAssay(immune.combined.sct)<-"integrated"

levels(Idents(immune.combined.sct))<-gsub("8.0","8",levels(Idents(immune.combined.sct)))
levels(Idents(immune.combined.sct))<-gsub("10.0","10",levels(Idents(immune.combined.sct)))
levels(Idents(immune.combined.sct))<-gsub("19.0","19",levels(Idents(immune.combined.sct)))

clabs<-read.table("cluster_classifications_120724.txt",sep="\t",header=T)
colnames(clabs)[2]<-"ident"
clabs$Group<-ifelse(clabs$Group=="Monocytes / Macrophages","Monocytes /\nMacrophages",clabs$Group) # for plotting purposes

level_order<-levels(immune.combined.sct)
clabs <- clabs[match(level_order, clabs$cluster), ]
new.cluster.ids <- as.character(clabs$ident)
names(new.cluster.ids) <- levels(immune.combined.sct)
immune.combined.sct1 <- RenameIdents(immune.combined.sct, new.cluster.ids)

# MT markers
mtgenes<-read.table("mitochondiral_gene_ENSEMBL_IDs.txt")
# add module for MT genes
immune.combined.sct<-AddModuleScore(object = immune.combined.sct1,
                                                features = list(mtgenes$V1),
                                                ctrl = 10, name = 'MT')
# Warning: The following features are not present in the object: ENSSTUG00000000007, ENSSTUG00000000022,
# ENSSTUG00000000026, ENSSTUG00000000028, ENSSTUG00000000029, ENSSTUG00000000033, ENSSTUG00000000034,
# not searching for symbol synonyms

# generate base plots
MTplot<-VlnPlot(immune.combined.sct,features="MT1")
nftplot<-VlnPlot(immune.combined.sct,features="nFeature_RNA")
nctplot<-VlnPlot(immune.combined.sct,features="nCount_RNA")

# nice plot of MT gene module
MTplotdat<-merge(MTplot$data,clabs,by="ident")
colnames(MTplotdat)[c(2)]<-"expression"
MTplotdat$Group = factor(MTplotdat$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))
MTplot<-ggplot(MTplotdat,aes(x=Shorthand,y=expression,fill=ident))+
  theme_classic()+
  geom_point(position=position_jitter(0.25),size=0.1,alpha=0.2)+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  geom_hline(yintercept=0,lty="dashed")+NoLegend()+
  labs(x="",y="Module expression score",title="A",subtitle="Mitochondrial genes (6)")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~Group,scales="free_x",space="free")

# nice plot of number of expressed genes
nftplotdat<-merge(nftplot$data,clabs,by="ident")
colnames(nftplotdat)[c(2)]<-"n_feature_RNA"
nftplotdat$Group = factor(nftplotdat$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))
nftplot<-ggplot(nftplotdat,aes(x=Shorthand,y=n_feature_RNA,fill=ident))+
  theme_classic()+
  geom_point(position=position_jitter(0.25),size=0.1,alpha=0.2)+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  #geom_hline(yintercept=0,lty="dashed")+
  NoLegend()+
  labs(x="",y="N expressed genes",title="B",subtitle="")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~Group,scales="free_x",space="free")

# nice plot of gene expression levels
nctplotdat<-merge(nctplot$data,clabs,by="ident")
colnames(nctplotdat)[c(2)]<-"n_count_RNA"
nctplotdat$Group = factor(nctplotdat$Group, levels=c("Neutrophils","Monocytes /\nMacrophages","T-cells","B-cells","Other","Unclassified"))
nctplot<-ggplot(nctplotdat,aes(x=Shorthand,y=log(n_count_RNA),fill=ident))+
  theme_classic()+
  geom_point(position=position_jitter(0.25),size=0.1,alpha=0.2)+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  #geom_hline(yintercept=0,lty="dashed")+
  NoLegend()+
  labs(x="",y="Log2(read count of expressed genes)",title="C",subtitle="")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~Group,scales="free_x",space="free")

pdf("FigureS_MT_FULL_191124.pdf", width = 10, height = 9)
grid.arrange(MTplot,nftplot,nctplot,ncol=1)
dev.off()
