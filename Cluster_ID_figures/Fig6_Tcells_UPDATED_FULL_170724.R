# FIGURE 6 T-CELLS UPDATED 170724
#getwd()
library(dplyr)
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
`%not_in%`<-Negate(`%in%`)

# Load in dataset
load("../trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.rds")
DefaultAssay(immune.combined.sct)<-"integrated"

# read in clustermarker data
clustmarkers<-read.delim("all_clustermarkers_annotated_180224.txt",header=T,sep="\t")%>%subset(avg_log2FC>=1)
fshmk_small<-read.table("all_potential_markers_long_290124.txt",sep="\t",header=T)
fshmk_small<-subset(fshmk_small, gene_id%in%VariableFeatures(immune.combined.sct))

# read in cluster labels
clabs<-read.table("cluster_classifications_120724.txt",sep="\t",header=T)
clustmarkers<-merge(clustmarkers,clabs[c(1:3)],by="cluster")

# keep one version with all clustermarkers
clustmarkers_all<-clustmarkers
# but mostly look at high abundance markers
clustmarkers<-subset(clustmarkers,high_abundance==1)

nrow(clustmarkers[!duplicated(clustmarkers$gene_id),])
nrow(clustmarkers_all[!duplicated(clustmarkers_all$gene_id),])


#####################################################
# PRIOR MARKERS THAT ALSO APPEAR AS CLUSTER MARKERS #
#####################################################

subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="T-cells")$gene_id&Shorthand!="T7")%>%distinct(gene_id,name_or_symbol,description)
# 13 --> makes no difference whether T7 is included or not

##################
# 6B) UPSET PLOT #
##################

Tcells<-subset(clustmarkers,Shorthand%like%"T"&Shorthand!="TH")
# subset for those NOT markers of other clusters (EVEN THOSE <60%! THIS IS VERY STRINGENT!)
Tcells<-subset(Tcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Tcells$cluster)$gene_id) 
# exclude T7 for the time being --> important as DEGs could be origin-specific
Tcells<-subset(Tcells,Shorthand!="T7")
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
Tcells<-merge(distinct(Tcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Tcells[!duplicated(Tcells$gene_id),])
# 53 genes are markers uniquely of the T-cells

# aggregated table
Tcellsagg<-aggregate(Shorthand~gene_id, Tcells, FUN=toString)
Tcellsagg<-merge(Tcellsagg,Tcells[c(1,8,10,9,14)],by="gene_id")%>%distinct()

# UPSET PLOT CODE:
Tcells_upset<-as.data.frame(Tcells[c(13,1)] %>% mutate(value = 1) %>%
                              pivot_wider(
                                names_from = Shorthand,
                                values_from = value, values_fill = 0
                              ))
clusters<-levels(as.factor(Tcells$Shorthand))
Tupset<-as.ggplot(upset(Tcells_upset, clusters))

F6B<-as.ggplot(grid.arrange(Tupset,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))

# do same but for all markers
Tcells2<-subset(clustmarkers_all,Shorthand%like%"T"&Shorthand!="TH")
Tcell2s<-subset(Tcells2,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Tcells2$cluster)$gene_id) 
Tcells2<-subset(Tcells2,Shorthand!="T7")
Tcells2<-merge(distinct(Tcells2[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Tcells2[!duplicated(Tcells2$gene_id),])
Tcells2agg<-aggregate(Shorthand~gene_id, Tcells2, FUN=toString)
Tcells2agg<-merge(Tcells2agg,Tcells2[c(1,8,10,9,14)],by="gene_id")%>%distinct()
colnames(Tcells2agg)[2]<-"Shorthand_nonHA"
# the following makes an aggregated table of high abundance markers, but lists also clusters for which those genes are markers regardless of abundance
Tcellsagg<-merge(Tcellsagg,Tcells2agg[c(1,2)],by="gene_id",all.x=T,all.y=F)

colnames(Tcellsagg)[c(3,7)]<-c("clusters_HA","clusters_all")
write.table(Tcellsagg,file="Tcell_markers.tsv",sep="\t",row.names = F,col.names = T,quote=F)

# See the table Tcellsagg to see which markers are pan-T and which mark specific clusters

##################
# 6C) PAN-T-CELL #
##################

# two genes are shared between six:
# CD3E and TRBC1 TCR beta constant 1

p1<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000020148",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD3E-like*",subtitle="ENSSTUG00000020148 | T1-6")
p2<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000037152",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="TRBC1-like*",subtitle="ENSSTUG00000037152 | T1-6")

# markers of five; si:dkey-81j8.6 (STK10-like), def6b, il7r, sh2d1ab, fynb
# I point out those that are markers of a 6th cluster but <60%

# STK10-like
# "negatively regulate interleukin 2 expression in T-cells" (Wiki)
p3<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000008344",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="STK10-like",subtitle="ENSSTUG00000008344 | T1,T2,T4-6\nT3 (<60%)")

# def6b: human homolog has role in T-cells (https://www.nature.com/articles/s41467-019-10812-x)
p4<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000033778",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="def6b",subtitle="ENSSTUG00000033778 | T1,T2,T4-6\nT3 (<60%)")

# il7r (interleukin 7 receptor)
p5<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000036602",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="il7r",subtitle="ENSSTUG00000036602 | T1,T3-6\nT2 (<60%)")

# sh2d1ab: involved in SLAM signalling (in turn involved in T-cell proliferation) (Wiki)
p6<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000037504",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="sh2d1ab",subtitle="ENSSTUG00000037504 | T1-4,T6")

# PLOT
F6C<-as.ggplot(grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# honorable mention: fynb
# fynb: fyn is associated with T-cell signalling (Wiki)
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000037621",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   #theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="fynb",subtitle="T1,T2,T4-6\nT3 (<60%)")

################################
# 6D: CLUSTER-SPECIFIC MARKERS #
################################

# See again Tcellsagg: markers that are specific to certain clusters both with and without high abundance filtering were prioritised

# cd4-2b: T4 high abundance, two others low abundance
p7<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000002855",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cd4-2b*",subtitle="ENSSTUG00000002855 | T4\nT3,T4,T6 (<60%)")

# CD28: T4 high abundance, T6 low abundance
p8<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000013913",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cd28*",subtitle="ENSSTUG00000013913 | T4\nT6 (<60%)")

# another cd28-like
p9<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000006325",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD28-like*",subtitle="ENSSTUG00000006325 | T1,T4,T6\nT2,T3 (<60%)")

# CD5-like
p10<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000043767",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD5-like*",subtitle="ENSSTUG00000043767 | T1,T4,T6\nT3,T5 (<60%)")

# CD8A: T1 high abundance, two others low abundance
p11<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000018765",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cd8a*",subtitle="ENSSTUG00000018765 | T1\nT3,T6 (<60%)")
# MENTION ALSO CD8B (same expression pattern)

# CD2-like: T5 high abundance, two others low abundance
p12<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000007876",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD2-like*",subtitle="ENSSTUG00000007876 | T5\nT1,T6 (<60%)")

# ilr2b: T2 and T5 high abundance
p13<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000015644",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="il2rb",subtitle="ENSSTUG00000015644 | T2,T5")
# found on memory, activated, and regulatory T-cells
# ("promotes the differentiation of T cells into effector T cells and into memory T cells")

# SLAM5-like: T5  only
p14<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000014088",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="SLAM5-like",subtitle="ENSSTUG00000014088 | T5")
# AKA CD84; has been shown to interact with sh2d1a (Wiki)

# plk3 is high abundance exclusively in T2, but is expressed at low abundance in all other clusters
# not so good
p15<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000015614",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="plk3",subtitle="ENSSTUG00000015614 | T2\nT1,T3-6 (<60%)")

# PLOT
F6D<-as.ggplot(grid.arrange(p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=3,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# Also of interest, but not plotted:

# # unnamed CCL3-like: T5 high abundance, T2 low abundance
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000005355",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   #theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="unnamed CCL3-like",subtitle="ENSSTUG00000005355 | T5\nT2 (<60%)")

# # unnamed IGSF3-like: T5 only
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000006882",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   #theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="unnamed IGSF3-like",subtitle="ENSSTUG00000006882 | T5")

# # dnase1l4.1: T5 only
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000017594",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   #theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="dnase1l4.1",subtitle="T5")

# # TCF7: T-cell marker expressed at high abundance in T3 but also others
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000002790",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   #theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="TCF7-like",subtitle="T1,T3,T4,T6\nT2,T5 (<60%)")

###################
# 6A SUMMARY UMAP #
###################

## SUMMARY UMAP FIGURE
uplot<-UMAPPlot(immune.combined.sct)
uplotdat<-uplot$data
colnames(uplotdat)[3]<-"cluster"
# need to manually alter some numbers for compatibility
clabs$cluster<-ifelse(clabs$cluster=="8","8.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="10","10.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="19","19.0",clabs$cluster)
uplotdat<-merge(uplotdat,clabs,by="cluster")
uplotdat<-subset(uplotdat,UMAP_1<12&UMAP_1> -1&UMAP_2>4.5&UMAP_2<12)
uplotdat$Shorthand<-ifelse(uplotdat$Name%like%"T",uplotdat$Shorthand,"Other")
levels(as.factor(uplotdat$Shorthand))

# had to set these manually :/
Tcolors<-c("grey","blue","cornflowerblue","darkslateblue","deepskyblue3","royalblue","lightskyblue","steelblue4")

Tumap<-ggplot(uplotdat,aes(x=UMAP_1,y=UMAP_2))+
  theme_classic()+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  geom_point(aes(color=Shorthand),size=0.1)+
  scale_color_manual(values=Tcolors)+
  NoLegend()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="")+
  # ANNOTATIONS
  annotate(geom="segment",x=4,xend=4,y=11.5,yend=10.6)+
  annotate(geom="richtext",x = 4, y = 11.5, label = "**T4**<br><i>cd4-2b</i><sup>+</sup>, cd28<sup>+</sup><br>mRNA binding activity")+
  annotate(geom="segment",x=9,xend=6.8,y=8,yend=8)+
  annotate(geom="richtext",x = 9, y = 8, label = "**T1**<br><i>cd8a/b</i><sup>+</sup><br>Kinase activity")+
  annotate(geom="segment",x=2,xend=2.75,y=6,yend=7)+
  annotate(geom="richtext",x = 2, y = 6, label = "**T6**<br>Proliferating<br>(<i>mki67</i><sup>+</sup>, <i>pcna</i><sup>+</sup>),<br>Chemotaxis")+
  annotate(geom="segment",x=1.25,xend=1.25,y=7.25,yend=8.25)+
  annotate(geom="richtext",x = 1, y = 7.25, label = "**T5**<br>CD2-like<sup>+</sup>,SLAM5-like<sup>+</sup>,<br><i>il2rb</i><sup>+</sup><br>mRNA binding activity")+
  annotate(geom="segment",x=6,xend=5.25,y=10,yend=9.5)+
  annotate(geom="richtext",x = 6, y = 10, label = "**T3**<br>CD3E-like<sup>+</sup><br>Mitochondrial activity")+
  annotate(geom="segment",x=-0.25,xend=-0.25,y=10.5,yend=8.25)+
  annotate(geom="richtext",x = 0, y = 10.5, label = "**T7**<br>Over-rep. 'Wild' group<br>Chemokine activity")+
  annotate(geom="segment",x=6.5,xend=5,y=6.5,yend=6.5)+
  annotate(geom="richtext",x = 6.5, y = 6.5, label = "**T2**<br><i>il2rb</i><sup>+</sup>,<i>plk3</i><sup>+</sup>")

F6A<-as.ggplot(grid.arrange(Tumap,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))

lay <- rbind(c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(3,3,4,4,4))

grobz <- lapply(list(F6A,F6B,F6C,F6D), ggplotGrob)

pdf("Figure6_UPDATED_FULL_170724.pdf", width = 19, height = 19)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()

# DONE