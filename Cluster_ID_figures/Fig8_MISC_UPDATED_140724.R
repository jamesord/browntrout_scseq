# FIGURE 8 MISCELLANEOUS
# RBCs, thrombocytes, NK-cells, and others
# UPDATED 14/07/24

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
`%not_in%`<-Negate(`%in%`)

# Load in dataset
load("C:/Users/James/Documents/R/trout_local/trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.downsamp.rds")
DefaultAssay(immune.combined.sct.downsampled)<-"integrated"

# read in clustermarker data
clustmarkers<-read.delim("Seurat_cluster_ID/current/all_clustermarkers_annotated_180224.txt",header=T,sep="\t")%>%subset(avg_log2FC>=1)
fshmk_small<-read.table("cell_marker_compilation/fish_compilation/all_potential_markers_long_290124.txt",sep="\t",header=T)
fshmk_small<-subset(fshmk_small, gene_id%in%VariableFeatures(immune.combined.sct.downsampled))

# read in cluster labels
clabs<-read.table("Seurat_cluster_ID/current/cluster_classifications_120724.txt",sep="\t",header=T)
clustmarkers<-merge(clustmarkers,clabs[c(1:3)],by="cluster")

# keep one version with all clustermarkers
clustmarkers_all<-clustmarkers
# but mostly look at high abundance markers
clustmarkers<-subset(clustmarkers,high_abundance==1)

nrow(clustmarkers[!duplicated(clustmarkers$gene_id),])
nrow(clustmarkers_all[!duplicated(clustmarkers_all$gene_id),])

subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="NK_cells")$gene_id)%>%distinct(gene_id,name_or_symbol,description)
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Erythroid_cells")$gene_id)%>%distinct(gene_id,name_or_symbol,description)
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Thrombocytes")$gene_id)%>%distinct(gene_id,name_or_symbol,description)

NKcells<-subset(clustmarkers,Shorthand%like%"NK")
# subset for those NOT markers of other clusters (EVEN THOSE NOT HIGH ABUNDANCE! THIS IS VERY STRINGENT!)
NKcells<-subset(NKcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%NKcells$cluster)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
NKcells<-merge(distinct(NKcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(NKcells[!duplicated(NKcells$gene_id),])
# 4 genes are markers uniquely of the NK-cells (nice)

# perforin (prior marker)
p1<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000019827",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="prf1*",subtitle="ENSSTUG00000019827 | NK")
# NK-lysin-like
p2<-FeaturePlot(immune.combined.sct.downsampled,features="mikado.33G809",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="NK-lysin-like",subtitle="mikado.33G809 | NK")
# Note that NK-lysin is a type of saposin, so we can call it NK-lysin-like

Thcells<-subset(clustmarkers,Shorthand%like%"TH")
# subset for those NOT markers of other clusters (EVEN THOSE NOT HIGH ABUNDANCE! THIS IS VERY STRINGENT!)
Thcells<-subset(Thcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Thcells$cluster)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
Thcells<-merge(distinct(Thcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Thcells[!duplicated(Thcells$gene_id),])
# 41 genes are markers uniquely of the Thrombocytes

# thrombopoietin receptor-like (prior marker)
# "Thrombopoietin was shown to be the major regulator of megakaryocytopoiesis and platelet formation" (Wiki)
# also important for HSC maintenance: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8292222/
p3<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000005928",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="mpl*",subtitle="ENSSTUG00000005928 | TH")
# thrombospondin-1-like
# role in platelet aggregation: https://pubmed.ncbi.nlm.nih.gov/6501568/
p4<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000013788",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="thbs1b",subtitle="ENSSTUG00000013788 | TH")

RBcells<-subset(clustmarkers,Shorthand%like%"RBC")
# subset for those NOT markers of other clusters (EVEN THOSE NOT HIGH ABUNDANCE! THIS IS VERY STRINGENT!)
RBcells<-subset(RBcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%RBcells$cluster)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
RBcells<-merge(distinct(RBcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(RBcells[!duplicated(RBcells$gene_id),])
# 11 genes are markers uniquely of the RBCs
# important to identify these because Hemoglobin is expressed elsewhere

# canonical RBC marker but NOT unique to RBC cluster:
p5<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000013812",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="hbaa2*",subtitle="ENSSTUG00000013812 | RBC\nNA1 (<60%)")
# one of 11 unique to RBC cluster:
# 5-aminolevulinate synthase, erythroid-specific, mitochondrial-like
p6<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000037820",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="alas2",subtitle="ENSSTUG00000037820 | RBC")

#subset(clustmarkers_all,gene_id=="ENSSTUG00000013812")

# MP --> putative myeloid precursors
MPcells<-subset(clustmarkers,Shorthand%like%"MP")
MPcells<-subset(MPcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%MPcells$cluster)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
MPcells<-merge(distinct(MPcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T)
# 4 unique markers:
# ptprga, kita, egr1, NXPH4
# could not find anything on ptprga relating to myleoid precursors or hematopoiesis / stem cells
# Kit / a.k.a c-kit is expressed highly in HSCs and progenitors: https://pubmed.ncbi.nlm.nih.gov/17350321/
# EGR1 positively regulates myeloid differentiation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9340249/
# Neurexophilin: expressed in neural cells but also implicated in suppression of proliferation in HSCs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3142901/
# see also the list prior to filtering for markers of other clusters:
# includes proliferation marker PCNA, H2AZ, thrombocyte marker mpl, MHCII B (antigen), CD74 (monocyte)
# also three 'fos' genes and JUND: these along with egr1 are all implicated in early myeloid cell function: https://www.sciencedirect.com/science/article/pii/S0006497120519794
# so I think the evidence is there that these are myeloid progenitors
p7<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000010274",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="kita",subtitle="ENSSTUG00000010274 | MP")
p8<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000012364",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="egr1",subtitle="ENSSTUG00000012364 | MP")

# RN1 OR RN2 --> putative renal cells
RNcells<-subset(clustmarkers,Shorthand%like%"RN1"|Shorthand%like%"RN2")
# subset for those NOT markers of other clusters (EVEN THOSE NOT HIGH ABUNDANCE! THIS IS VERY STRINGENT!)
RNcells<-subset(RNcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%RNcells$cluster)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
RNcells<-merge(distinct(RNcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(RNcells[!duplicated(RNcells$gene_id),])
# 17 genes are markers uniquely of RN1/2 (nice); note GAPDH is there

# collectrin: most abundantly expressed in kidneys https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4164340/
# GST: "in humans, renal proximal tubular cells contain high concentrations of alpha GST", i.e. kidneys (Wiki)
# CD59: kidney endothelial cells; https://pubmed.ncbi.nlm.nih.gov/9555661/
# So basically these are cells of the kidney

#subset(clustmarkers_all,gene_id=="ENSSTUG00000049645")

# CD59-like
p9<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000040811",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="ly97.3",subtitle="ENSSTUG00000012672 | RN1\nRN2 (<60%)")
# NOTE these are the only clusters to express GAPDH
p10<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000049645",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="gapdh",subtitle="ENSSTUG00000049645 | RN1\nRN2 (<60%)")
# NOTE same expression pattern for collectrin (ENSSTUG00000012672) and GST (ENSSTUG00000037965)

# NA2 strange patterns: for supplementary if anything
# NA2cells<-subset(clustmarkers,Shorthand%like%"NA2")
# # NA2 has one high abundance marker, collagenase 3
# # subset for those NOT markers of other clusters (EVEN THOSE NOT HIGH ABUNDANCE! THIS IS VERY STRINGENT!)
# NA2cells<-subset(NA2cells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%NA2cells$cluster)$gene_id) 
# # However it is not unique to NA2, so it has no unique markers
# FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000018003",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="MMP13-like",subtitle="ENSSTUG00000018003 | NA2,N1")
# # Note also however the elevated expression of B-cell markers, with IGLC1-like being a prime example:
# FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000015148",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="IGLC1-like",subtitle="ENSSTUG00000015148 | B1,B2,B5,B6\nNA2(logFC<1)")

# NA3
NA3cells<-subset(clustmarkers,Shorthand%like%"NA3")
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
NA3cells<-merge(distinct(NA3cells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(NA3cells[!duplicated(NA3cells$gene_id),])
# only one of these is unique to NA3: NFIL3
# Mention some observations but maybe no need to plot:
# has high abundance markers of:
# NK cells (PRF1-like, shared with cluster NK)
# T-cells (CD2-like, shared expression with T-cell clusters)
# Monocytes (CD74-like, shared expression with macrophages and B-cells)
# see also TOX2: shared expression with T-cell clusters

# so we would just make 10 plots
F8A<-as.ggplot(grid.arrange(p1,p2,nrow=2,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
F8B<-as.ggplot(grid.arrange(p3,p4,nrow=2,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
F8C<-as.ggplot(grid.arrange(p5,p6,nrow=2,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
F8D<-as.ggplot(grid.arrange(p7,p8,nrow=2,top=textGrob("D", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
F8E<-as.ggplot(grid.arrange(p9,p10,nrow=2,top=textGrob("E", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))

pdf("Figure8_UPDATED_140724.pdf", width = 17, height = 7)
grid.arrange(F8A,F8B,F8C,F8D,F8E,nrow=1)
dev.off()