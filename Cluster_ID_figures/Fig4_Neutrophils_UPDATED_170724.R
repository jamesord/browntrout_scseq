# NEUTROPHIL FIGURE UPDATED 170724

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

#####################################################
# PRIOR MARKERS THAT ALSO APPEAR AS CLUSTER MARKERS #
#####################################################

subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Neutrophils")$gene_id)%>%distinct(gene_id)
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Granulocytes")$gene_id)%>%distinct(gene_id)
# There are five: mpx, LYZ, mmp9, HCE and CPA5 (named CPA1 in ensembl)
# ADDITIONALLY: NCF1, which we have as a GRANULOCYTE maker
# The two MMP9 homologs have similar expression patterns so in the figure we can show just one of them

##################
# 4B) UPSET PLOT #
##################

Ncells<-subset(clustmarkers,cluster==1|cluster==2|cluster==3|cluster==6|cluster==20|cluster==27)
# subset for those NOT markers of other clusters # NOTE that we exclude 23 from the exclusion (if that makes sense)
# also note that we exclude genes found as markers as other clusters in 'clustermarkers_all' (including non-high abundance markers)
Ncells<-subset(Ncells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Ncells$cluster&cluster!=23)$gene_id)
# some rows are duplicated as genes can be markers of either lineages or functions (e.g. antigen presenting cells):
Ncells<-merge(distinct(Ncells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Ncells[!duplicated(Ncells$gene_id),])
# 45 genes are markers of the 'neutrophil' clusters

# aggregated table
Ncellsagg<-aggregate(Shorthand~gene_id, Ncells, FUN=toString)
Ncellsagg<-merge(Ncellsagg,Ncells[c(1,8,10,9,14)],by="gene_id")%>%distinct()

# IN CASE WE WANT TO LOOK ALSO AT LOW ABUNDANCE MARKERS: do same but including markers that are not high abundance
Ncells2<-subset(clustmarkers_all,cluster==1|cluster==2|cluster==3|cluster==6|cluster==20|cluster==27)
Ncells2<-subset(Ncells2,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Ncells2$cluster&cluster!=23)$gene_id) 
Ncells2<-merge(distinct(Ncells2[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Ncells2[!duplicated(Ncells2$gene_id),])
# 70 genes when including low abundance ones

# UPSET PLOT CODE:
Ncells_upset<-as.data.frame(Ncells[c(13,1)] %>% mutate(value = 1) %>%
                              pivot_wider(
                                names_from = Shorthand,
                                values_from = value, values_fill = 0
                              ))
clusters<-c("N1","N2","N3","N4","N5","N6")
Nupset<-as.ggplot(upset(Ncells_upset, clusters))
F4B<-as.ggplot(grid.arrange(Nupset,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))

# do same but for all markers
# IN CASE WE WANT TO LOOK ALSO AT LOW ABUNDANCE MARKERS: do same but including markers that are not high abundance
Ncells2<-subset(clustmarkers_all,cluster==1|cluster==2|cluster==3|cluster==6|cluster==20|cluster==27)
Ncells2<-subset(Ncells2,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Ncells2$cluster&cluster!=23)$gene_id)
Ncells2<-merge(distinct(Ncells2[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Ncells2[!duplicated(Ncells2$gene_id),])
Ncells2agg<-aggregate(Shorthand~gene_id, Ncells2, FUN=toString)
Ncells2agg<-merge(Ncells2agg,Ncells2[c(1,8,10,9,14)],by="gene_id")%>%distinct()
colnames(Ncells2agg)[2]<-"Shorthand_nonHA"
# the following makes an aggregated table of high abundance markers, but lists also clusters for which those genes are markers regardless of abundance
Ncellsagg<-merge(Ncellsagg,Ncells2agg[c(1,2)],by="gene_id",all.x=T,all.y=F)

# See the table Ncellsagg to see which markers are pan-neutrophil and which mark specific clusters

colnames(Ncellsagg)[c(2,7)]<-c("clusters_HA","clusters_all")
write.table(Ncellsagg,file="Ncell_markers.tsv",sep="\t",row.names = F,col.names = T,quote=F)

# No genes shared by all clusters

# The patters we then show are:
# - expression specific to upper clusters
# - expression specific to lower cluster
# - expression shared by upper and lower clusters (pan-neutrophil)

# Genes to plot
# C) pan-neutrophil (expression common to upper and lower clusters):
#    - taldo1
#    - mpx*
#    - alox5ap
#    - gpx1b
#    - FBP1-like
#    - tktb
# D) expression specific to upper clusters:
#    - hce*
#    - CFD-like
#    - cpa1*
#    - tmsb2
#    - npdc1a
#    - unknown N4-specific
# E) expression specific to lower cluster:
#    - gadl1
#    - mmp9*
#    - adam28

# NOTE 190324: LYZ is NOT unique to neutrophils as it appears in macrophage clusters; replaced in figure with GADL1
# another figure should be prepared for such markers that appear in multiple cell types

################################################################
# 4C) GENES WITH EXPRESSION SHARED BY UPPER AND LOWER CLUSTERS #
################################################################

# Four genes are markers of N1,N2,N3 and N5:
distinct(Ncells[c(13,1,8,9,14,3,4,5)]) %>% add_count(gene_id) %>%
  subset(n==4&(Shorthand=="N1"|Shorthand=="N2"|Shorthand=="N3"|Shorthand=="N5"))%>% add_count(gene_id)%>%
  subset(nn==4)%>%distinct(gene_id,name_or_symbol,description)

# gpx1b,taldo1 (x2), and PLAC8-like mark these lower clusters
p10<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000017863",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="gpx1b",subtitle="ENSSTUG00000017863 | N1-3,N5")

# taldo1 is involved with oxidative burst in neutrophils https://www.nature.com/articles/s42255-022-00550-8
p11<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000044076",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="taldo1",subtitle="ENSSTUG00000044076 | N1-3,N5")

# N1,N2, N5 and N6: mpx
distinct(Ncells[c(13,1,8,9,14,3,4,5)]) %>% add_count(gene_id) %>%
  subset(n==4&(Shorthand=="N1"|Shorthand=="N2"|Shorthand=="N5"|Shorthand=="N6"))%>%
  add_count(gene_id)%>% subset(nn==4)%>%distinct(gene_id,name_or_symbol,description)

p12<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000007163",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="mpx*",subtitle="ENSSTUG00000007163 | N1,N2,N5,N6",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# or of just 1 and 2 (marking the largest compartments):
cN12<-distinct(Ncells[c(13,1,8,9,14,3,4,5)]) %>% add_count(gene_id) %>%
  subset(n==2&(Shorthand=="N1"|Shorthand=="N2"))%>% add_count(gene_id) %>%
  subset(nn==2)#%>%distinct(gene_id,name_or_symbol,description)
# Focusing on genes that were shared across the largest two neutrophil clusters
# Note that *all* of these genes have stronger expression in N1 and most show a pattern of increasing expression towards cluster N1
mean(subset(cN12,Shorthand=="N1")$avg_log2FC)-mean(subset(cN12,Shorthand=="N2")$avg_log2FC)
# we can say that the mean was 0.76 higher for N1 than N2

# transkelotase; see this paper also referring to taldo1 https://www.nature.com/articles/s42255-022-00550-8
# expression throughout N2
p13<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000007986",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.y = element_blank())+theme(axis.line.x = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="tktb",subtitle="ENSSTUG00000007986 | N1,N2")

# FBP (fructose-1,6-bisphosphatase) is expressed in cultured neutrophils upon LPS exposure, and also plays a role in oxidative burst:
# https://www.cell.com/cell-metabolism/pdfExtended/S1550-4131(20)30651-3
# https://www.nature.com/articles/s42255-022-00550-8
# expression throughout N2
p14<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000010287",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="FBP1-like",subtitle="ENSSTUG00000010287 | N1,N2")
# alox5ap = expression throughout N2
p15<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000012256",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="alox5ap",subtitle="ENSSTUG00000012256 | N1,N2")

F4C<-as.ggplot(grid.arrange(p10,p11,p12,p13,p14,p15,nrow=3,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# Interestingly, expression is often lacking in the right side of the mass, suggesting either:
# - a continuation of the progenitor lineage
# - active vs inactive states

##############################################################
# 4D) GENES WITH EXPRESSION SPECIFIC TO THE UPPER CLUSTER(S) #
##############################################################

# ONLY HCE SHARED BY ALL UPPER CLUSTERS INCLUDING N6!
p1<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000011854",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="HCE1-like*",subtitle="ENSSTUG00000011854 | N2-6",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank(),axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Check also the one shared by N2-4
# tmsb1,CFD-like (x2),CPA1,

p2<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000045318",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="CFD-like",subtitle="ENSSTUG00000045318 | N2-5",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
# Note that this is one of two homologs marking the same cluster.
# Note however that the other homolog ENSSTUG00000009474 shows slightly stronger expression in N3
# NOTE complement factor D has been shown to inhibit degranulation:
# https://pubmed.ncbi.nlm.nih.gov/7556615/
p3<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000021190",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="cpa1*",subtitle="ENSSTUG00000021190 | N2-5",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# beta thymosin 1
p4<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000006838",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="tmsb1",subtitle="ENSSTUG00000006838 | N2-5",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
# see review on beta thymosin 1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7747025/

# additional
# other markers specific to upper clusters

# N2: npdc1a --> a gene associated with neural tissue
p5<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000039310",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="npdc1a",subtitle="ENSSTUG00000039310 | N2\nN1,N2 (<60%)",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# N2: kcnq5a potassium voltage-gated channel subfamily KQT member 5-like
# FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000025248",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   labs(title="kcnq5a",subtitle="N2",x="",y="")+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())

# N4: unnamed gene, scattered expression
p6<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000030567",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="unnamed IG-like",subtitle="ENSSTUG00000030567 | N4\nN2,N3,N5 (<60%)",x="",y="")+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
# most blast hits to uncharacterised proteins, some to T-cell receptor V-alpha
# also contains IG domain, see:
# https://www.ensembl.org/Salmo_trutta/Transcript/ProteinSummary?db=core;g=ENSSTUG00000030567;r=21:22339998-22418327;t=ENSSTUT00000074026

F4D<-as.ggplot(grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,top=textGrob("D", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# no unique to N5 or N6 even amongst low abundance markers

###########################################################
# 4E) GENES WITH EXPRESSION SPECIFIC TO THE LOWER CLUSTER #
###########################################################

# Cluster 1 is the most distinct with 14 'unique' 'cluster markers'
c1ms<-distinct(Ncells[c(13,1,8,9,14,3,4,5)]) %>% add_count(gene_id) %>% subset(n==1&Shorthand=="N1")
c1ms <- c1ms[order(-c1ms$avg_log2FC),]
# the top ones are enzymes, specifically metalloproteinases.

# We would have mmp9, GADL1, adam28 (another metalloprotease)
p7<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000015832",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  labs(title="mmp9*",subtitle="ENSSTUG00000015832 | N1",x="",y="")+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())
# NOTE that collagenase 3-like (MMP13) has very similar expression
# acidic amino acid decarboxylase GADL1-like
p8<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000019948",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="GADL1",subtitle="ENSSTUG00000019948 | N1")
# adam28: another metalloprotease. In humans this gene is expressed largely in B-cells as well as basophils and eosinophils (see protein atlas)
p9<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000002143",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="adam28",subtitle="ENSSTUG00000002143 | N1")
# a similar gene, ADAM8, is highly expressed in mature granulocytes:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8855804/

F4E<-as.ggplot(grid.arrange(p7,p8,p9,ncol=1,top=textGrob("E", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

#######################
# SUMMARY UMAP FIGURE #
#######################

uplot<-UMAPPlot(immune.combined.sct.downsampled)
uplotdat<-uplot$data
colnames(uplotdat)[3]<-"cluster"
# need to manually alter some numbers for compatibility
clabs$cluster<-ifelse(clabs$cluster=="8","8.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="10","10.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="19","19.0",clabs$cluster)
uplotdat<-merge(uplotdat,clabs,by="cluster")
uplotdat<-subset(uplotdat,UMAP_1<5&UMAP_2<7)
uplotdat$Shorthand<-ifelse(uplotdat$Name%like%"Neutrophil",uplotdat$Shorthand,"Other")
Ncolors<-c(subset(clabs,Name%like%"Neutrophil")$plot_colour,"grey")

Numap<-ggplot(uplotdat,aes(x=UMAP_1,y=UMAP_2))+
  theme_classic()+
  expand_limits(x = -15,y=-13)+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  geom_point(aes(color=Shorthand),size=0.1)+
  scale_color_manual(values=Ncolors)+
  annotate(geom="richtext",x = -8, y = 5.5, label = "**N4**<br>Lacking enzymatic\activity")+
  annotate(geom="richtext",x = 1, y = 0.5, label = "**N3**<br>Proliferating<br>(<i>mki67</i><sup>+</sup>, <i>pcna</i><sup>+</sup>)")+
  annotate(geom="richtext",x = -4, y = -12, label = "**N1**<br><i>lyz</i><sup>+</sup>, <i>mpx</i><sup>+</sup><br>Metallopeptidase activity")+
  annotate(geom="richtext",x = -12, y = -8, label = "**N2**<br>HCE1-like<sup>+</sup>, <i>mpx</i><sup>+</sup><br>Carbohydrate metabolism")+
  annotate(geom="segment",x=-0,xend=-4,y=-2,yend=-2)+
  annotate(geom="richtext",x = 0.5, y = -2.5, label = "**N6**<br>HCE1-like<sup>+</sup>, <i>mpx</i><sup>+</sup><br>Mitochondrial activity")+
  annotate(geom="segment",x=-12,xend=-6,y=1,yend=1)+
  annotate(geom="richtext",x = -12, y = 1, label = "**N5**<br>HCE1-like<sup>+</sup>, <i>mpx</i><sup>+</sup><br>No unique markers")+
  NoLegend()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="")

F4A<-as.ggplot(grid.arrange(Numap,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))

lay <- rbind(c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(3,3,4,4,5),
             c(3,3,4,4,5),
             c(3,3,4,4,5))

grobz <- lapply(list(F4A,F4B,F4C,F4D,F4E), ggplotGrob)

pdf("Figure4_UPDATED_170724.pdf", width = 19, height = 19)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()

# OK good. There is enough here.