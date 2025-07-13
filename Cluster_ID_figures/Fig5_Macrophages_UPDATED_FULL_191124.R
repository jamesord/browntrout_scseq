# Macrophages / Monocytes figure UPDATED 191124

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

#####################################################
# PRIOR MARKERS THAT ALSO APPEAR AS CLUSTER MARKERS #
#####################################################

# Two macrophage markers: MPEG1, irf8
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Macrophages")$gene_id)%>%distinct(gene_id,name_or_symbol,description)
# Two monocyte markers, both CD74 homologs (aka histocompatibility antigen gamma chain)
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="Monocytes")$gene_id)%>%distinct(gene_id,name_or_symbol,description)
# Check also antigen presenting cell markers
subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="APCs")$gene_id)%>%distinct(gene_id,name_or_symbol,description) # all marking MA1,2,3,4

# NOTE THAT CD74 IS A MONOCYTE MARKER IN ZEBRAFISH:
# https://pubmed.ncbi.nlm.nih.gov/29229905/ https://pubmed.ncbi.nlm.nih.gov/35027916/
FeaturePlot(immune.combined.sct,feature="ENSSTUG00000008754",min.cutoff = "q10", max.cutoff = "q90",raster=T)
FeaturePlot(immune.combined.sct,feature="ENSSTUG00000006024",min.cutoff = "q10", max.cutoff = "q90",raster=T)

##################
# 5B) UPSET PLOT #
##################

Mcells<-subset(clustmarkers,Shorthand=="M1"|Shorthand=="M2"|Shorthand=="M3"|Shorthand=="M4"|Shorthand=="M5")
# subset for those NOT markers of other non-macrophage clusters (EVEN THOSE <60%! THIS IS VERY STRINGENT!)
Mcells<-subset(Mcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Mcells$cluster)$gene_id) 
# NOTE that M2 does not have markers that are not appearing in non-macro/mono clusters
# some rows are duplicated as genes can be markers of either macrophage/myeloid or antigen presenting cells:
Mcells<-merge(distinct(Mcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Mcells[!duplicated(Mcells$gene_id),])
# 25 genes are markers uniquely of the monocyte/macrophage clusters

# aggregated table
Mcellsagg<-aggregate(Shorthand~gene_id, Mcells, FUN=toString)
Mcellsagg<-merge(Mcellsagg,Mcells[c(1,8,10,9,14)],by="gene_id")%>%distinct()

# UPSET PLOT CODE:
Mcells_upset<-as.data.frame(Mcells[c(13,1)] %>% mutate(value = 1) %>%
                              pivot_wider(
                                names_from = Shorthand,
                                values_from = value, values_fill = 0
                              ))
clusters<-c("M1","M5","M3","M4")
Mupset<-as.ggplot(upset(Mcells_upset, clusters))

F5B<-as.ggplot(grid.arrange(Mupset,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))

# do same but for all markers
# IN CASE WE WANT TO LOOK ALSO AT LOW ABUNDANCE MARKERS: do same but including markers that are <60%
Mcells2<-subset(clustmarkers_all,Shorthand=="M1"|Shorthand=="M2"|Shorthand=="M3"|Shorthand=="M4"|Shorthand=="M5")
Mcells2<-subset(Mcells2,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Mcells2$cluster)$gene_id) 
Mcells2<-merge(distinct(Mcells2[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Mcells2[!duplicated(Mcells2$gene_id),]) # 73
Mcells2agg<-aggregate(Shorthand~gene_id, Mcells2, FUN=toString)
Mcells2agg<-merge(Mcells2agg,Mcells2[c(1,8,10,9,14)],by="gene_id")%>%distinct()
colnames(Mcells2agg)[2]<-"Shorthand_nonHA"
# the following makes an aggregated table of high abundance markers, but lists also clusters for which those genes are markers regardless of abundance
Mcellsagg<-merge(Mcellsagg,Mcells2agg[c(1,2)],by="gene_id",all.x=T,all.y=F)

# See the table Mcellsagg to see which markers are pan-mono/macro and which mark specific clusters

colnames(Mcellsagg)[c(3,7)]<-c("clusters_HA","clusters_all")
#write.table(Mcellsagg,file="Mcell_markers.tsv",sep="\t",row.names = F,col.names = T,quote=F)

nrow(Mcells[!duplicated(Mcells$gene_id),])
nrow(Mcells2[!duplicated(Mcells$gene_id),])

######################
# 5C: PAN MACROPHAGE #
######################

# No genes were markers of all macrophage+mono clusters
# however, four are shared between M1 and M3:

# c1qa
# https://en.wikipedia.org/wiki/Complement_component_1q
# very highly expressed in macrophages: https://www.proteinatlas.org/ENSG00000159189-C1QC/single+cell+type
# also some in M2
p1<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000022041",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="c1qa",subtitle="ENSSTUG00000022041 | M1,M3\nM2 (<60%)")
# c1qb and c have the same expression profile
# they are also all markers of M2 at 33-41%

# olfactomectin: 38% in M2
# facilitates cell adhesion (wiki) --> also expressed in M2 (MARKER at 38%)
p2<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000007928",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="OLFM4-like",subtitle="ENSSTUG00000007928 | M1,M3\nM2 (<60%)")
# https://ashpublications.org/blood/article/114/22/1356/110498/Olfactomedin-4-Is-Essential-for-Superoxide
# NOTE some isolated but strong expression in the lower cluster(s)

# SEE ALSO HIGH ABUNDANCE MARKERS OF M1 THAT ARE LOW ABUNDANCE MARKERS OF OTHER CLUSTERS
# MPEG1 # note also a marker of M2,3,and 4 but <60% (53-57%)
# is clearly pan-macrophage
p3<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000031848",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="MPEG1*",subtitle="ENSSTUG00000031848 | M1\nM2,M3,M4 (<60%)")
# bactericidal permeability-increasing protein-like --> also in M3 at 36%
p4<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000001658",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="BPI-like",subtitle="ENSSTUG00000001658 | M1\nM3 (<60%)")

# See also:
# # CD22-like 41% in M2
# FeaturePlot(immune.combined.sct,features="mikado.1G2286",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD22-like",subtitle="mikado.1G2286 | M1\nM2 (<60%)")
# 
# Not representative (low abundance in other clusters)
# # keratin, cytoskeletal 18-like
# # https://www.uniprot.org/uniprotkb/P05783/entry
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000044613",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="K1C18-like",subtitle="ENSSTUG00000044613 | M1\nM3,M4 (<60%)")

# PLOT 4 of ABOVE (4/10)

F5C<-as.ggplot(grid.arrange(p1,p2,p3,p4,nrow=2,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# MENTION BUT DON'T PLOT?
# high affinity immunoglobulin gamma Fc receptor I-like: found on surface of macrophages and other immune cells --> M1 and M3
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000008409",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="FCGR1A-like",subtitle="M1\nM3 (<60%)")
# specific to macrophages in mice:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4890338/

# receptor-type tyrosine-protein phosphatase S-like --> expression similar to above CD22-like
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000024661",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="ptprsa",subtitle="M1\nM2 (<60%)")

# GOOD so far
# that now leaves those markers that are unique to particular clusters

#############################################
# 5D: MARKERS UNIQUE TO PARTICULAR CLUSTERS #
#############################################

# CLUSTER M1

# the following were found to be highly specific to M1:
# M1 expresses lots of keratin --> also in M3 and 4 but basically unique to M1
p5<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000044613",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="KRT18-like",subtitle="ENSSTUG00000044613 | M1\nM3,M4 (<60%)")
# macrophage mannose receptor 1 (CD206)
p6<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000033316",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD206-like*",subtitle="ENSSTUG00000033316 | M1")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3963702/

# FCGR1-like (high affinity immunoglobulin gamma Fc receptor)
# Mediates IgG effector functions on monocytes triggering antibody-dependent cellular cytotoxicity (ADCC) of virus-infected cells
# https://www.uniprot.org/uniprotkb/P12314/entry
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000008409",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="FCGR1-like",subtitle="ENSSTUG00000008409 | M1\nM3 (<60%)")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3963702/

# s100a10a --> marker of M1
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000042289",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="s100a10a",subtitle="ENSSTUG00000042289 | M1")
# TLR signalling; constitutively expressed on macrophages:
# https://www.nature.com/articles/s41423-019-0278-1

# PLOT ABOVE (6/10)

# M3
distinct(Mcells[c(13,1,8,9,14)]) %>% add_count(gene_id) %>% subset(n==1&Shorthand=="M3")
# only two
# fatty acid binding protein 7
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000013745",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="FABP7",subtitle="ENSSTUG00000013745 | M3")
# # ependymin-like: # a glycoprotein https://en.wikipedia.org/wiki/Ependymin
# # note it also seems to be in M2 but not a marker
# FeaturePlot(immune.combined.sct,features="ENSSTUG00000013542",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="epdl1",subtitle="ENSSTUG00000013542 | M3")
# NOTE that M3 also expressed MARCO as a marker in low abundance (55%) --> note also some expression in part of M2
p7<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000016758",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="MARCO*",subtitle="ENSSTUG00000013542\nM3 (<60%)")

# PLOT THE ABOVE (7/10)

# M4
distinct(Mcells[c(13,1,8,9,14)]) %>% add_count(gene_id) %>% subset(n==1&Shorthand=="M4")
# just one
p8<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000013713",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="flt3",subtitle="ENSSTUG00000013713 | M4")
# see here: https://en.wikipedia.org/wiki/FMS-like_tyrosine_kinase_3_ligand
# essentially a DENDRITIC CELL marker, is similar to csf-1

# PLOT IT (8/10)

# M5
distinct(Mcells[c(13,1,8,9,14)]) %>% add_count(gene_id) %>% subset(n==1&Shorthand=="M5")
# just one --> complement factor d
p9<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000040640",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cfd",subtitle="ENSSTUG00000040640 | M5")
# M5 low abundance
# mrc1a is macrophage marker (54%)
p10<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000000974",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="mrc1a*",subtitle="ENSSTUG00000000974\nM5 (<60%)")

# PLOT THE ABOVE (10/10)
F5D<-as.ggplot(grid.arrange(p5,p6,p7,p8,p9,p10,nrow=2,top=textGrob("D", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

#######################
# SUMMARY UMAP FIGURE #
#######################

uplot<-UMAPPlot(immune.combined.sct)
uplotdat<-uplot$data
colnames(uplotdat)[3]<-"cluster"
# need to manually alter some numbers for compatibility
clabs$cluster<-ifelse(clabs$cluster=="8","8.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="10","10.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="19","19.0",clabs$cluster)
uplotdat<-merge(uplotdat,clabs,by="cluster")
uplotdat<-subset(uplotdat,UMAP_1<12&UMAP_1> -5&UMAP_2>1)
uplotdat$Shorthand<-ifelse(uplotdat$Group=="Monocytes / Macrophages",uplotdat$Shorthand,"Other")
# had to set these manually :/
Mcolors<-c("darkorange3","orange","orange3","darkorange1","rosybrown","grey")

Mumap<-ggplot(uplotdat,aes(x=UMAP_1,y=UMAP_2))+
  theme_classic()+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  geom_point(aes(color=Shorthand),size=0.1)+
  scale_color_manual(values=Mcolors)+
  NoLegend()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="")+
  # ANNOTATIONS
  annotate(geom="richtext",x = 6, y = 12, label = "**M1**<br>CD206-like<sup>+</sup>, <i>MPEG1</i><sup>+</sup><br>Antigen-presenting")+
  annotate(geom="segment",x=3,xend=3,y=7,yend=4.5)+
  annotate(geom="richtext",x = 4, y = 7, label = "**M3**<br><i>FABP7</i><sup>+</sup><br><i>MARCO</i><sup>+</sup> cell-containing<br>Mitochondrial activity<br>Antigen-presenting")+
  annotate(geom="segment",x=0,xend=1.5,y=2.5,yend=3)+
  annotate(geom="richtext",x = -2, y = 2, label = "**M2**<br><i>MPEG1</i><sup>+</sup> cell-containing<br>Mitochondrial activity<br>Antigen-presenting")+
  annotate(geom="segment",x=6,xend=4.1,y=2.1,yend=2.1)+
  annotate(geom="richtext",x = 7.5, y = 2, label = "**M4**<br><i>flt3</i><sup>+</sup> (dendritic cells)<br><i>MPEG1</i><sup>+</sup> cell-containing<br>APC-containing")+
  annotate(geom="segment",x=-1.5,xend=-1.5,y=8,yend=6.5)+
  annotate(geom="richtext",x = -2, y = 8, label = "**M5**<br><i>cfd</i><sup>+</sup>; <i>mrc1a</i><sup>+</sup> cell-containing<br>Mitochondrial activity<br>Non-antigen-presenting")

F5A<-as.ggplot(grid.arrange(Mumap,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))

################################################################################
# PLOT OF APC MARKERS (formerly supp. figure, but works well with macrophages) #
################################################################################

nlevels(as.factor(subset(fshmk_small,primary_ct=="APCs")$gene_id)) # 5 APC markers present in the data
# Add APC marker module score
immune.combined.sct<-AddModuleScore(object = immune.combined.sct,
                                                 features = list(subset(fshmk_small,primary_ct=="APCs")$gene_id),
                                                 ctrl = 10, name = 'APCs')
# make initial violin plot --> we will grab the data from the object
Aplot<-VlnPlot(immune.combined.sct,features="APCs1")

# rename cluster column in clabs for merging purposes
colnames(clabs)[1]<-"ident"
# There are some small inconsistencies in cluster number formatting between the cluster labelling metadata and the Seurat object because of the subclustering
# Not necessary if the code prior to the above UMAP plot has been run
# clabs$ident<-as.character(clabs$ident)
# clabs$ident<-ifelse(clabs$ident=="8","8.0",clabs$ident)
# clabs$ident<-ifelse(clabs$ident=="10","10.0",clabs$ident)
# clabs$ident<-ifelse(clabs$ident=="19","19.0",clabs$ident)
clabs$ident<-as.factor(clabs$ident)

# Merge cluster label metadata with the violin plot data
Aplotdat<-merge(Aplot$data,clabs,by="ident")
# match up level orders so that colours can be matched up correctly
level_order<-levels(Aplotdat$ident)
clabs <- clabs[match(level_order, clabs$ident), ]

# Generate plot...
Aplotdat$Group = factor(Aplotdat$Group, levels=c("Neutrophils","Monocytes / Macrophages","T-cells","B-cells","Other","Unclassified"))
Aplot<-ggplot(Aplotdat,aes(x=Shorthand,y=APCs1,fill=ident))+
  theme_classic()+
  geom_point(position=position_jitter(0.25),size=0.1,alpha=0.2)+
  geom_boxplot(colour="black",outlier.shape=NA,notch=TRUE)+
  geom_hline(yintercept=0,lty="dashed")+NoLegend()+
  labs(x="",y="Module expression score",subtitle="Antigen-presenting cell markers (5):\nCD40, CD74 (x2), MHCII beta-like, MPEG1")+
  theme(plot.subtitle = element_text(hjust = 0.5,size=12))+
  scale_fill_manual(values=clabs$plot_colour)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~Group,scales="free_x",space="free")
# This will be Fig. 5E:
F5E<-as.ggplot(grid.arrange(Aplot,top=textGrob("E", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

lay <- rbind(c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(5,5,5,5,5),
             c(5,5,5,5,5),
             c(5,5,5,5,5))

grobz <- lapply(list(F5A,F5B,F5C,F5D,F5E), ggplotGrob)

pdf("Figure5_UPDATED_FULL_191124.pdf", width = 20, height = 20)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()

# DONE
