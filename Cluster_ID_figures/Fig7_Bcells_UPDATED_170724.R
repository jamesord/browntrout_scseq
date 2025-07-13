# FIGURE 7 B-CELLS UPDATED 14/11/24
# Updated 14/11/24: The code was changed to allow NA2 markers among unique Bcell markers (i.e. DO NOT exclude NA2 markers when filtering for uniqueness)
# However, this did not change the number of unique Bcell markers, so there is no need to regenerate the figures
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

#####################################################
# PRIOR MARKERS THAT ALSO APPEAR AS CLUSTER MARKERS #
#####################################################

subset(clustmarkers,gene_id%in%subset(fshmk_small,primary_ct=="B-cells")$gene_id)%>%distinct(gene_id,name_or_symbol,description)
# 9 including both cd79s

##################
# 7B) UPSET PLOT #
##################

Bcells<-subset(clustmarkers,Shorthand%like%"B"&Shorthand!="RBC")
# subset for those NOT markers of other clusters (EVEN THOSE <60%! THIS IS VERY STRINGENT!)
# also NOTE that we exclude 23 from the exclusion (if that makes sense) because it is an ambiguous cluster with both Bcell and Neutrophil markers
Bcells<-subset(Bcells,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Bcells$cluster&cluster!=23)$gene_id) 
# some rows are duplicated as genes can be markers of lineage or, e.g. antigen presenting cells:
Bcells<-merge(distinct(Bcells[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Bcells[!duplicated(Bcells$gene_id),])
# 113 genes are markers uniquely of the B-cells

# aggregated table
Bcellsagg<-aggregate(Shorthand~gene_id, Bcells, FUN=toString)
Bcellsagg<-merge(Bcellsagg,Bcells[c(1,8,10,9,14)],by="gene_id")%>%distinct()

# UPSET PLOT CODE:
Bcells_upset<-as.data.frame(Bcells[c(13,1)] %>% mutate(value = 1) %>%
                              pivot_wider(
                                names_from = Shorthand,
                                values_from = value, values_fill = 0
                              ))
clusters<-levels(as.factor(Bcells$Shorthand))
Bupset<-as.ggplot(upset(Bcells_upset, clusters))

F7B<-as.ggplot(grid.arrange(Bupset,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))

# do same but for all markers
Bcells2<-subset(clustmarkers_all,Shorthand%like%"B"&Shorthand!="RBC")
Tcell2s<-subset(Bcells2,gene_id%not_in%subset(clustmarkers_all,cluster%not_in%Bcells2$cluster)$gene_id) 
Bcells2<-merge(distinct(Bcells2[-c(10)]),(aggregate(primary_ct~gene_id, fshmk_small, FUN=toString)),by="gene_id",all.x=T) # add aggregated primary_ct column
nrow(Bcells2[!duplicated(Bcells2$gene_id),])
Bcells2agg<-aggregate(Shorthand~gene_id, Bcells2, FUN=toString)
Bcells2agg<-merge(Bcells2agg,Bcells2[c(1,8,10,9,14)],by="gene_id")%>%distinct()
colnames(Bcells2agg)[2]<-"Shorthand_nonHA"
# the following makes an aggregated table of high abundance markers, but lists also clusters for which those genes are markers regardless of abundance
Bcellsagg<-merge(Bcellsagg,Bcells2agg[c(1,2)],by="gene_id",all.x=T,all.y=F)

# See the table Bcellsagg to see which markers are pan-B and which mark specific clusters

colnames(Bcellsagg)[c(3,7)]<-c("clusters_HA","clusters_all")
write.table(Bcellsagg,file="Bcell_markers.tsv",sep="\t",row.names = F,col.names = T,quote=F)

##################
# 7C) PAN-B-CELL #
##################

# No genes are high abundance markers of all B-cell clusters but 9 genes are high abundance markers of at least 6/7

# of interest:
# BCR heavy IgM ENSSTUG00000010446: B6, B4, B5, B7, B1, B2
# BCR light IgL kappa 1 IGKC1 ENSSTUG00000027283: B4, B1, B2, B6, B7, B5
# CD79A ENSSTUG00000031393: B4, B1, B3, B5, B2, B6 (note that CD79B is HA only in five...lacking from B2)
# mef2cb ENSSTUG00000006634: B1, B2, B3, B4, B5, B6 (required for B-cell prolif..see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518613/)

# NOTE: Immunoglobulins to be plotted separately later
# FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000010446",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(axis.line.x = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="unnamed IGHM-like*",subtitle="ENSSTUG00000010446 | B1,B2,B4-7")
# reduced IGHM in B1 and B2, absence in B3
# FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000027283",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
#   theme(plot.title = element_text(hjust = 0))+
#   theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
#   theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="unnamed IGKC1-like*",subtitle="ENSSTUG00000027283 | B1,B2,B4-7")
# reduced IGKC1 in B1 and B2, absence in B3

# CD79 A and B
p1<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000031393",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cd79a*",subtitle="ENSSTUG00000031393 | B1-6")
p2<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000027061",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cd79b*",subtitle="ENSSTUG00000027061 | B1,B3-6")
# NOTE the disparity between cd79a and b
# CD79a seems to be a more robust marker!
# Note reduced expression of CD79a in B7 and of CD79 B in 1,2 and 7
# reduced CD79 expression is a marker of highly differentiated B-cells (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563169/)

# hipk2: KINASE (ENSSTUG00000017457)
# scattered expression throughout B1-6
# involved in cell proliferation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6496359/
# HIPK1 is needed for the function of certain B-cell types in mice: https://pubmed.ncbi.nlm.nih.gov/22545114/
p3<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000017457",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="hipk2",subtitle="ENSSTUG00000017457 | B3-6\nB1,B2 (<60%)")
# erythrocyte protein 4.1-like ENSSTUG00000001843: B1, B2, B3, B4, B5, B6 (involved in B-cell differentiation...see: https://onlinelibrary.wiley.com/doi/full/10.1111/imm.13250)
# suppressor of WNT signalling: https://www.sciencedirect.com/science/article/abs/pii/S0304383520306224
p4<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000001843",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="EPB41-like",subtitle="ENSSTUG00000001843 | B1-6")

# mef2cb
p5<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000006634",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="mef2cb",subtitle="ENSSTUG00000006634 | B1-6")
# NOTE however low expression in B3
# regulator of B-cell homeostasis: https://pubmed.ncbi.nlm.nih.gov/19211936/

# CD37: highly abundant on mature B-cells
p6<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000006678",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD37-like*",subtitle="ENSSTUG00000006634 | B1,B2,B4-6")

# PLOT
F7C<-as.ggplot(grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

# Mention but don't plot:

# the others do not seem to have annotation. So I will blast them...
# ENSSTUG00000000897: hits to Ig-lambda and Ig-kappa chains, but per. ident not so high (60-76%)
# ENSSTUG00000042457: mostly unknown proteins but a couple of hits to Ig-kappa light chain
# ENSSTUG00000042537: mostly unknown proteins but a couple of hits to Ig-kappa
# ENSSTUG00000044159: good hit to Ig-kappa light chain

################################
# 7D: CLUSTER-SPECIFIC MARKERS #
################################

# Plot genes representative of certain expression patterns
# Tricky to choose only nine
# B1 essentially exclusive: zgc:194275 (aka plaat1l)
# Shared expression in B1 and B2: arhgap24
# B2 essentially exclusive: EPB41-like
# B3 (pre-B-cells?): rag1
# B4: carhsp1
# B5: cxcr5
# B6: nfic, pkdr2
# B7 (highly differentiated): DKK4-like

###############
# B1-specific #
###############

# Quite B1-specific...but mostly expressed in upper part
# also scseq Bcell marker in zebrafish: https://elifesciences.org/reviewed-preprints/92424v1
p7<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000017243",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="zgc:194275",subtitle="ENSSTUG00000017243 | B1\nB2 (<60%)")
# aka plaat1:"It is an enzyme that resides in a class of enzymes called phospholipase that hydrolyze phospholipids into fatty acids" (Wiki)
# So this is maybe an enzyme with some digestive function

# See also: unknown protein ENSSTUG00000007017
# top blast hit to: adhesive plaque matrix protein-like [Pelobates fuscus], 47% query cover

######################################
# mirrored expression pattern B1/B2: #
######################################

# B1 with some in B2: arhgap24
p8<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000010099",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="arhgap24",subtitle="ENSSTUG00000010099 | B1\nB2 (<60%)")
# enriched in memory and naive B-cells: https://www.proteinatlas.org/ENSG00000138639-ARHGAP24/immune+cell
# Inhibits cell proliferation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7197742/

# See also:
# B1 and some in B2: MAST4: microtubule-associated serine/threonine protein kinase ENSSTUG00000014200

# B1 and B2: calcium/calmodulin-dependent protein kinase type 1D-like ENSSTUG00000017654
# overexpressed in B-cell lymphoma
# regulation of neutrophil activation / chemotaxis (NCBI)
# may be involved in regulation of granulocyte function (Wiki)

# B2 and some in B1: atp2a3 ENSSTUG00000017361
# "catalyzes the hydrolysis of ATP coupled with the translocation of calcium" (source?)

# B2 and some in B1/6: ctnna1 ENSSTUG00000037275
# Recognised role in cell adhesion; DNA damage response (WNT signalling): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5450192/

###############
# B2-specific #
###############

# B2: band 4.1-like protein 1 (aka EPB41) --> quite specific to B2 but some in B1
p9<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000043834",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="EPB41-like",subtitle="ENSSTUG00000043834 | B2\nB1,B6 (<60%)")
# suppressor of WNT signalling: https://www.sciencedirect.com/science/article/abs/pii/S0304383520306224
# involved in cell migration / proliferation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7085237/
# NOTE the other EPB41-like that is a supposed pan-B-cell marker

# See also:
# mikado gene mikado.19G1247 with top blast hit to NHSL2: B2, B6
# involved in cell migration: https://www.nature.com/articles/s41467-023-39276-w
# "Cell migration is a critical process for animal cells, especially in the embryo, but also in the adult.
# For example, immune cells constantly patrol the organism to fight infections."

###############
# B3-specific #
###############

# B3: RAG1 --> VERY IMPORTANT GENE!!!
# expressed in developing B-cells https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2676217/
p10<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000005374",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="rag1",subtitle="ENSSTUG00000005374 | B3\nB5,B6 (<60%)")
# Note that RAG2 was not detected amongst cluster markers

# See also:
# B3: akap12: KINASE-anchoring protein ENSSTUG00000018953
# Role in B-cell development https://ashpublications.org/blood/article/120/21/855/87490/RNAi-Screening-Identifies-A-Novel-Role-for-A

# B3: itga5 (aka CD49e) ENSSTUG00000030861: upregulated in pre-b-cells: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5066240/
# however also involved with b-cell activation: https://pubmed.ncbi.nlm.nih.gov/36016494/
# marker of MARGINATING B-cells: https://pubmed.ncbi.nlm.nih.gov/34313733/

# B3: tek (aka CD202B) ENSSTUG00000017124: tyrosine KINASE receptor, in mammals expressed in epithelial cells
# in humans also dendritic cell maturation marker?
# Role of tyrosine kinases in pre-b-cell development: https://www.sciencedirect.com/science/article/abs/pii/B9780128002674000043

# B3: stk35l ENSSTUG00000031221: KINASE involved in cell cycle and cell migration
# B3: unnamed, dedicator of cytokinesis protein 1-like (DOCK1) ENSSTUG00000040345
# B3/4: mikado.22G195 unnamed Ig-like (IGLL1-like)

# B3 and B5: histone H1.0-B-like: suggests cell proliferation?
# B3: ENSSTUG00000040345 dedicator of cytokinesis protein 1-like

# Note also absence of CD37

###############
# B4-specific #
###############

# B4: carhsp1: calcium-regulated heat stable protein 1-like
# not much else that is very specific to B4 (i.e. not also expressed in other clusters)
# we know that this is a highly proliferative cluster
p11<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000026177",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="carhsp1",subtitle="ENSSTUG00000026177 | B4\nB5 (<60%)")

###############
# B5-specific #
###############

# B5: cxcr5...also a lot in B1 (but <60%)
# essential role in B-cell migration (Wiki)
# B-cell activation: https://pubmed.ncbi.nlm.nih.gov/21659539/
# B-helper cells https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2193094/
# rapidly activated in T-cells upon activation: https://onlinelibrary.wiley.com/doi/pdf/10.1002/eji.200425478
p12<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000007343",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="cxcr5",subtitle="ENSSTUG00000007343 | B5\nB1 (<60%)")
# If cxcr5 is a marker of activation, expression patterns suggests B1 and B2 contain a mixture of activated and quiescent cells

# See also:
# B5: unnamed ASGR1-like ENSSTUG00000020630; could not find anything related to B-cells
# B5 and B1: IGLC1-like ENSSTUG00000015148: note pattern in B1 and B2 as well

# dapp1 ENSSTUG00000034936: expression shared between B1,4, and 5 (with some in 2)
# role in B-cell proliferation: https://pubmed.ncbi.nlm.nih.gov/14588241/

######
# B6 #
######

# B3,B5,B6: xpa (DNA repair)
p13<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000026759",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="xpa",subtitle="ENSSTUG00000026759 | B3,B5,B6\nB4 (<60%)")

# B6: pkdr2 is the only one which is high abundance exclusively in B6, but it is also expressed in 1,2,and 5-
# important for regulatory B-cell function: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8761795/
p14<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000001609",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="pkdr2",subtitle="ENSSTUG00000001609 | B6\nB1,B2,B5 (<60%)")

# See also:
# B4 & 6: nfic ENSSTUG00000044896

###############
# B7-specific #
###############

# B7: dickkopf-related protein 4
# wnt signalling
p15<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000032962",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="DKK4-like",subtitle="ENSSTUG00000032962 | B7")
# see also DKK3: https://journals.aai.org/jimmunol/article/194/6/2624/104857/Dickkopf-3-Acts-as-a-Modulator-of-B-Cell-Fate-and

# See also:
# B7, some in B1/2: BCR light IgL lambda IGLC1 ENSSTUG00000001951
# B7, some in B1: uncharacterised gene, lots of top BLAST hits to immunoglobulin light chain (IGLC) ENSSTUG00000027299
# Note also reduced expression of CD79a&b in B7
# reduced CD79 expression is a marker of highly differentiated B-cells (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563169/)
# NOTE however also absence of CD37

F7D<-as.ggplot(grid.arrange(p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=3,top=textGrob("D", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

######################
# 7E IMMUNOGLOBULINS #
######################

# Showing representative expression patterns of Ig-like genes
# PATTERN1: ENSSTUG00000010446 BCR heavy IgM
# PATTERN2: ENSSTUG00000027283 BCR light IgL kappa 1 IGKC1
# PATTERN3: ENSSTUG00000040813
# PATTERN4: ENSSTUG00000040711
# PATTERN5: mikado.22G195

# PATTERN1: Expression throughout B-cells except for B3
# ENSSTUG00000010446 BCR heavy IgM: B6, B4, B5, B7, B1, B2 (weak expression THROUGHOUT B1 and 2)
p16<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000010446",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="IGHM-like*",subtitle="ENSSTUG00000010446 | B1-2,B4-7")
# PATTERN2: Expression tapering in B1 and 2 in parallel fashion
# ENSSTUG00000027283 BCR light IgL kappa 1 IGKC1: B4, B1, B2, B6, B7, B5 (TAPERING expression B1 and 2)
p17<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000027283",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="IGKC1-like*",subtitle="ENSSTUG00000027283 | B1-2,B4-7")
# PATTERN3: Expression THROUGHOUT B1 but tapering or absent in B2
# ENSSTUG00000040813: expressed THROUGHOUT B1 (but also 4,5,7); ENSEMBL lists immunoglobulin-like domains (specifically V-set) in detected protein domains
p18<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000040813",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="V-set domain-containing",subtitle="ENSSTUG00000040813 | B1,B4-5,B7")
# PATTERN4:
# ENSSTUG00000040711 concentrated in B7, tapering / bimodal expression in B1, some also in B4/5: ENSEMBL lists immunoglobulin-like domains (specifically C1-set) in detected protein domains; top blast hits to immunoglobulin lambda-like; top blast hit to Ig light chain (O.mykiss, 92.34% ident)
p19<-FeaturePlot(immune.combined.sct.downsampled,features="ENSSTUG00000040711",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="C-set domain-containing",subtitle="ENSSTUG00000040711 | B1,B4-5,B7")
# PATTERN5: B3/4
# mikado.22G195 Ig-like domain-containing protein (IGLL1-like): B3/4
p20<-FeaturePlot(immune.combined.sct.downsampled,features="mikado.22G195",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="IGLL1-like",subtitle="mikado.22G195 | B3,B4") 

F7E<-as.ggplot(grid.arrange(p16,p17,p18,p19,p20,nrow=1,top=textGrob("E", x = 0, hjust = 0,gp=gpar(cex=2))))+
  theme(panel.border=element_rect(fill=NA))

#####################
# 7A SUMMARY FIGURE #
#####################

## SUMMARY UMAP FIGURE
uplot<-UMAPPlot(immune.combined.sct.downsampled)
uplotdat<-uplot$data
colnames(uplotdat)[3]<-"cluster"
# need to manually alter some numbers for compatibility
clabs$cluster<-ifelse(clabs$cluster=="8","8.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="10","10.0",clabs$cluster)
clabs$cluster<-ifelse(clabs$cluster=="19","19.0",clabs$cluster)
uplotdat<-merge(uplotdat,clabs,by="cluster")
uplotdat<-subset(uplotdat,UMAP_1> 0&UMAP_2<10)
uplotdat$Shorthand<-ifelse(uplotdat$Name%like%"B-",uplotdat$Shorthand,"Other")
levels(as.factor(uplotdat$Shorthand))

# had to set these manually :/
Bcolors<-c("darkred","violetred2","mediumvioletred","palevioletred1","orangered3","palevioletred4","indianred","grey")

Bumap<-ggplot(uplotdat,aes(x=UMAP_1,y=UMAP_2))+
  theme_classic()+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  geom_point(aes(color=Shorthand),size=0.1)+
  scale_color_manual(values=Bcolors)+
  NoLegend()+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="")+
  # ANNOTATIONS
  annotate(geom="segment",x=13,xend=10,y=-5,yend=-5)+
  annotate(geom="segment",x=13,xend=13,y=-5,yend=0)+
  annotate(geom="richtext",x = 13, y = -5, label = "**< B1 | B2 ^**<br>Mature B-cells<br>CD37-like<sup>+</sup><br>reduced <i>cd79b</i><br><i>cxcr5</i><sup>+</sup> cell-containing<br>APC-containing<br>Actin binding")+
  annotate(geom="segment",x=6,xend=9.25,y=8,yend=7)+
  annotate(geom="richtext",x = 6, y = 8, label = "**B3**<br>Pro/pre-B-cells<br>CD37-like<sup>-</sup><br><i>rag1</i><sup>+</sup>")+
  annotate(geom="segment",x=4,xend=6.5,y=3,yend=2.5)+
  annotate(geom="richtext",x = 4, y = 3, label = "**B4**<br>Proliferating<br>(<i>mki67</i><sup>+</sup>, <i>pcna</i><sup>+</sup>)<br>CD37-like<sup>+</sup>")+
  annotate(geom="segment",x=4,xend=8.5,y=0,yend=0)+
  annotate(geom="richtext",x = 4, y = 0, label = "**B5**<br>CD37-like<sup>+</sup><br><i>cxcr5</i><sup>+</sup><br>Actin binding")+
  annotate(geom="segment",x=13,xend=10.25,y=4,yend=3)+
  annotate(geom="richtext",x = 13, y = 4, label = "**B6**<br>Possible regulatory<br><i>pkdr2</i><sup>+</sup>, CD37-like<sup>+</sup><br>Actin binding")+
  annotate(geom="segment",x=3,xend=6,y=-3,yend=-2)+
  annotate(geom="richtext",x = 3, y = -3, label = "**B7**<br>CD37-like<sup>-</sup><br>reduced <i>cd79</i>")

F7A<-as.ggplot(grid.arrange(Bumap,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))

lay <- rbind(c(1,1,2,2,2),
             c(1,1,2,2,2),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(3,3,4,4,4),
             c(5,5,5,5,5))

grobz <- lapply(list(F7A,F7B,F7C,F7D,F7E), ggplotGrob)

pdf("Figure7_UPDATED_170724.pdf", width = 20, height = 25)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()