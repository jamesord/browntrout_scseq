# MARKER GENES FUNCTIONING IN MULTIPLE CELL TYPES

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

# lyz: N1 (neutrophil) plus two macrophage clusters (low abundance)
p1<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000025205",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="lyz*",subtitle="N1\nM3 (<60%)")

# CD74: a marker of monocytes / APCs that is also expressed in B-cells
p2<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000006024",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "plain"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="CD74-like*",subtitle="ENSSTUG00000006024 | multiple")

# two genes that show a similar pattern and expressed in multiple specific clusters
# (central cluster, T-cells, some B-cells, peripheral neutrophil)

# ribosomal protein L12
p3<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000013074",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="rpl12",subtitle="ENSSTUG00000013074 | multiple")
# see role of ribosomal diversity in differentiation:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7265316/
# see also rpl12 has been found to be enriched in human embryonic stem cells:
# https://www.nature.com/articles/s41419-022-04635-w
# also see lots of activity in T-cells and some B-cells
# protein synthesis in T-cells: https://www.nature.com/articles/s41423-021-00792-8

# upf3: nonsense mediated decay
p4<-FeaturePlot(immune.combined.sct,features="ENSSTUG00000019183",min.cutoff = "q10", max.cutoff = "q90",raster=T)+
  theme(plot.title = element_text(hjust = 0))+theme(plot.title = element_text(face = "italic"))+
  theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+labs(x="",y="",title="upf3b",subtitle="ENSSTUG00000019183 | multiple")
# see here, role in T-cells: https://www.embopress.org/doi/full/10.1038/sj.emboj.7601628
# see also role of upf/NMD in hematopoiesis: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2377192/
# note that the above two appear to be inversely correlated with proliferation markers

FSA<-as.ggplot(grid.arrange(p1,nrow=1,top=textGrob("A", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
FSB<-as.ggplot(grid.arrange(p2,nrow=1,top=textGrob("B", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
FSC<-as.ggplot(grid.arrange(p3,nrow=1,top=textGrob("C", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))
FSD<-as.ggplot(grid.arrange(p4,nrow=1,top=textGrob("D", x = 0, hjust = 0,gp=gpar(cex=2))))+theme(panel.border=element_rect(fill=NA))

pdf("FigureS_multicluster_markers_FULL_170724.pdf", width = 15, height = 4)
grid.arrange(FSA,FSB,FSC,FSD,nrow=1)
dev.off()