# plotting per-individual cell counts for each cluster
# for supplementary figure
# UPDATED 29/12/24

library(ggplot2)

clustcounts<-read.table("Seurat_cluster_ID/current/cluster_individual_counts_subclusters_040224.txt",sep="\t",header=T)[-c(5)]
head(clustcounts)
clabs<-read.table("Seurat_cluster_ID/current/cluster_classifications_120724.txt",sep="\t",header=T)
head(clabs)

table(clustcounts$individual)
# A6, B6 and C8 are all missing from one cluster
levels(as.factor(clustcounts$cluster))
levels(as.factor(subset(clustcounts,individual=="A6")$cluster))
levels(as.factor(subset(clustcounts,individual=="B6")$cluster))
levels(as.factor(subset(clustcounts,individual=="C8")$cluster))
# in each case it is cluster 28 from which they are missing
# add the missing '0' counts to a new dataframe and merge that with the main one
clustcounts<-merge(rbind(clustcounts,
                   data.frame(
                     cluster=rep(28,3),
                     individual=c("A6","B6","C8"),
                     group=c("mix","mix","frm"),
                     batch=c("Batch_A","Batch_B","Batch_C"),
                     cluster_individual_count=rep(0,3),
                     individual_count=c(8822,9633,12421)
                   )),clabs,by="cluster")
clustcounts$cluster_individual_frac<-clustcounts$cluster_individual_count/clustcounts$individual_count

# reorder factor levels
clustcounts$group<-factor(clustcounts$group,levels=c("wld","mix","frm"))

head(clustcounts)

# code for two  sets of plots: one for fraction of cells in cluster
pdf("SupFigure_cell_individual_fractions_291224.pdf", width = 12, height = 10)
ggplot(clustcounts,aes(x=group,y=cluster_individual_frac,shape=batch,fill=group))+
  theme_classic()+
  #geom_line(aes(x = group, group=batch), position = "identity",colour="black")+
  geom_point(size=3,position=position_jitter(0.25)) +
  scale_fill_manual(values=c("darkblue","deeppink","yellow"),guide = "none")+
  scale_shape_manual("Batch:",values=c(22,21,24),labels=c("A","B","C"))+
  labs(x="",y="Frac. cells in cluster")+
  scale_x_discrete(labels=c("Wild","Mix","Farm"))+
  facet_wrap(.~Shorthand,scales="free")
dev.off()

# second plot for total cell counts
pdf("SupFigure_cell_individual_counts_291224.pdf", width = 12, height = 10)
ggplot(clustcounts,aes(x=group,y=cluster_individual_count,shape=batch,fill=group))+
  theme_classic()+
  #geom_line(aes(x = group, group=batch), position = "identity",colour="black")+
  geom_point(size=3,position=position_jitter(0.25)) +
  scale_fill_manual(values=c("darkblue","deeppink","yellow"),guide = "none")+
  scale_shape_manual("Batch:",values=c(22,21,24),labels=c("A","B","C"))+
  labs(x="",y="N cells in cluster")+
  scale_x_discrete(labels=c("Wild","Mix","Farm"))+
  facet_wrap(.~Shorthand,scales="free")
dev.off()

# Looking at cell counts in candidate clusters for pseudobulk analyses
subset(clustcounts,Shorthand=="B3")
subset(clustcounts,Shorthand=="T5") # select T5 for T-cells because numbers are relatively consistent between individuals
subset(clustcounts,Shorthand=="M1")
subset(clustcounts,Shorthand=="N1")