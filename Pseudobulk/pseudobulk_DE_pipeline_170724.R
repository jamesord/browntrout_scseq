# DIFFERENTIAL EXPRESSION TESTING WITH PSEUDOBULK COUNTS 17/07/24
# EDITED 17/11/24 to fix labeling issue in heatmap

rm(list = ls()) # rm is the command to remove the list

library(edgeR)
library("data.table")
library("ComplexHeatmap")
library("colorRamps")
library(ggplot2)
library(gridExtra)
library(stringr)
library(RColorBrewer)
library("ggthemes")
library(scales)
library("pals")
library("VennDiagram")
library(ggplotify)
library(biomaRt)
library(data.table)
library(tidyr)
library(clusterProfiler)
library(Seurat)
library("enrichplot")
library(ggh4x)
# z-score function
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
getwd()
# As we want to run the following for pseudobulk data from FOUR clusters, the below code comprises one very long function
# Within this function:
# - Process count data from input cluster
# - Generate PCA plot
# - Run differential expression testing and GO term enrichment
#   -> Generate output files containing the results of the above
# - Generate heatmap
# - Generate a mini UMAP plot highlighting the cluster of interest

# function has three input parameters: cluster number (integer), colour to use on UMAP (string) and UMAP title / cluster name (string)
do_pseudobulk_DE_stuff<-function(cluster_num,cluster_colour,UMAP_title)
{
  ############################
  # PROCESS COUNT (META)DATA #
  ############################
  
  # read in pseudobulk counts as DGElist object and add group information
  group <- factor(rep(c("Wild","Mix","Farm"),3))
  RG <- DGEList(counts=read.delim(paste("Pseudobulk_DE/pseudobulk_counts_cluster",cluster_num,".txt",sep="")),group=group)
  
  # filter to retain only genes with at least 10 reads in some samples
  keep.exprs <- filterByExpr(RG, group=group, min.count = 10)
  dgList <- RG[keep.exprs,, keep.lib.sizes=FALSE]
  
  # normalization
  dgList <- calcNormFactors(dgList, method="TMM")
  
  ############
  # PCA plot #
  ############
  
  pca<-plotMDS(dgList,gene.selection="common"); # base plot
  
  # then make it nicer, first set up dataframe: lots of information to fit on the plot...
  pca <- data.frame(x=pca$x,y=pca$y,
                    group=rep(c("Wild","Mix","Farm"),3),
                    batch=c(rep("A",3),rep("B",3),rep("C",3)))
  pca$group<-factor(pca$group,levels=c("Wild","Mix","Farm"))
  ploT_pca<-ggplot(pca, aes(x = x, y=y,fill=group,shape=batch)) +
    geom_point(size=4,color="black") + 
    scale_fill_manual("Origin:",values=c("darkblue","deeppink","yellow"))+
    scale_shape_manual("Batch:",values=c(22,21,24))+
    labs(x = "PC 1", y = "PC 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
    theme(legend.position="right")+
    scale_x_continuous(labels = label_number(accuracy = 0.1))+
    scale_y_continuous(labels = label_number(accuracy = 0.1))+
    guides(fill=guide_legend(override.aes=list(shape=21))) # necessary to overcome a bug preventing colours displaying on the legend
  
  ###################################
  # DIFFERENTIAL EXPRESSION TESTING #
  ###################################
  
  # set the design matrix
  batch<-factor(c(rep("A",3),rep("B",3),rep("C",3)))
  gtype<-factor(rep(c("W","M","F"),3),levels=c("W","M","F"))
  design <- model.matrix(~batch+gtype);rownames(design) <- colnames(dgList)
  
  # estimate dispersion
  y <- estimateDisp(dgList, design, verbose=TRUE)
  fit <- glmFit(y, design)
  fit$coefficients
  
  # loop for two DE contrasts (wild v mix and wild v farm)
  cpm_log_all<-NULL
  toptags_all<-NULL
  for (i in c(4,5)){
    # perform likelihood ratio test on the selected coefficient (4 or 5)
    lrt <- glmLRT(fit,coef=i)
    # get the numbers of DEGs
    desum<-summary(decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05))
    # put DEGs into a dataframe
    toptags<-as.data.frame(topTags(lrt, n=desum[1]+desum[3], adjust.method="fdr", sort.by="PValue"))
    # assign either farm or mix, depending on the coefficient
    toptags$contrast<-ifelse(i==4,"mix","farm")
    # add gene id column
    toptags$gene_id<-rownames(toptags);rownames(toptags)<-NULL
    # append results from two coefficients
    toptags_all<-rbind(toptags_all,toptags)
    # get counts per million of DEGs
    top <- rownames(topTags(lrt,desum[1]+desum[3]))
    cpm_log <- cpm(y, log = TRUE)
    cpm_log_all<-rbind(cpm_log_all,cpm_log)
  }
  # deduplicate CPM
  cpm_log_all<-cpm_log_all[!duplicated(cpm_log_all),]
  
  # get gene names and descriptions
  
  trutta_names<-read.delim("annotation_stuff/all_gnames_RNAassay.txt",sep="\t",header=T)
  
  # prepare the table of differentially expressed genes
  toptags_all<-merge(toptags_all,trutta_names,by="gene_id",all.x=TRUE)
  # if still no gene name for ENSEMBL () or mikado gene (NA), use gene id instead
  toptags_all$name_or_symbol<-ifelse(toptags_all$name_or_symbol=="",toptags_all$gene_id,toptags_all$name_or_symbol)
  toptags_all$name_or_symbol<-ifelse(is.na(toptags_all$name_or_symbol),toptags_all$gene_id,toptags_all$name_or_symbol)
  
  # sort by group and FDR
  toptags_all<-toptags_all[with(toptags_all, order(contrast, FDR)),]
  
  # To this list of DEGs we want to add some extra overlaps:
  # --> are any of the DEGs prior immune cell markers?
  # --> do any of the DEGs have ohnologs that are also differentially expressed?
  
  fshmk<-read.delim("cell_marker_compilation/fish_compilation/all_potential_markers_long_290124.txt",sep="\t",header=T)
  # we will ONLY want to highlight IMMUNE CELL markers, which means removing some of these...
  fshmk<-subset(fshmk,primary_ct!="Thrombocytes")
  fshmk<-subset(fshmk,primary_ct!="Erythroid_cells")
  fshmk<-subset(fshmk,primary_ct!="Endothelial_cells")
  fshmk<-subset(fshmk,primary_ct!="HSCs")
  fshmk<-subset(fshmk,primary_ct!="Proliferating_cells")
  fshmk<-subset(fshmk,primary_ct!="All_cells")
  # add information (simple YES or NO if prior immune marker) to dataframe
  toptags_all$prior_immune_marker<-ifelse(toptags_all$gene_id %in% fshmk$gene_id,"YES","NO")
  
  # get ohnolog pair information
  ohno<-read.delim("Neofunctionalisation/RNAassay_genes_ohnolog_pairs.txt",sep="\t",header=T)
  colnames(ohno)[1]<-"gene_id"
  # add ohnolog pair ID to dataframe
  toptags_all<-merge(toptags_all,ohno,by="gene_id",all.x=T,all.y=F)
  # count the occurrences of a pair member within a given experimental contrast
  toptags_all<-toptags_all%>%add_count(contrast,elucius_homolog_ensembl_gene,name="n_ohnologs_DE_contrast")
  # not pretty but does the job...
  toptags_all$n_ohnologs_DE_contrast<-ifelse(is.na(toptags_all$elucius_homolog_ensembl_gene),NA,toptags_all$n_ohnologs_DE_contrast)
  toptags_all$both_ohnologs_DE_contrast<-ifelse(toptags_all$n_ohnologs_DE_contrast==2,"YES","NO")
  toptags_all$n_ohnologs_DE_contrast<-NULL
  # it seems to be quite rate that both ohnologs in a pair are DE in the same contrast, so I don't think we need to come up with some elaborate labelling for that
  # we can point it out in the results, however.
  
  # write out
  #write.table(toptags_all,file=paste("Pseudobulk_DE/cluster",cluster_num,"_DEGs_170723.txt",sep=""),sep="\t",quote=F,row.names = F)
  
  ###############################
  # GO term enrichment analysis #
  ###############################
  
  # get GO terms from file, term to gene and term to name mappings
  T2G2N<-subset(read.delim("annotation_stuff/all_GOterms_RNAassay.txt",sep="\t",header=T),ontology=="MF")[c(2,1,4)]
  T2G<-T2G2N[c(1,2)]
  T2N<-T2G2N[c(1,3)]
  
  # define universe
  universe<-row.names(dgList$counts)
  
  # run enricher using four input gene lists (two, up and downregualted, per contrast)
  farmgenesup_enrch<-(enricher(gene=subset(toptags_all,contrast=="farm"&logFC>0)$gene_id,universe=universe,TERM2GENE=T2G,TERM2NAME=T2N))@result;farmgenesup_enrch$DEG_group<-"farm_up"
  farmgenesdown_enrch<-(enricher(gene=subset(toptags_all,contrast=="farm"&logFC<0)$gene_id,universe=universe,TERM2GENE=T2G,TERM2NAME=T2N))@result;farmgenesdown_enrch$DEG_group<-"farm_down"
  mixgenesup_enrch<-(enricher(gene=subset(toptags_all,contrast=="mix"&logFC>0)$gene_id,universe=universe,TERM2GENE=T2G,TERM2NAME=T2N))@result;mixgenesup_enrch$DEG_group<-"mix_up"
  mixgenesdown_enrch<-(enricher(gene=subset(toptags_all,contrast=="mix"&logFC<0)$gene_id,universe=universe,TERM2GENE=T2G,TERM2NAME=T2N))@result;mixgenesdown_enrch$DEG_group<-"mix_down"
  # concatenate GO term results from each input gene list and filter according to:
  # - FDR < 0.05
  # - At least three genes per enriched GO term
  enrichment_results_all<-subset(rbind(farmgenesup_enrch,farmgenesdown_enrch,mixgenesup_enrch,mixgenesdown_enrch),p.adjust<0.05&Count>2)[c(1:4,6,8:10)]
  # the above all assumes that the 'results' slot of the enricher object is not empty...otherwise there will be an error when trying to retrieve it.
  # as long as there are even non-significant GO terms in there it is not a problem.
  
  # convert results to long format
  enrichment_results_long<-enrichment_results_all %>% separate_rows(geneID,sep = "/")
  
  # write results out to a file
  #write.table(enrichment_results_all,file=paste("Pseudobulk_DE/cluster",cluster_num,"_GOterms_170723.txt",sep=""),sep="\t",quote=F,row.names = F)
  
  ###########
  # HEATMAP #
  ###########
  
  # Heatmap row annotation:
  # - add a star to the gene name if it is an immune marker
  # - embolden and colour text if the gene belongs to an enriched GO term
  
  # Get CPMs of DEGs with cpm < 0.01 (for plotting)
  toptagsdedup<-toptags_all[!duplicated(toptags_all$gene_id),]
  toptags_subset<-subset(toptagsdedup,FDR<0.01)
  #toptags_subset<-toptagsdedup # or just use toptagsdedup if we don't want to filter
  cpm_log_subset<-cpm_log_all[toptags_subset$gene_id,]
  
  data_subset_norm <- as.data.frame(t(apply(cpm_log_subset, 1, cal_z_score)))
  # ADD GENE NAMES TO CPM Z SCORES
  data_subset_norm$gene_id<-rownames(data_subset_norm)
  data_subset_norm<-merge(data_subset_norm,trutta_names[c(1,2)],by="gene_id",all.x=T)
  # if no gene name (empty cell or NA) use id instead
  data_subset_norm$name_or_symbol<-ifelse(data_subset_norm$name_or_symbol=="",data_subset_norm$gene_id,data_subset_norm$name_or_symbol)
  data_subset_norm$name_or_symbol<-ifelse(is.na(data_subset_norm$name_or_symbol),data_subset_norm$gene_id,data_subset_norm$name_or_symbol)
  # in case of duplicate gene names...
  # https://stackoverflow.com/questions/72701949/iterate-through-a-column-for-duplicates-append-a-number-to-duplicate-column-val
  data_subset_norm <- data_subset_norm |> 
    mutate(dupl = ifelse(duplicated(name_or_symbol), 1, 0)) |> 
    group_by(name_or_symbol) |> 
    mutate(dupl = cumsum(dupl),
           name_or_symbol = paste(name_or_symbol, dupl, sep = "_d")) |> 
    select(-dupl)
  data_subset_norm$name_or_symbol <- gsub("_d0","",data_subset_norm$name_or_symbol)
  # convert back to data.frame
  data_subset_norm<-as.data.frame(data_subset_norm)
  # add the star to the name to signify prior immune marker
  data_subset_norm$name_or_symbol<-ifelse(data_subset_norm$gene_id%in%subset(toptags_all,prior_immune_marker=="YES")$gene_id,paste(data_subset_norm$gene_id,"*",sep=""),data_subset_norm$name_or_symbol)
  
  # row names need to be the names of the genes
  row.names(data_subset_norm)<-data_subset_norm$name_or_symbol
  
  # Add enriched GO term annotations to the dataframe - we will highlight genes associated with three GO terms
  # It does not matter if the GO term is not enriched in the cluster being run in the function.
  data_subset_norm$GO_annot<-NA
  data_subset_norm$GO_annot<-ifelse(data_subset_norm$gene_id %in% subset(enrichment_results_long,ID=="4930")$geneID,"Gprotein_coupled_receptor_activity",data_subset_norm$GO_annot)
  data_subset_norm$GO_annot<-ifelse(data_subset_norm$gene_id %in% subset(enrichment_results_long,ID=="20037"|ID=="5506")$geneID,"iron_or_heme_binding",data_subset_norm$GO_annot)
  data_subset_norm$GO_annot<-ifelse(data_subset_norm$gene_id %in% subset(enrichment_results_long,ID=="5509")$geneID,"calcium_ion_binding",data_subset_norm$GO_annot)
  
  # select columns
  data_subset_norm_GOannots<-data_subset_norm[12] # separate dataframe for the annotations
  data_subset_norm<-data_subset_norm[c(2:10)] # main dataframe with only row names and CPM Z-scores for the heatmap
  
  # Establish row annotation metadata for the heatmap
  # https://stackoverflow.com/questions/65966708/complexheatmap-highlight-specific-rows
  # Rows to highlight
  ironrows <- rownames(subset(data_subset_norm_GOannots,GO_annot=="iron_or_heme_binding"))
  calcrows <- rownames(subset(data_subset_norm_GOannots,GO_annot=="calcium_ion_binding"))
  gprorows <- rownames(subset(data_subset_norm_GOannots,GO_annot=="Gprotein_coupled_receptor_activity"))
  # Set stylings for row names and make our selected rows unique
  row_idx1 <- which(rownames(data_subset_norm) %in% ironrows)
  row_idx2 <- which(rownames(data_subset_norm) %in% calcrows)
  row_idx3 <- which(rownames(data_subset_norm) %in% gprorows)
  fontsizes <- rep(10, nrow(data_subset_norm))
  fontcolors <- rep('black', nrow(data_subset_norm))
  fontcolors[row_idx1] <- 'red'
  fontcolors[row_idx2] <- 'purple'
  fontcolors[row_idx3] <- 'darkorange3'
  fontfaces <- rep('plain',nrow(data_subset_norm))
  fontfaces[row_idx1] <- 'bold'
  fontfaces[row_idx2] <- 'bold'
  fontfaces[row_idx3] <- 'bold'
  rowAnno <- rowAnnotation(rows = anno_text(rownames(data_subset_norm), gp = gpar(fontsize = fontsizes,fontface = fontfaces, col = fontcolors)))
  
  # generate heatmap (ComplexHeatmap)
  hm1<-Heatmap(as.matrix(data_subset_norm),col=viridis(n=500),
               rect_gp = gpar(col = "grey", lwd = 1),row_names_gp = gpar(fontsize = 6),
               heatmap_legend_param = list(direction = "horizontal",title="Z-score"),
               right_annotation = rowAnno, show_row_names = FALSE,
               top_annotation = HeatmapAnnotation(Origin = rep(c("Wild","Mix","Farm"),3),gp = gpar(col = "grey"),
                                                  col=list(Origin = c(Wild="darkblue", Mix="deeppink", Farm="yellow")),
                                                  annotation_legend_param = list(Origin = list(nrow = 1,labels=c("F","M","W"),border="grey"))))
  # convert to ggplot object
  heatmap<-as.ggplot(grid.grabExpr(draw(hm1,merge_legend = T, heatmap_legend_side = "bottom",annotation_legend_side = "bottom")))
  
  #############################################
  # UMAP plot showing the cluster in question #
  #############################################
  
  # load Seurat data (subsample is sufficient as it is only a mini subfigure)
  load("C:/Users/James/Documents/R/trout_local/trout_SCTprocessed_PC30k30_var.ft.filt_finalclusters_180224.downsamp.rds")
  DefaultAssay(immune.combined.sct.downsampled)<-"integrated"
  
  # generate base UMAP and extract the coordinate and expression data from the object
  uplot<-UMAPPlot(immune.combined.sct.downsampled)
  uplotdat<-uplot$data
  # replace cluster IDs with binary variable: either a cell is from the cluster of interest, or from some other cluster
  colnames(uplotdat)[3]<-"cluster"
  uplotdat$cluster<-ifelse(uplotdat$cluster==cluster_num,"Cluster","Other")
  # generate UMAP highlighting only the cluster of interest
  plot_UMAP<-ggplot(uplotdat,aes(x=UMAP_1,y=UMAP_2))+
    theme_classic()+
    theme(axis.line.x = element_blank())+theme(axis.line.y = element_blank())+
    geom_point(aes(color=cluster),size=0.1)+
    scale_color_manual(values=c(cluster_colour,"grey"))+
    NoLegend()+
    theme(axis.text = element_blank(), axis.ticks = element_blank())+
    labs(x="",y="",title=UMAP_title)
  
  # save outputs as list
  outlist<-list()
  outlist$ploT_pca<-ploT_pca
  outlist$toptags_all<-toptags_all
  outlist$heatmap<-heatmap
  outlist$enrichment_results_all<-enrichment_results_all
  outlist$enrichment_results_long<-enrichment_results_long
  outlist$plot_UMAP<-plot_UMAP
  outlist$data_subset_norm<-data_subset_norm
  outlist$cpm_log_all<-cpm_log_all
  
  return(outlist)
}

# Run the function on count data from four clusters
pseudo_test_cl1<-do_pseudobulk_DE_stuff(1,"darkgreen","Neutrophils 1")
pseudo_test_cl5<-do_pseudobulk_DE_stuff(5,"darkorange3","Macrophages 1")
pseudo_test_cl12<-do_pseudobulk_DE_stuff(12,"mediumvioletred","B-cells 3")
pseudo_test_cl21<-do_pseudobulk_DE_stuff(21,"royalblue","T-cells 5")

pseudo_test_cl1$enrichment_results_all$cluster<-"N1"
pseudo_test_cl5$enrichment_results_all$cluster<-"M1"
pseudo_test_cl12$enrichment_results_all$cluster<-"B3"

enrichment_results_all<-rbind(pseudo_test_cl1$enrichment_results_all,pseudo_test_cl5$enrichment_results_all,pseudo_test_cl12$enrichment_results_all)
write.table(enrichment_results_all,file="Pseudobulk_DE/pseudobulk_GO_enricher_compiled.txt",sep="\t",quote=F,row.names = F)

nrow(subset(pseudo_test_cl1$toptags_all,contrast=="mix"))
nrow(subset(pseudo_test_cl1$toptags_all,contrast=="farm"))

nrow(subset(pseudo_test_cl5$toptags_all,contrast=="mix"))
nrow(subset(pseudo_test_cl5$toptags_all,contrast=="farm"))

nrow(subset(pseudo_test_cl12$toptags_all,contrast=="mix"))
nrow(subset(pseudo_test_cl12$toptags_all,contrast=="farm"))

nrow(subset(pseudo_test_cl21$toptags_all,contrast=="mix"))
nrow(subset(pseudo_test_cl21$toptags_all,contrast=="farm"))


##############################
# GENERATE COMPOSITE FIGURES #
##############################

# Main figure to include clusters 1,5, and 12:
# these have many more DEGs than 21 so it ends up looking weird to show heatmap of 21 next to the other three

# Cluster 1 / Neutrophils 1
MFA1<-as.ggplot(pseudo_test_cl1$ploT_pca)+labs(title="A")
MFA2<-as.ggplot(pseudo_test_cl1$plot_UMAP)
MFB<-as.ggplot(pseudo_test_cl1$heatmap)+labs(title="B")
# Cluster 5 / Macrophages 1
MFC1<-as.ggplot(pseudo_test_cl5$ploT_pca)+labs(title="C")
MFC2<-as.ggplot(pseudo_test_cl5$plot_UMAP)
MFD<-as.ggplot(pseudo_test_cl5$heatmap)+labs(title="D")
# Cluster 12 / B-cells 3
MFE1<-as.ggplot(pseudo_test_cl12$ploT_pca)+labs(title="E")
MFE2<-as.ggplot(pseudo_test_cl12$plot_UMAP)
MFF<-as.ggplot(pseudo_test_cl12$heatmap)+labs(title="F")

# grid set up so that for each cluster we show the heatmap with the UMAP and PCA on top, and have three of these combinations in a row
grobz <- lapply(list(MFA1,MFA2,MFB,MFC1,MFC2,MFD,MFE1,MFE2,MFF), ggplotGrob)
lay <- rbind(c(1,1,1,2,2,4,4,4,5,5,7,7,7,8,8),
             c(3,3,3,3,3,6,6,6,6,6,9,9,9,9,9),
             c(3,3,3,3,3,6,6,6,6,6,9,9,9,9,9),
             c(3,3,3,3,3,6,6,6,6,6,9,9,9,9,9),
             c(3,3,3,3,3,6,6,6,6,6,9,9,9,9,9),
             c(3,3,3,3,3,6,6,6,6,6,9,9,9,9,9))
# write out figure PDF
pdf("Pseudobulk_DE/PseudobulkDE_main_composite_figure_171124.pdf", width = 21, height = 20)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()

# A smaller subsidiary figure shows the results of Cluster 21 / T-cells 5
MFG1<-as.ggplot(pseudo_test_cl21$ploT_pca)+labs(title="A")
MFG2<-as.ggplot(pseudo_test_cl21$plot_UMAP)
MFH<-as.ggplot(pseudo_test_cl21$heatmap)+labs(title="B")
# Slightly different layout: blank plots are used to create buffer space around the mini UMAP
blankplot<-ggplot()+theme_minimal()

grobz <- lapply(list(MFG2,MFG1,MFH,blankplot,blankplot), ggplotGrob)
lay <- rbind(c(4,4,1,1,5,5,3,3,3,3,3,3),
             c(2,2,2,2,2,2,3,3,3,3,3,3),
             c(2,2,2,2,2,2,3,3,3,3,3,3))
# write out figure PDF
pdf("Pseudobulk_DE/PseudobulkDE_subsid_composite_figure_171124.pdf", width = 10, height = 6)
grid.arrange(grobs = grobz, layout_matrix = lay)
dev.off()