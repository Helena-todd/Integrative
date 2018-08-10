library(FlowSOM)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(openxlsx)

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/fsom_recip.Rdata")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/fsom_meta.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/ff_agg_recip.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/pctgs_and_MFIs.RData")

#####################################################
#########   find metaclusters cell types   ##########
#####################################################

PlotStars(fsom_meta$FlowSOM)

## table of associations btw markers and populations ------
pop_mark <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_3:04:18/CYTOF_David_Michonneau_excel/Populations and markers_filtered.xlsx",
                check.names=F) %>%
  dplyr::filter(!is.na(Markers)) %>%
  select(-c(CD45, GranzymeB)) # rm CD45: they are all positive (pre-gated by laetitia)
                              # rm granzymeB: functional marker, not phenotypic

# extract info of which markers are expressed in which cell types:
cellTypes <- markers_of_cellTypes(marker_table = pop_mark)

# identify cell type for each fsom cluster
celllabels <- identify_fsom_cellPopulations(fsom = fsom_meta, prettyMarkerNames, cellTypes,
                              pdf_name = "pops_filt_marks_grid.pdf", view="MST")


## QueryStarPlot doesn't work so well -> I set the labels manually
PlotStars(fsom_meta$FlowSOM, markers = c("CD11a","CD16","CD27","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                         "CD38","CD161","CCR5"))

my_labels <- c("Monocytes","mDCs","Monocytes","Unknown","B","B","B","B","B","Plasmo","CD8","CD4","CD4","Unknown","B",
               "Unknown","CD4","Monocytes","Unknown","CD8","CD4","CD8","CD8","CD4","CD4","CD4",
               "CD4","CD8","CD8","CD4-CD8 TSCM")

my_labels2 <- c("NK?","mDC?","Monocytes","CD25","B memory?","B naive","B naive","B naive","B",
                "Plasmoblasts","CD8 Temra","CD4 TSCM","TH1","Myeloid?","B naive","CD3 CD45RA CD11a","CD4 TEM",
                "Monocytes","NKT? Treg?","TFH reg?","CD4 TEM","CD8 TEM","CD8 TEM","CD4 TCM?","CD4 naive",
                "TH17","TH2","CD8 TCM","CD8 TCM","CD4 and CD8 TSCM")

PlotStars(fsom_meta$FlowSOM,
          markers=c("CD11a","CD16","CD27","CD3","CD4","CD45RA","CD8a","HLADR","CD19","CD38","CD161","CCR5"),
          backgroundValues = my_labels,
          backgroundColor = c(rainbow(n = 7, alpha = 0.2), "#FFFFFF00"))

plot.new()
graphics::legend("center", legend = levels(as.factor(my_labels)),
                 fill = rainbow(11, alpha = 0.3), cex = 0.7, ncol = 2, bty = "n")




################################################################################
#############    clustering on pctgs_meta + functional markers    ##############
################################################################################

## On pctgs_meta only :

big_mat <- merge.data.frame(pctgs_meta_recip, samp_recip, by="row.names")
rownames(big_mat) <- big_mat$Row.names
dist_mat <- dist(big_mat[,2:31])
library(dendextend)
hclust_meta<-hclust(dist_mat, method = "ward.D2")
dend <- as.dendrogram(hclust_meta)
grop<-rep(1,49)
grop[which(big_mat$GROUP=="Secondary_tolerant")]<-2
grop[which(big_mat$GROUP=="Non_Tolerant")]<-3
labels_colors(dend) <-
  c("red","blue","green")[sort_levels_values(
    as.numeric(grop)[order.dendrogram(dend)]
  )]
dend <- set(dend, "labels_cex", 0.7)

plot(dend,
     main = "Clustering on metaclusters only",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = c("non_tolerant","primary_tol","secondary_tol"), fill = c("red","green","blue"),cex=0.75)



# On pctgs_meta and functional markers :

metacluster_ids <- FlowSOM::GetMetaclusters(fsom_recip)
# functional_marks<-markers[which(markers[,5]==1),1]
# functional_marks <- functional_marks[-c(4,6,9)] # Lag3, PD1 and CD24 already used as pheno_marks
# names(functional_marks) <- names(prettyMarkerNames)[which(prettyMarkerNames%in%functional_marks)]
data("functional_marks")
markersToCluster <- names(functional_marks)

# heterogeneous expression of funct markers in metaclusters?
clusters <- lapply(seq_along(fsom$metaclustering), function(i){
  clust <- fsom$FlowSOM$data[which(fsom$FlowSOM$map$mapping[,1]==i),]
})

metaclusters <- lapply(seq_along(1:length(table(fsom$metaclustering))), function(i){
  meta_clust <- do.call("rbind",clusters[which(fsom$metaclustering==i)])
  meta_clust
})

lapply(seq_along(metaclusters), function(j){
  png(paste0("funct_markers_heterogeneity_cluster_",j,".png"),
      width = 4000,
      height = 2000)
  par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
  layout(matrix(1:9, nrow=3, byrow = TRUE))
  lapply(seq_along(functional_marks), function(i){
    hist(metaclusters[[j]][,names(prettyMarkerNames[which(prettyMarkerNames==functional_marks[i])])],
         breaks=1000,
         main = functional_marks[i],
         cex.main=3,
         xlab="")
  })
  dev.off()
})



pctgs_and_MFIs <- compute_marker_MFIs_per_metacluster(ff_agg = ff_agg_recip, pctgs_meta = pctgs_meta_recip,
                                                      fsom = fsom_recip, functional_marks)

# je merge pctgs_and_MFIs avec sample_recip
big_mat <- merge(pctgs_and_MFIs, samp_recip, by = "row.names")
rownames(big_mat) <- big_mat$Row.names
# rm columns containing NAs :
big_mat <- dplyr::select(big_mat,which(colSums(is.na(big_mat))==0))

pheatmap(big_mat[,2:241])
annot = as.data.frame(big_mat$GROUP)
rownames(annot) <- rownames(big_mat)
pheatmap::pheatmap(as.matrix(big_mat[,2:241]), annotation_row = annot, cex=.7)

big_mat[,2:241] <- apply(big_mat[,2:241],2,scale) # scale matrix
big_mat[,2:31] <- big_mat[,2:31 * length(functional_marks)] # rebalance weights so that weights pctgs = weights MFIs
dist_mat <- dist(big_mat[,2:241])
library(dendextend)
hclust_meta<-hclust(dist_mat, method = "ward.D2")
dend <- as.dendrogram(hclust_meta)
grop<-rep(1,49)
grop[which(big_mat$GROUP=="Secondary_tolerant")]<-2
grop[which(big_mat$GROUP=="Non_Tolerant")]<-3
labels_colors(dend) <-
  c("blue","red","green")[sort_levels_values(
    as.numeric(grop)[order.dendrogram(dend)]
  )]
dend <- set(dend, "labels_cex", 0.7)

plot(dend,
     main = "Clustering on fsom metaclusters and functional marker MFIs",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = c("non_tolerant","primary_tol","secondary_tol"), fill = c("red","green","blue"),cex=0.75)

pheatmap::pheatmap(as.matrix(big_mat[,2:241]), cluster_rows = hclust_meta, annotation_row = annot, cex=.7)
pheatmap::pheatmap(as.matrix(big_mat[,2:241]),
                   cluster_rows = hclust_meta,
                   annotation_row = big_mat[,c("GROUP","DATEOFCYTOFEXPERIMENT")],
                   cex=.7)


## clustering only on metadata :
filtered_metadata <- samp_recip[,c(2,3,12:24,26)]
bla <- filtered_metadata
for ( i in 1:16){
  bla[,i] <- as.numeric(as.factor((bla[,i])))
}



####################################
####### concensus clustering #######
####################################

library(ConsensusClusterPlus)
title=tempdir()
pdf("concensus_clust.pdf")
set.seed(seed = 1)
results = ConsensusClusterPlus(t(big_mat[,2:241]),maxK=10,reps=200,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1)
dev.off()
save(results, file = "results_consensus_clustering.RData")

# load consensus clustering results :
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/results_consensus_clustering.RData")

heatmap_k7 <- results[[7]]$consensusMatrix
rownames(heatmap_k7) <- rownames(big_mat)
colnames(heatmap_k7) <- rownames(big_mat)
pheatmap::pheatmap( heatmap_k7, annotation_row = big_mat[,c(243,244,266)],
                    fontsize_row = 5, fontsize_col = 5)
pheatmap::pheatmap( heatmap_k7, annotation_row = big_mat[,c(252,243)],
                    fontsize_row = 5, fontsize_col = 5)
pheatmap::pheatmap( heatmap_k7, annotation_row = big_mat[,c(266:267)])

res_dend <- results[[7]][["consensusTree"]]
res_dend$labels <- rownames(big_mat)
dend <- as.dendrogram(res_dend)
plot(dend,
     main = "Clustering on fsom metaclusters and functional marker MFIs",
     horiz =  F,  nodePar = list(cex = .007))

pheatmap(results[[7]]$ml)



#########################################################
#####  clustering on 3 dimensions at the same time  #####
#########################################################

## how do I find the info about my patients afterwards?

library(reshape2)   # for melt(...)
library(rgl)        # for plot3d(...)

set.seed(1)         # to create reproducible sample

# 3D matrix, values clustered around -2 and +2
m      <- c(rnorm(500,-2),rnorm(500,+2))
dim(m) <- c(10,10,10)
v      <- melt(m, varnames=c("x","y","z"))  # 4 columns: x, y, z, value
# interactive 3D plot, coloring based on value
plot3d(v$x,v$y,v$z, col=1+round(v$value-min(v$value)),size=5)
# identify clusters
v      <- scale(v)                          # need to scale or clustering will fail
v      <- data.frame(v)                     # need data frame for later
d  <- dist(v)                               # distance matrix
km <- kmeans(d,centers=2)                   # kmeans clustering, 2 clusters
v$clust <- km$cluster                       # identify clusters
# plot the clusters
plot(v[1:4],col=v$clust)                    # scatterplot matrix
plot3d(v$x,v$y,v$z, col=v$clust,size=5)




