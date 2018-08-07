suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("FlowSOM"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("flowCore"))
suppressPackageStartupMessages(library("flowWorkspace"))
suppressPackageStartupMessages(library("Rtsne"))
options("scipen"=100)

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/ff_agg_donor")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/sample_donor.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/donor_names.RData")

################## Markers to plot ###################
######################################################

markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx", 
                     check.names = FALSE) ## excel file with info about markers
prettyMarkerNames <- ff_agg_donor@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <- 
  ff_agg_donor@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_donor@exprs)

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
markersToPlot<- markersToPlot[-1]
colsToUse<-markersToPlot






################## PCA directly on MFIs #####################
#############################################################

library(data.table)
mfi<-ff_agg_donor@exprs[,markersToPlot]
patients<-rep(names(donor_names), each=10000)
mfis<-as.data.table(mfi)
mfis<-cbind(mfis,patients)

patients_mfis<-mfis[,lapply(.SD, median), by=patients]
patients_mfis<-as.data.frame(patients_mfis)
rownames(patients_mfis)<-patients_mfis[,1]
patients_mfis<-patients_mfis[,-1]

pca_mfi <- prcomp(patients_mfis)
tsne_mfi <- Rtsne(patients_mfis, perplexity = 8) 
samples_mfi <- sample_donor %>% dplyr::mutate(pca_1 = pca_mfi$x[,1],
                                              pca_2 = pca_mfi$x[,2],
                                              tsne_1 = tsne_mfi$Y[,1],
                                              tsne_2 = tsne_mfi$Y[,2])

ggplot(samples_mfi) + 
  geom_point(aes(x = pca_1, y = pca_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = pca_1, y = pca_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

ggplot(samples_mfi) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(GROUP)), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))






##################     FlowSOM      ###################
#######################################################

seed <- 1
set.seed(seed)
fsom <- FlowSOM(ff_agg_donor, 
                colsToUse = colsToUse,
                scale = FALSE,
                xdim = 15, ydim = 15, # larger grid because ++ markers
                nClus = 30,
                seed = seed)
fsom$FlowSOM$prettyColnames <- prettyMarkerNames

PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          legend = FALSE)
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = c("Ce142Di","Nd144Di","Nd145Di","Nd146Di","Tm169Di","Er170Di","Yb172Di",
                      "Pr141Di","Sm147Di"))
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = c("Gd158Di","Yb174Di","Dy161Di","Eu151Di")) # marqueurs fonctionnels
save(fsom, file="fsom_donor.RData")




## matrix of counts with fSOM clusters --------------------------

counts <- matrix(0,
                 length(donor_names),
                 ncol = fsom$FlowSOM$map$nNodes,
                 dimnames = list(donor_names,
                                 as.character(1:fsom$FlowSOM$map$nNodes)))
i <- 1
pdf(file = "Plot_Stars_donors_32_marks.pdf")
for (i in seq_along(donor_names)){
  file<-donor_names[[i]]
  message(file)
  ff <- read.FCS(file.path(fcs_dir,file))
  ff <- transform(ff,
                  transformList(colnames(ff_agg_donor)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
  fsom_tmp <- NewData(fsom$FlowSOM, ff)
  name<-names(donor_names)[i]
  PlotStars(fsom_tmp,main = name)
  
  #PlotPies(UpdateNodeSize(fsom_tmp, reset= TRUE, maxNodeSize = 8),
  #          backgroundValues = fsom$metaclustering,
  #          main = file)
  
  t <- table(fsom_tmp$map$mapping[,1])
  counts[file, names(t)] <- t
}
dev.off()

pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
save(fsom, counts, pctgs, file = "FlowSOM_49donors.Rdata")


## visualising cluster percentages
library(pheatmap)
toPlot_2 <- pctgs[,c(24,39)]
rownames(toPlot_2) <- rownames(sample_donor)
pheatmap(toPlot_2[myord,] ,
         annotation_row = sample_donor[,"GROUP", drop = FALSE],
         cluster_rows = FALSE
)






## PCA/ tSNE on cluster pctgs --------------------------------------

pca <- prcomp(pctgs) 
samp_donor <- sample_donor %>%
  slice(match(names(donor_names), rownames(sample_donor)))
samp_bis <- samp_donor %>% dplyr::mutate(pca_1 = pca$x[,1],
                                         pca_2 = pca$x[,2])
ggplot(samp_bis) + 
  geom_point(aes(x = pca_1, y = pca_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = pca_1, y = pca_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

library(Rtsne)
tsne<-Rtsne(pctgs,perplexity = 8)

samp_donor <- samp_donor %>% dplyr::mutate(tsne_1 = tsne$Y[,1],
                                           tsne_2 = tsne$Y[,2])
ggplot(samp_donor) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

pctgs_meta <- t(apply(pctgs, 1, function(x){tapply(x, fsom$metaclustering, sum)}))
rownames(pctgs_meta)<-sample_donor$Id.Cryostem.R
save(pctgs_meta, file="pctgs_meta.RData")


annot_cols <- as.data.frame(meta_size)



bla<-apply(pctgs_meta,2,scale)
rownames(bla)<-rownames(pctgs_meta)
bla <- as.data.frame(bla) %>%
  mutate(id = rownames(pctgs_meta),
         group = sample_donor$GROUP, 
         aGVHD = sample_donor$aGVHD,
         day = sample_donor$DATEOFCYTOFEXPERIMENT) %>%  
  arrange(group) %>% 
  column_to_rownames('id') 
bla[,33] <- as.factor(bla[,33])

pheatmap::pheatmap(as.matrix(bla[,-c(31,32,33)]),
                   cluster_rows = T,
                   scale = "none",
                   cluster_cols = F,
                   annotation_row = bla[,31:32],
                   #annotation_col = annot_cols,
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")

pheatmap::pheatmap(as.matrix(bla[which(bla$group != "non_tolerant"),-c(31,32)]),
                   cluster_rows = T,
                   scale = "none",
                   cluster_cols = T,
                   annotation_row = bla[which(bla$group != "non_tolerant"),31:32],
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")




#colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(255)
#pheatmap(my_matrix,cluster_cols = FALSE,cellwidth = 30,fontsize = 7,height = 40,show_rownames = FALSE, col=colors)


spca_meta <- prcomp(pctgs_meta)$x
rownames(pca_meta) <- names(donor_names)
samples_ter <- sample_donor %>% 
  #dplyr::arrange(DATEOFCYTOFEXPERIMENT) %>% 
  dplyr::mutate(pca_1_meta = pca_meta[, 1],
                pca_2_meta = pca_meta[, 2])

ggplot(samples_ter) + 
  #geom_point(aes(x = pca_1_meta, y = pca_2_meta), size = 4) +
  geom_point(aes(x = pca_1_meta, y = pca_2_meta, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  #geom_text(aes(x = pca_1_meta, y = pca_2_meta, label = seq_along(fcs_names))) +
  geom_text(aes(x = pca_1_meta, y = pca_2_meta,label=Id.Cryostem.R))+
  theme_minimal() +
  theme(text = element_text(size = 12))

tsne<-Rtsne(pctgs_meta,perplexity = 8)
samples_ter <- sample_donor %>% dplyr::mutate(tsne_1 = tsne$Y[,1],
                                              tsne_2 = tsne$Y[,2])
ggplot(samples_ter) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))


## plotStars metaclusters instead of clusters -----------------------

metacluster_MFIs <- t(apply(fsom$FlowSOM$data[, colsToUse], 2, function(x){
  tapply(x, fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]], median)
}))

metacluster_pctgs <- t(apply(pctgs, 1, function(x){ 
  tapply(x, fsom$metaclustering, sum) 
}))

fsom_test <- fsom
fsom_test$FlowSOM$map$codes <- t(metacluster_MFIs)
fsom_test$FlowSOM$map$medianValues <- t(metacluster_MFIs)
fsom_test$FlowSOM$map$colsUsed <- seq_along(rownames(metacluster_MFIs))
fsom_test$FlowSOM <- BuildMST(fsom_test$FlowSOM)

# Size depending on 1 file
fsom_test$FlowSOM$MST$size <- metacluster_pctgs[1,]*100 
# Size equal for all nodes
fsom_test$FlowSOM$MST$size <- rep(15, ncol(metacluster_MFIs))
fsom_test$FlowSOM$prettyColnames <- fsom_test$FlowSOM$prettyColnames[fsom$FlowSOM$map$colsUsed]

PlotStars(fsom_test$FlowSOM, 
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))])
PlotNumbers(fsom_test$FlowSOM)



# Anova ------------------------------------------------------------------------

gr_res<- sample_donor$GROUP
gr_res[which((sample_donor$aGVHD==1)&(sample_donor$GROUP=="non_tolerant"))] <- "non_tol_GVHD"

p_v <- rep(NA, ncol(pctgs))
for (i in seq_len(ncol(pctgs))) {
  data_tmp <- data.frame(Var = pctgs[, i],
                         #Visit = group_res$groups)
                         Visit = group_res$groups)
  
  library(ggplot2)
  ggplot(data_tmp) + geom_boxplot(aes(x = Visit, y = Var)) + theme_minimal()
  
  fit <- aov(Var ~ Visit, data = data_tmp)
  p_v[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
}

means <- apply(pctgs, 2, function(x){tapply(x, group_res$groups, mean)})
means_norm <- (means - min(means))/(max(means) - min(means))

type <- factor(c("--", "Decreased", "Increased")[1 + 
                                                   (p_v < 0.05) + 
                                                   (p_v < 0.05 & means[3,] > means[1,])],
               levels = c("--", "Decreased", "Increased"))

plot_list = list()
for (i in seq_along(which(type!="--"))){
  data_tmp <- data.frame(Var = pctgs[, which(type!="--")[i]],
                         Visit = group_res$groups)
  p=ggplot(data_tmp) + geom_boxplot(aes(x = Visit, y = Var)) + theme_minimal() +
    ggtitle(paste0("Cluster ", which(type!="--")[i]))
  plot_list[[i]] = p
}

pdf(file="differing_clusters.pdf")
for (i in seq_along(which(type!="--"))) {
  print(plot_list[[i]])
}
dev.off()

PlotStars(UpdateNodeSize(fsom$FlowSOM, reset = TRUE, maxNodeSize = 5),
          markers=c("Ce142Di","Nd144Di","Nd145Di","Nd146Di","Tm169Di","Er170Di","Yb172Di",
                    "Pr141Di","Sm147Di"),
          backgroundValue = type,
          backgroundColor = c("#FFFFFF00", "#00FFFF55", "#FF000055"))
PlotNumbers(UpdateNodeSize(fsom$FlowSOM, reset = TRUE, maxNodeSize = 0),
            backgroundValue = type,
            backgroundColor = c("#FFFFFF00", "#00FFFF55", "#FF000055")) # change font size!

PlotVariable(fsom$FlowSOM, as.numeric(fsom$metaclustering))


## plot each cluster of interest separately:
gr<-rep(1, dim(ff_agg_donor@exprs)[1])
dat<-as.data.frame(cbind(ff_agg_donor@exprs, gr))
dat$gr<-rep(sample_donor$GROUP, each=10000)
colnames(dat)[1:73]<-prettyMarkerNames
dat<-dat[which(fsom$FlowSOM$map$mapping[,1]==24),]


pdf("Cluster_24.pdf", width = 20, height = 20)
plots <- list()
for(i in seq_along(prettyMarkerNames[markersToPlot])){
  marker<-prettyMarkerNames[markersToPlot][i]
  plots[[i]] <- ggplot(dat, aes_string(x = as.character(marker), fill = "gr")) + geom_density(alpha = 0.5)
  #tidy_df <- dat %>% gather(col, val, one_of(prettyMarkerNames[markersToPlot])))
  #ggplot(tidy_df) + geom_point(aes(x, val)) + facet_wrap(~col)
}
multiplot(plotlist = plots, 
          layout = matrix(1:36, nrow = 6, byrow = TRUE)[6:1,])
dev.off()





## comparing marker expressions in clusters tSNE ----------------------

apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395","R212")
#apart<- c("R1152","R2618","R2794","R709","R1131","R1267","R997","R773","09R","R874","R370","R297")
colnames(patients_mfis) <- as.character(prettyMarkerNames[colnames(patients_mfis)])

mark<-"CD45"
patients_mfis$clustered<- rep("other", dim(patients_mfis)[1])
patients_mfis$clustered[which(rownames(patients_mfis)%in%apart)]<- "clustered"


png("boxplots.png",
    width = 4000,
    height = 2000)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:30, nrow=5, byrow = TRUE))
lapply(seq_along(colnames(patients_mfis)[-27]), function(i){
  boxplot(patients_mfis[,i]~clustered, data=patients_mfis, 
          main=colnames(patients_mfis)[i], xaxt="n", cex.main=4,
          cex.axis=3, cex.sub=3, col=c("lightgreen","gray"))
  axis(side=1, at=c(1:2), labels=c("clustered","other"), las=2,cex=1)
})
dev.off()




## RandomForest ---------------------------------------

library(randomForest)
# I will analyse Receptors separately:
donor_data<-cbind(sample_donor[,-c(1,4:12,27:49)],pctgs)
donor2<-data.frame(lapply(donor_data[,1:16], as.factor))
donor<-cbind(donor2,donor_data[,-c(1:16)])
colnames(donor)[17:241]<-paste0("cluster_",c(1:225))
set.seed(seed)
rf_clust<-randomForest(GROUP~., donor, ntree= 5000, mtry=100)

set.seed(seed)
rf_r_metadata<-randomForest(GROUP~., donor2)

#on metaclusters and metadata
donor_meta<-cbind(sample_donor[,-c(1,4:12,27:49)],pctgs_meta)
donor2<-data.frame(lapply(donor_meta[,1:16], as.factor))
donor<-cbind(donor2,donor_meta[,-c(1:16)])
colnames(donor)[17:46]<-paste0("cluster_",c(1:30))
rf_r<-randomForest(GROUP~., donor)

#on metaclusters only
donor_meta<-as.data.frame(cbind(sample_donor[,2],pctgs_meta))
donor_meta[,1]<-as.factor(sample_donor[,2])
colnames(donor_meta)<-c("GROUP",paste0("cluster_",c(1:30)))
set.seed(seed)
rf_r<-randomForest(GROUP~., donor_meta, ntree=5000, mtry=25)

#on mfis
donor_data<-cbind(sample_donor[,-c(1,4:12,27:49)],patients_mfis[,-27])
donor2<-data.frame(lapply(donor_data[,1:16], as.factor))
donor<-cbind(donor2,donor_data[,-c(1:16)])
colnames(donor)[17:42]<-colnames(patients_mfis)[-27]
set.seed(seed)
rf_mfis<-randomForest(GROUP~., donor, ntree= 50000, mtry=20)

vector2col<-rep(1,225)
vector2col[c(15,121,182,157,24,219,12)]<-2
PlotNumbers(fsom$FlowSOM,backgroundValues = vector2col, backgroundColor = c("red","green"))



