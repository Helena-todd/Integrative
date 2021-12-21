# Libraries and functions -----------------------------------------------------

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("flowCore"))
#suppressPackageStartupMessages(library("flowAI"))
suppressPackageStartupMessages(library("flowWorkspace"))
options("scipen"=100)
suppressPackageStartupMessages(library("MetaCyto"))
suppressPackageStartupMessages(library("FlowSOM"))

# PCA directly on MFIs -------------------------------------------------------

# rm CD19 and CD20, and CD45 (all positive)
markersToPlot<- markersToPlot[-1]
ff_BE_5pat<-ff_BE
ff_BE_5pat@exprs<-ff_BE@exprs[-which(exprs(ff_agg)[,72]%in%which(names(fcs_names)%in%c("12D","12R","D1071","D369","18R"))),]

library(data.table)
mfi<-ff_BE_5pat@exprs[,markersToPlot]
patients<-rep(samples$Id.Cryostem.R, each=10000)
mfis<-as.data.table(mfi)
mfis<-cbind(mfis,patients)

patients_mfis<-mfis[,lapply(.SD, median), by=patients]
patients_mfis<-as.data.frame(patients_mfis)
rownames(patients_mfis)<-patients_mfis[,1]
patients_mfis<-patients_mfis[,-1]

pca_mfi <- prcomp(patients_mfis)
tsne_mfi <- Rtsne(patients_mfis, perplexity = 8) 
samples_mfi <- samples %>% dplyr::mutate(pca_1 = pca_mfi$x[,1],
                                         pca_2 = pca_mfi$x[,2],
                                         tsne_1 = tsne_mfi$Y[,1],
                                         tsne_2 = tsne_mfi$Y[,2])

ggplot(samples_mfi) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

# FlowSOM on everything -------------------------------------------------------

#colsToUse <- c(1, 4, 7:12,15:18)
colsToUse<-markersToPlot[-1]
set.seed(seed)
fsom <- FlowSOM(ff_BE_5pat, 
                colsToUse = colsToUse,
                scale = FALSE,
                xdim = 15, ydim = 15, # larger grid because ++ markers
                nClus = 30,
                seed = seed)
fsom$FlowSOM$prettyColnames <- prettyMarkerNames
# fsom$FlowSOM <- BuildMST(fsom$FlowSOM, tSNE=TRUE)

PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          legend = FALSE)
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = c("Ce142Di","Nd144Di","Nd145Di","Nd146Di","Tm169Di","Er170Di","Yb172Di",
                      "Pr141Di","Sm147Di"))
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers=which(colnames(fsom$FlowSOM$data)%in%markersToPlot))
PlotStars(fsom$FlowSOM,backgroundValues = fsom$metaclustering,
          legend = FALSE)

## comparison between non tolerant and others:
groups <- tapply(fcs_names[-c(82:85,88)], samples$GROUP[-c(82:85,88)], 
                 function(x){x})
groups$tolerant <- c(groups$primary_tolerant,groups$secondary_tolerant) # group all tolerant together, flowSOM can compare only two groups
group2<-list(non_tolerant=counts[groups$non_tolerant,], 
             tolerant=counts[groups$tolerant,])
group_res <- CountGroups(fsom$FlowSOM, group2, plot = FALSE)

PlotGroups(fsom$FlowSOM,group_res,tresh = NULL,p_tresh = 0.05,
           markers = c("Ce142Di","Nd144Di","Nd145Di","Nd146Di","Tm169Di","Er170Di","Yb172Di",
                       "Pr141Di","Sm147Di"))
PlotGroups(fsom$FlowSOM,group_res,tresh = 2,p_tresh = NULL)

PlotNumbers(UpdateNodeSize(fsom$FlowSOM, reset = TRUE, maxNodeSize = 5))





boxplot(pctgs[gsub(".*/", "", groups$tolerant), 129], 
        pctgs[gsub(".*/", "", groups$non_tolerant),129])


wilcox.test(pctgs[gsub(".*/", "", groups$tolerant), 129], 
            pctgs[gsub(".*/", "", groups$non_tolerant),129])

## matrix of counts with fSOM clusters --------------------------

counts <- matrix(0,
                 length(fcs_names[-c(82:85,88)]),
                 ncol = fsom$FlowSOM$map$nNodes,
                 dimnames = list(fcs_names[-c(82:85,88)],
                 as.character(1:fsom$FlowSOM$map$nNodes)))
i <- 1
#pdf(file = "Plot_Stars_patients.pdf")
#for(i in seq_len(length(fcs_names))){
#  file <- fcs_names[i]
for (file in fcs_names[-c(82:85,88)]){
  message(file)
  ff <- read.FCS(file.path(fcs_dir,file))
  ff <- transform(ff,
            transformList(colnames(ff)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
  if(file %in% fcs_names[c("D1073","D1502")]){
    for (marker in markersToPlot){
      ff@exprs[,marker] <- (((ff@exprs[,marker] - min(ff@exprs[,marker])) / 
                                                                     (max(ff@exprs[,marker]) - min(ff@exprs[,marker]))) * (max_goal[marker] - min_goal[marker]) ) + min_goal[marker]
    } 
  } else {ff<-ff}
  fsom_tmp <- NewData(fsom$FlowSOM, ff)
  #PlotStars(fsom_tmp)

  # PlotPies(UpdateNodeSize(fsom_tmp, reset= TRUE, maxNodeSize = 8),
  #          backgroundValues = fsom$metaclustering,
  #          main = file)

  t <- table(fsom_tmp$map$mapping[,1])
  counts[file, names(t)] <- t
}
#dev.off()

pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))

save(fsom, counts, pctgs, file = "FlowSOM_88samples.Rdata")


# samples$EXP <- as.factor(samples$EXP)
# annotation_row <- samples[,c("Genotype","EXP")]
# rownames(annotation_row) <- paste0(names(workspaces)[samples$EXP], 
#                                    "_", 
#                                    samples$File)
# annotation_col <- data.frame(Metacluster = as.numeric(fsom$metaclustering))
# rownames(annotation_col) <- colnames(pctgs)
# 
# pheatmap::pheatmap(pctgs,
#                    scale = "column",
#                    annotation_row = annotation_row,
#                    annotation_col = annotation_col,
#                    show_rownames = FALSE,
#                    show_colnames = FALSE,
#                    main = "Percentages")

#rownames(pctgs)<-names(fcs_names)
pca <- prcomp(pctgs) 
samples_bis <- samples[-c(82:85,88),] %>% dplyr::mutate(pca_1 = pca$x[,1],
                                     pca_2 = pca$x[,2])
ggplot(samples_bis) + 
  geom_point(aes(x = pca_1, y = pca_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = pca_1, y = pca_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

library(Rtsne)
tsne<-Rtsne(pctgs,perplexity = 8)

samples_tsne <- samples[-c(82:85,88),] %>% dplyr::mutate(tsne_1 = tsne$Y[,1],
                                         tsne_2 = tsne$Y[,2])
ggplot(samples_tsne) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

pctgs_meta <- t(apply(pctgs, 1, function(x){tapply(x, fsom$metaclustering, sum)}))
rownames(pctgs_meta)<-samples[-c(82:85,88),]$Id.Cryostem.R
# annot<-as.data.frame(as.factor(samples$GROUP))
# rownames(annot)<-samples$Id.Cryostem.R
pheatmap::pheatmap(pctgs_meta,
                   #cluster_rows = FALSE,
                   scale = "column",
                   #annotation_row = annot,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = FALSE,
                   main = "Percentages")

pca_meta <- prcomp(pctgs_meta)$x
rownames(pca_meta) <- names(fcs_names)[-c(82:85,88)]
samples_ter <- samples[-c(82:85,88),] %>% 
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

samples_tsne <- samples[-c(82:85,88),] %>% dplyr::mutate(tsne_1 = tsne$Y[,1],
                                          tsne_2 = tsne$Y[,2])
ggplot(samples_tsne) + 
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(DATEOFCYTOFEXPERIMENT), shape = GROUP), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= Id.Cryostem.R)) +
  theme_minimal() +
  theme(text = element_text(size = 12))

# pca_manual <- prcomp(scale(apply(samples[,4:19],2,function(x){
#   as.numeric(gsub(" %", "", as.character(x)))})))
# samples <- samples %>% dplyr::mutate(pca_1_manual = pca_manual$x[,1],
#                                      pca_2_manual = pca_manual$x[,2])
# ggplot(samples) + 
#   geom_point(aes(x = pca_1_manual, y = pca_2_manual, col = EXP, shape = Genotype), size = 4) +
#   theme_minimal() +
#   theme(text = element_text(size = 12))

# -----------------------------------------------------------------------------
library(randomForest)
# I will analyse Receptors and Donors separately:
complete_data<-cbind(samples[-c(82:85,88),],pctgs)
useful_data<-complete_data[,-c(1,4:12,27:49)]
recip<-useful_data[grep("R",rownames(useful_data)),]
recip2<-data.frame(lapply(recip[,1:16], as.factor))
recip<-cbind(recip2,recip[,-c(1:16)])
colnames(recip)[17:241]<-paste0("cluster_",c(1:225))
rf_r<-randomForest(GROUP~., recip)

rf_r_metadata<-randomForest(GROUP~., recip2)

vector2col<-rep(1,225)
vector2col[c(13,63,83,55,12,150)]<-2
PlotNumbers(fsom$FlowSOM,backgroundValues = vector2col, backgroundColor = c("red","green"))


# -----------------------------------------------------------------------------

MetaCyto_dir <- "180211_MetaCyto"
dir.create(MetaCyto_dir)

fcs_info <- data.frame(fcs_files = file.path(preprocessed_dir,
                                             paste0(names(workspaces)[samples$EXP], 
                                                    "_", 
                                                    samples$File)),
                       study_id = samples$EXP)
write.csv(fcs_info, file.path(MetaCyto_dir, "fcs_info.csv"))

sample_info <- data.frame(fcs_files = file.path(preprocessed_dir,
                                                paste0(names(workspaces)[samples$EXP], 
                                                       "_", 
                                                       samples$File)),
                          Genotype = samples$Genotype)
write.csv(sample_info, file.path(MetaCyto_dir, "sample_info.csv"))

# Preprocessing - Can this be skipped? ----------------------------------------
#                 Creates an aggregate file per experiment
#                 Also compensates and transforms...
#                 Set to CyTOF to avoid compensation,
#                 and only transform scatters which shouldn't be used any way

preprocessing.batch(inputMeta = fcs_info,
                    assay = rep("CyTOF", nrow(fcs_info)),
                    b = rep(1/150, nrow(fcs_info)),
                    outpath = MetaCyto_dir,
                    excludeTransformParameters = colnames(ff)[-c(1:6)])

files=list.files(MetaCyto_dir,
                 pattern="processed_sample",
                 recursive=TRUE,
                 full.names=TRUE)
panel_info = collectData(files,longform=FALSE)

PS=panelSummary(panelInfo = panel_info,
                folder = MetaCyto_dir,
                cluster=FALSE,
                width=30,
                height=20)
knitr::kable(head(PS))
sort(rownames(PS))

# Clustering ------------------------------------------------------------------

excludeClusterParameters=c("FSC-A","FSC-W","FSC-H",
                           "SSC-A","SSC-W","SSC-H",
                           "TIME",
                           "VIABILITY#EFLUOR780",
                           "CD45#AF700",
                           "SAMPLE_ID")

cluster_label=autoCluster.batch(preprocessOutputFolder= MetaCyto_dir,
                                excludeClusterParameters = excludeClusterParameters,
                                labelQuantile = 0.95)

searchCluster.batch(preprocessOutputFolder = MetaCyto_dir,
                    outpath = MetaCyto_dir,
                    clusterLabel = cluster_label)

# Statistical analysis --------------------------------------------------------
pdf("MetaCytoTest.pdf",
    width = 14,
    height = 21)

library(dplyr)
files=list.files(MetaCyto_dir,
                 pattern="cluster_stats_in_each_sample",
                 recursive=TRUE,
                 full.names=TRUE)

fcs_stats=collectData(files,longform=TRUE)
sample_info=read.csv("180211_MetaCyto/sample_info.csv",
                     stringsAsFactors=FALSE,
                     check.names=FALSE)

all_data=inner_join(fcs_stats,sample_info,by="fcs_files")

GA=glmAnalysis(value="value",variableOfInterst="Genotype",parameter="fraction",
               otherVariables=c(),studyID="study_id",label="label",
               data=all_data,CILevel=0.95,ifScale=c(TRUE,FALSE))
GA=GA[order(GA$Effect_size),]

GA$label=as.character(GA$label)
GA$label <- gsub("#.*[A-Z0-9]([\\+\\-])$", "\\1", gsub("#[A-Z0-9\\-]+([\\+\\-][\\|])","\\1",GA$label))
w = which(nchar(GA$label)<30)
plotGA(GA[w,])
plotGA(GA[abs(GA$Effect_size) > 0.9,])

all_data$label <- gsub("#.*[A-Z0-9]([\\+\\-])$", "\\1", gsub("#[A-Z0-9\\-]+([\\+\\-][\\|])","\\1",all_data$label))

L="CD11C+|CD172A+|CD3CD19-|CD64+"
dat=subset(all_data,all_data$parameter_name=="fraction"&
             all_data$label==L)
MA=metaAnalysis(value="value",variableOfInterst="Genotype",main=L,
                otherVariables=c(),studyID="study_id",
                data=dat,CILevel=0.95,ifScale=c(TRUE,FALSE))

dat$Genotype <- factor(dat$Genotype,
                       levels = c( "WT", "SKO", "DKO", "TKO","DKO PERK", "DKO PERK HEZ"))

print(ggplot(dat,
             aes(x = Genotype, y = value)) +
        geom_boxplot() +
        theme_minimal() +
        ggtitle(L))

print(ggplot(dat,
             aes(x = Genotype, y = value, col=as.factor(study_id))) +
        geom_boxplot() +
        theme_minimal() +
        ggtitle(L) + facet_wrap("study_id"))

print(ggplot(dplyr::filter(fcs_stats, parameter_name == "fraction"),
             aes(x=fcs_names, y=label)) + 
        geom_tile(aes(fill = value), colour = "white") + 
        scale_fill_gradient(low = "white", high = "steelblue")+ 
        theme_grey(base_size = 9) + 
        labs(x = "", y = "") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)))


# opts(legend.position = "none",
#      axis.ticks = theme_blank(), 
#      axis.text.x = theme_text(size = 9 * 0.8, 
#                               angle = 330, 
#                               hjust = 0, 
#                               colour = "grey50")) 


frequencies <- spread(dplyr::filter(fcs_stats, parameter_name == "fraction") %>% 
                        dplyr::select(fcs_names, label, value),
                      label, value) 
rownames(frequencies) <- frequencies[,"fcs_names"]
frequencies <- frequencies[,-1]
colnames(frequencies) <- gsub("#[^|]*([+-])","\\1",colnames(frequencies))

annotation <- samples[,"Genotype", drop = FALSE]
rownames(annotation) <- gsub("/", "_", samples$filePath)

pheatmap::pheatmap(t(frequencies), scale = "none",
                   main = "MetaCyto Cluster frequencies",
                   annotation_col = annotation)

pheatmap::pheatmap(t(frequencies), scale = "row",
                   main = "MetaCyto Cluster frequencies - scaled",
                   annotation_col = annotation)


pheatmap::pheatmap(t(frequencies)[,order(samples$Genotype)], scale = "row",
                   main = "MetaCyto Cluster frequencies - scaled",
                   cluster_cols = FALSE,
                   gaps_col = which(sort(samples$Genotype) != 
                                      sort(samples$Genotype)[c(2:nrow(samples),nrow(samples))]),
                   annotation_col = annotation)




#### TMP ------------------
files <- file.path(preprocessed_dir,
                   paste0(names(workspaces)[samples$EXP], 
                          "_", 
                          samples$File))

counts <- matrix(0,
                 ncol = length(table(manual$manual)),
                 nrow = length(files),
                 dimnames = list(files, names(table(manual$manual))))

for(file in files){
  manual <- readRDS(gsub(".fcs","_manual.RDS",file))
  levels(manual$manual) <- gsub("M.*/","",levels(manual$manual))
  t <- table(manual$manual)
  counts[file, names(t)] <- t
}


## further analysis ----------------------------
## with CYTOF workflow ----------------------------

source("https://bioconductor.org/biocLite.R")
biocLite("cytofWorkflow")
