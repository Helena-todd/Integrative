suppressPackageStartupMessages({
  library("openxlsx")
  library("FlowSOM")
  library("tidyverse")
  library("magrittr")
  library("flowCore")
  library("flowWorkspace")
})
library(BioGVHD)
options("scipen"=100)

################################################
############  aggregate flowframes #############
################################################

# test

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
recip_names<-fcs_names[grep("R",names(fcs_names))]
recip_names<- recip_names[-which(names(recip_names)%in%c("12R","18R"))] # only CD19+ cells: removed

ff_agg_recip <- fcs_to_agg(fcs_dir= fcs_dir, fcs_names= recip_names, seed = 1, cTotal = 10000*49,
           output_name = "aggregate_recip.fcs")

prettyMarkerNames <- ff_agg_recip@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_recip@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_recip)

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/ff_agg_recip.RData")

plot_aggregate_markers(patient_names = recip_names, samp_patients=samp_recip, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, functional_marks, png_name= "Aggregate_date_recip_funct.png",
                       ff_agg = ff_agg_recip )
plot_aggregate_markers(patient_names = recip_names, samp_patients=samp_recip, color_by = "GROUP",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_group_recip_all_marks.png",
                       ff_agg = ff_agg_recip )



################## Markers to plot ###################
######################################################

markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
pheno_marks <- pheno_marks[-which(pheno_marks=="CD45")] # all cells are already pregated on CD45+
markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
names(pheno_marks) <- markersToPlot
colsToUse<-markersToPlot


################## PCA directly on MFIs #####################
#############################################################

library(data.table)
mfi<-ff_agg_recip@exprs[,markersToPlot]
patients<-rep(names(recip_names), each=10000)
mfis<-as.data.table(mfi)
mfis<-cbind(mfis,patients)

patients_mfis<-mfis[,lapply(.SD, median), by=patients]
patients_mfis<-as.data.frame(patients_mfis)
rownames(patients_mfis)<-patients_mfis[,1]
patients_mfis<-patients_mfis[,-1]
colnames(patients_mfis) <- as.character(prettyMarkerNames[colnames(patients_mfis)])

ggplot_analysis_results("PCA", data_matrix = patients_mfis, metadata = samp_recip,
                        col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")




#######################################################
##################     FlowSOM      ###################
#######################################################

seed <- 1
set.seed(seed)
fsom <- FlowSOM(ff_agg_recip,
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
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD27","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR5"))])

PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = c("Gd158Di","Yb174Di","Dy161Di","Eu151Di")) # marqueurs fonctionnels
save(fsom, file="fsom_recip.RData")



####################################################################################
##### generate matrix of counts by matching all patient cells to fSOM clusters #####
####################################################################################

pctgs <- generate_pctgs(
  recip_names = recip_names,
  fsom = fsom_recip,
  pdf_name = "Plot_Stars_recipients_32_marks.pdf",
  fcs_dir = fcs_dir,
  output_dir = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/annotated_fcs_files/"
)
#pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
save(fsom, counts, pctgs, file = "FlowSOM_49recipients.Rdata")
pctgs_recip <- pctgs
save(pctgs_recip, file = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/pctgs_recip.RData")

#### Visualising cluster pctgs :
## PCA/ tSNE on cluster pctgs

rownames(pctgs_recip) <- names(rownames(pctgs_recip))
ggplot_analysis_results("tSNE", data_matrix = pctgs_recip, metadata = samp_recip,
                        col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")



####################################################################################
#######  generate meta_pctgs by matching patient cells to fSOM metaclusters  #######
####################################################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/fsom_recip.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/pctgs_recip.RData")

pctgs_meta_recip <- t(apply(pctgs_recip, 1, function(x){tapply(x, fsom_recip$metaclustering, sum)}))
rownames(pctgs_meta_recip) <- names(rownames(pctgs_meta_recip))
save(pctgs_meta_recip, file="pctgs_meta_recip.RData")

big_mat <- merge.data.frame(pctgs_meta_recip, samp_recip, by="row.names")
rownames(big_mat) <- big_mat$Row.names
big_mat[,2:31]<-apply(big_mat[,2:31],2,scale)
big_mat[,57] <- as.factor(big_mat[,57])
big_mat[,43] <- as.factor(as.Date(big_mat[,43],format = "%d.%m.%Y"))
big_mat[which(big_mat$Id.Cryostem.R=="09R"),"GROUP"] <- "Secondary_tolerant"

date_colors <- RColorBrewer::brewer.pal(7, "YlOrRd")
names(date_colors) <- levels(big_mat[,43])

pheatmap::pheatmap(as.matrix(big_mat[,2:31]),
                   #cluster_rows = hclust_meta, #(ward.D2 from ideas Yvan and Sofie, line 103)
                   cluster_rows = F,
                   scale = "none",
                   cluster_cols = T,
                   # annotation_row = big_mat[,c(33,57,43)],
                   # annotation_colors = list(
                   #   "GROUP"=c("Non_Tolerant"="#e31a1c90",
                   #           "Primary_tolerant"="#00FF0090",
                   #           "Secondary_tolerant"="#0000FF90"),
                   #   "aGVHD" = c("0" = "#bfd3e6",
                   #               "1" = "#88419d"),
                   #   "DATEOFCYTOFEXPERIMENT"= date_colors),

                   annotation_row = big_mat[,c(33,57)],
                   annotation_colors = list("GROUP"=c("Non_Tolerant"="#e31a1c90",
                                                      "Primary_tolerant"="#00FF0090",
                                                      "Secondary_tolerant"="#0000FF90"),
                                            "aGVHD" = c("0" = "#bfd3e6",
                                                        "1" = "#88419d")),
                   #annotation_col = annot_cols,
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")

cbind(b)

#big_mat <- big_mat %>% arrange(GROUP)
#rownames(big_mat) <- big_mat$Id.Cryostem.R

pheatmap::pheatmap(as.matrix(big_mat[which(big_mat$GROUP != "Non_Tolerant"),2:31]),
                   cluster_rows = T,
                   scale = "none",
                   cluster_cols = T,
                   annotation_row = big_mat[which(big_mat$GROUP != "Non_Tolerant"),c(33,57)],
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")


ggplot_analysis_results("PCA", data_matrix = pctgs_meta_recip, metadata = samp_recip,
                        col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")
ggplot_analysis_results("tSNE", data_matrix = pctgs_meta_recip, metadata = samp_recip,
                        col_by = "CMVStatus", shape_by = "GROUP")



#################################################################
#######  plotStar fsom metaclusters instead of clusters  ########
#################################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/recip/fsom_recip.RData")
markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]

fsom_meta_recip <- fsom2fsom_meta(fsom = fsom_recip, colsToUse = markersToPlot,
                                  pctgs_patients = pctgs_recip, plot_size = 1)

PlotStars(fsom_meta_recip$FlowSOM,
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))])
PlotNumbers(fsom_meta_recip$FlowSOM)

fsom_meta_recip$FlowSOM$MST$size <- rep(15, 30)
PlotStars(fsom_meta_recip$FlowSOM,
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD27","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                         "CD38","CD161","CCR5"))],
          backgroundValues = my_labels,
          backgroundColor = c(rainbow(n = 7, alpha = 0.3), "#FFFFFF00"))

pdf("Plot_Stars_meta_recipients_9_marks.pdf")
for (i in 1:nrow(pctgs_recip)){
  fsom_meta_recip <- fsom2fsom_meta(fsom = fsom_recip, colsToUse = markersToPlot,
                                    pctgs_patients = pctgs_recip, plot_size = i)
  print(paste0("Plotting file ", names(rownames(pctgs_recip))[i]))
  PlotStars(fsom_meta_recip$FlowSOM,
            markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))],
            main = names(rownames(pctgs_recip))[i])
}
dev.off()




# Anova ------------------------------------------------------------------------

gr_res<- sample_recip$GROUP
gr_res[which((sample_recip$aGVHD==1)&(sample_recip$GROUP=="non_tolerant"))] <- "non_tol_GVHD"

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
gr<-rep(1, dim(ff_agg_recip@exprs)[1])
dat<-as.data.frame(cbind(ff_agg_recip@exprs, gr))
dat$gr<-rep(sample_recip$GROUP, each=10000)
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




###############################################################################
####### BOXPLOTS of marker expressions in different groups of patients ########
###############################################################################

## between tSNE clusters:
apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395","R212")
apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395")
#apart<- c("R1152","R2618","R2794","R709","R1131","R1267","R997","R773","09R","R874","R370","R297")
colnames(patients_mfis) <- as.character(prettyMarkerNames[colnames(patients_mfis)])

mark<-"CD45"
patients_mfis$clustered<- rep("other", dim(patients_mfis)[1])
patients_mfis$clustered[which(rownames(patients_mfis)%in%apart)]<- "clustered"

## between tolerant 1 and tolerant 2:
mfis_tol <- mfis %>% # remove non tolerant patients
  filter(!patients %in% samp_recip$Id.Cryostem.R[which(samp_recip$GROUP=="Non_Tolerant")])
colnames(mfis_tol)[-which(colnames(mfis_tol)=="patients")] <-
  as.character(prettyMarkerNames[which(names(prettyMarkerNames)%in%colnames(mfis_tol)[-which(colnames(mfis_tol)=="patients")])])

mfis_tol$clustered <- rep("Tolerant_1", nrow(mfis_tol))
tol_2 <- samp_recip$Id.Cryostem.R[which(samp_recip$GROUP=="Secondary_tolerant")]
mfis_tol$clustered[which(mfis_tol$patients%in%tol_2)]<- "Tolerant_2"

mat2plot <- mfis_tol
png("boxplots_tol1_tol2.png",
    width = 4000,
    height = 2000)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:30, nrow=5, byrow = TRUE))
lapply(seq_along(colnames(mat2plot)[-c(ncol(mat2plot),ncol(mat2plot)-1)]), function(i){
  boxplot(mat2plot[,i]~clustered, data=mat2plot,
          main=colnames(mat2plot)[i], xaxt="n", cex.main=4,
          cex.axis=3, cex.sub=3, col=c("lightgreen","blue"))
  axis(side=1, at=c(1:2), labels=c("Tolerant_1","Tolerant_2"), las=2,cex=1)
})
dev.off()




## RandomForest ---------------------------------------

library(randomForest)
# I will analyse Recipients separately:

samp_recip_filtered <- import_patient_info(data_synthesis_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx",
                                                      patient_names = recip_names)
samp_recip_filtered <- samp_recip_filtered[,-c(1,4,7,9,10,27)]
apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395")

# on metadata only, to understand weird group:
status <- rep ("normal", nrow(samp_recip_filtered))
status[which(rownames(samp_recip_filtered)%in%apart)] <- "strange"
status <- as.factor(status)
samp_recip_annot <- samp_recip_filtered %>%
  mutate(status = status)
rownames(samp_recip_annot) <- rownames(samp_recip_filtered)

rf_r<-randomForest(status~., samp_recip_annot, ntree=15000, mtry=20)
tree_func(final_model = rf_r)
rf_r
plot(rf_r)

# on pctgs only, to understand weird group :

annot_status <- as.data.frame(samp_recip_annot$status)
rownames(annot_status) <- rownames(samp_recip_annot)
colnames(annot_status) <- "status"

big_mat<- merge.data.frame(annot_status, pctgs_meta_recip,
                           by = "row.names") %>%
  column_to_rownames("Row.names")

colnames(big_mat) <- c(colnames(big_mat)[1], paste0("meta",colnames(big_mat[,2:ncol(big_mat)])))
set.seed(1)
rf_meta<-randomForest(status~., big_mat, ntree= 5000, mtry=20)
tree_func(final_model = rf_meta)
rf_meta
plot(rf_meta)

pheatmap::pheatmap(pctgs_meta_recip, annotation_row = annot_status,
                   cex=.8)


#on mfis
recip_data <- merge.data.frame(annot_status, patients_mfis, by = "row.names") %>%
  column_to_rownames("Row.names")

set.seed(1)
rf_mfis<-randomForest(status~., recip_data, ntree= 5000, mtry=20)
tree_func(final_model = rf_mfis)
rf_mfis
plot(rf_mfis)

pheatmap::pheatmap(recip_data[,-1], annotation_row = annot_status,
                   cex=.8)

vector2col<-rep(1,225)
vector2col[c(15,121,182,157,24,219,12)]<-2
PlotNumbers(fsom$FlowSOM,backgroundValues = vector2col, backgroundColor = c("red","green"))



