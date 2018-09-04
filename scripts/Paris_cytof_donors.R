suppressPackageStartupMessages({
  library("openxlsx")
  library("FlowSOM")
  library("tidyverse")
  library("magrittr")
  library("flowCore")
  library("flowWorkspace")
  library("Rtsne")
})
library(BioGVHD)
options("scipen"=100)



################################################
############  aggregate flowframes #############
################################################

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
donor_names<-fcs_names[grep("D",names(fcs_names))]
donor_names <- donor_names[-which(names(donor_names)%in%c("12D","D1071","D369"))]

ff_agg_donor <- fcs_to_agg(fcs_dir= fcs_dir,
                           fcs_names= donor_names,
                           seed = 1,
                           cTotal = 10000*length(donor_names),
                           output_name = "aggregate_donor.fcs")

prettyMarkerNames <- ff_agg_donor@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_donor@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_donor)

save(ff_agg_donor, file = "ff_agg_donor.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/ff_agg_donor.RData")


## import metadat info :
dataPath <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx"
samp_donor <- import_patient_info(data_synthesis_file = dataPath,
                                  patient_names = donor_names,
                                  patient_type = "donor")

plot_aggregate_markers(patient_names = donor_names, samp_patients=samp_donor, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_date_34donor.png",
                       ff_agg = ff_agg_donor )
plot_aggregate_markers(patient_names = donor_names, samp_patients=samp_donor, color_by = "GROUP",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_group_donor_all_marks.png",
                       ff_agg = ff_agg_donor )

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/ff_agg_donor.RData")


######################################################
#######   Rescale values of 2 first patients   #######
######################################################

dim(ff_agg_donor)
files2rescale <- which(names(donor_names) %in% c("D1073", "D1502"))
ref_file <- which(names(donor_names) %in% c("D2031"))

min_ref <- apply(ff_agg_donor@exprs[which(ff_agg_donor@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                  function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_donor@exprs[which(ff_agg_donor@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                  function(x) quantile(x, 0.999))

for (file_nb in files2rescale){
  for (marker in colnames(ff_agg_donor@exprs)[c(3,17,28:62,71)]){
    ff_agg_donor@exprs[which(ff_agg_donor@exprs[,"File"]==file_nb), marker] <-
      scales::rescale(ff_agg_donor@exprs[which(ff_agg_donor@exprs[,"File"]==file_nb), marker],
                      to = c(min_ref[marker], max_ref[marker]))
  }
}

plot_aggregate_markers(patient_names = donor_names, samp_patients=samp_donor, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_rescaled_date_34donor.png",
                       ff_agg = ff_agg_donor )


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
colnames(patients_mfis) <- as.character(prettyMarkerNames[colnames(patients_mfis)])

ggplot_analysis_results("PCA", data_matrix = patients_mfis, metadata = samp_donor,
                        col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")
ggplot_analysis_results("tSNE", data_matrix = patients_mfis, metadata = samp_donor,
                        col_by = "GROUP", shape_by = "CMVStatus")





##################     FlowSOM      ###################
#######################################################

seed <- 1
set.seed(seed)
fsom <- FlowSOM(ff_agg_donor,
                colsToUse = colsToUse,
                scale = FALSE,
                xdim = 15, ydim = 15, # larger grid because ++ markers
                nClus = 40,
                seed = seed)
fsom$FlowSOM$prettyColnames <- prettyMarkerNames

PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          legend = FALSE)
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = c("Ce142Di","Nd144Di","Nd145Di","Nd146Di","Tm169Di","Er170Di","Yb172Di",
                      "Pr141Di","Sm147Di"))
PlotStars(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))])
save(fsom, file="fsom_donor.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/fsom_donor.RData")




####################################################################################
##### generate matrix of counts by matching all patient cells to fSOM clusters #####
####################################################################################

pctgs <- generate_donor_pctgs(
  donor_names = donor_names,
  fsom = fsom,
  pdf_name = "Plot_Stars_donors_32_marks.pdf",
  fcs_dir = fcs_dir,
  output_dir = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/annotated_fcs_files/"
)
#pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
pctgs_donor <- pctgs
save(pctgs_donor, file = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/pctgs_donor.RData")


#### Visualising cluster pctgs :
## PCA/ tSNE on cluster pctgs

rownames(pctgs_donor) <- names(rownames(pctgs_donor))
ggplot_analysis_results("tSNE", data_matrix = pctgs_donor, metadata = samp_donor,
                        col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")




####################################################################################
#######  generate meta_pctgs by matching patient cells to fSOM metaclusters  #######
####################################################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/fsom_donor.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/pctgs_donor.RData")

pctgs_meta_donor <- t(apply(pctgs_donor, 1, function(x){tapply(x, fsom$metaclustering, sum)}))
#rownames(pctgs_meta_donor) <- rownames(pctgs_donor)
save(pctgs_meta_donor, file="/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/pctgs_meta_donor.RData")

big_mat <- merge.data.frame(pctgs_meta_donor, samp_donor, by="row.names")
rownames(big_mat) <- big_mat$Row.names
big_mat[,2:41]<-apply(big_mat[,2:41],2,scale)

date_colors <- RColorBrewer::brewer.pal(7, "YlOrRd")
names(date_colors) <- names(table(big_mat[,51]))
big_mat[,51] <- as.factor(big_mat[,51])

pheatmap::pheatmap(as.matrix(big_mat[,2:41]),
                   #cluster_rows = hclust_meta, #(ward.D2 from ideas Yvan and Sofie, line 103)
                   cluster_rows = T,
                   scale = "none",
                   cluster_cols = T,
                   annotation_row = big_mat[,c(43,51)],
                   annotation_colors = list(
                     "GROUP"=c("non_tolerant"="#e31a1c90",
                             "primary_tolerant"="#00FF0090",
                             "secondary_tolerant"="#0000FF90"),
                     "DATEOFCYTOFEXPERIMENT"= date_colors),

                   # annotation_row = big_mat[,c(43,56)],
                   # annotation_colors = list("GROUP"=c("non_tolerant"="#e31a1c90",
                   #                                    "primary_tolerant"="#00FF0090",
                   #                                    "secondary_tolerant"="#0000FF90"),
                   #                          "aGVHD" = c("0" = "#bfd3e6",
                   #                                      "1" = "#88419d")),
                   #annotation_col = annot_cols,
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")



big_mat <- big_mat %>% arrange(GROUP)
rownames(big_mat) <- big_mat$Id.Cryostem.R

pheatmap::pheatmap(as.matrix(big_mat[which(big_mat$GROUP != "non_tolerant"),2:41]),
                   cluster_rows = T,
                   scale = "none",
                   cluster_cols = T,
                   annotation_row = big_mat[which(big_mat$GROUP != "non_tolerant"),c(43,51)],
                   labels_col = my_labels,
                   #col = colors,
                   show_rownames = TRUE,
                   cex=0.7,
                   show_colnames = T,
                   main = "Percentages")




#################################################################
#######  plotStar fsom metaclusters instead of clusters  ########
#################################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/fsom_donor.RData")
markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]

fsom_meta_donors <- fsom2fsom_meta(fsom = fsom, colsToUse = markersToPlot,
                                  pctgs_patients = pctgs_donor, plot_size = "equal_size")

PlotStars(fsom_meta_donors$FlowSOM,
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))])
PlotNumbers(fsom_meta_donors$FlowSOM)

#fsom_meta_recip$FlowSOM$MST$size <- rep(15, 30)
PlotStars(fsom_meta_donors$FlowSOM,
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))])#,
          #backgroundValues = my_labels,
          #backgroundColor = rainbow(n = 22, alpha = 0.3))

plot.new()
graphics::legend("center", legend = levels(as.factor(my_labels)),
                 fill = rainbow(22, alpha = 0.3), cex = 0.6, ncol = 2, bty = "n")

pdf("Plot_Stars_meta_recipients_21_marks.pdf")
for (i in 1:nrow(pctgs_recip)){
  fsom_meta_recip <- fsom2fsom_meta(fsom = fsom_recip, colsToUse = markersToPlot,
                                    pctgs_patients = pctgs_recip, plot_size = i)
  print(paste0("Plotting file ", rownames(pctgs_recip)[i]))
  PlotStars(fsom_meta_recip$FlowSOM,
            markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                             "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                             "foxP3","CD24","CXCR5"))],
            main = rownames(pctgs_recip)[i])
}
dev.off()


## RandomForest ---------------------------------------

library(randomForest)
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/pctgs_meta_donor.RData")

### On groups:
### on metadata only:
samp_donor_filtered <- samp_donor[,-c(1,4,6:9)]
rf_r<-randomForest(GROUP~., samp_donor_filtered, ntree=5000, mtry=3)
rf_r
tree_func(final_model = rf_r)
plot(rf_r) # big error (> 80%)

### on pctgs meta :
group <- as.data.frame(samp_donor$GROUP)
rownames(group) <- rownames(samp_donor)
colnames(group) <- "GROUP"

big_mat<- merge.data.frame(group, pctgs_meta_donor,
                           by = "row.names") %>%
  column_to_rownames("Row.names")

colnames(big_mat) <- c(colnames(big_mat)[1], paste0("meta",colnames(big_mat[,2:ncol(big_mat)])))
set.seed(1)
rf_meta<-randomForest(GROUP~., big_mat, ntree= 15000, mtry=20)
rf_meta
tree_func(final_model = rf_meta)
plot(rf_meta) # huge error (100% for primary and secondary)



### on groups tol 1 and tol 2, to see if there is any difference:
samp_tol <- samp_donor_filtered[which(samp_donor_filtered$GROUP != "non_tolerant"),]
samp_tol$GROUP <- as.factor(as.character(samp_tol$GROUP))

### on metadata only:
rf_r<-randomForest(GROUP~., samp_tol, ntree=10000, mtry=3)
rf_r
tree_func(final_model = rf_r)
plot(rf_r) # huge errors remaining (> 85%)

### on pctgs meta :
group <- as.data.frame(samp_tol$GROUP)
rownames(group) <- rownames(samp_tol)
colnames(group) <- "GROUP"

big_mat<- merge.data.frame(group, pctgs_meta_donor,
                           by = "row.names") %>%
  column_to_rownames("Row.names")

colnames(big_mat) <- c(colnames(big_mat)[1], paste0("meta",colnames(big_mat[,2:ncol(big_mat)])))
set.seed(1)
rf_meta<-randomForest(GROUP~., big_mat, ntree= 50000, mtry=20)
rf_meta
tree_func(final_model = rf_meta)
plot(rf_meta) # huge error (30% for secondary, 70% for primary)






