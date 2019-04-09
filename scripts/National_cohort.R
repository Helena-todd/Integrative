suppressPackageStartupMessages({
  library("openxlsx")
  library("FlowSOM")
  library("tidyverse")
  library("magrittr")
  library("flowCore")
  library("flowWorkspace")
  library("ggraph")
  library("igraph")
  library("scales")
  library("dendextend")
  library("cowplot")
})
library(BioGVHD)
options("scipen"=100)


## Import cyto data from the national cohort :

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/National_cohort/CYTOF/CYTOF_David_Michonneau_fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)

## Import associated metadata :

samp_rd_all <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/National_cohort/Data synthesis national cohort Cryostem 22222019.xlsx")
samp_rd_all$DATEOFCYTOFEXPERIMENT <- as.Date(samp_rd_all$DATEOFCYTOFEXPERIMENT, format = "%d.%m.%Y")

## Filter out data :
samp_rd <- samp_rd_all %>%
  dplyr::filter(DELAY_SAMPLE >= 148) %>%
  dplyr::filter(!COUPLENUMBER %in% c(32, 71)) %>%
  arrange(DATEOFCYTOFEXPERIMENT)#arrange(HOSPITAL,DATEOFCYTOFEXPERIMENT)

save(samp_rd, file = "~/Documents/VIB/Projects/Integrative_Paris/National_cohort/samp_rd_national.RData")

rd_names <- fcs_names[samp_rd$Id.Cryostem.R] # rearrange per hospital for the plots


## plot aggregate :

ff_agg_rd <- fcs_to_agg(fcs_dir= fcs_dir,
                       fcs_names= rd_names,
                       seed = 1,
                       cTotal = 10000*length(rd_names),
                       output_name = "aggregate_rd_national_per_date.fcs")

prettyMarkerNames <- ff_agg_rd@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_rd@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_rd)

x_names <- paste0(samp_rd$Id.Cryostem.R, "_", samp_rd$DATEOFCYTOFEXPERIMENT)
plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_nat_date_140rd.png",
                       x_names = samp_rd$Id.Cryostem.R ,ff_agg = ff_agg_rd )

save(ff_agg_rd,
     file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_140_date.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_140RData.RData")


#####################################
##########   SHIFT DOWN   ###########
#####################################

ff_agg_rd_shifted <- ff_agg_rd

for (file_nb in seq_len(length(table(exprs(ff_agg_rd$File))))){
  print(paste0("Reading in file ", file_nb))
  for (marker in colnames(ff_agg_rd@exprs)[c(3,17,28:62,71)]){
    selection_file <- which(ff_agg_rd@exprs[,"File"]==file_nb)
    min_file <- quantile(ff_agg_rd@exprs[selection_file, marker], 0)
    if(min_file != 0){
      ff_agg_rd_shifted@exprs[selection_file, marker] <-
        ff_agg_rd@exprs[selection_file, marker] - min_file
    }
  }
}

plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_nat_date_shifted.png",
                       x_names = samp_rd$Id.Cryostem.R ,ff_agg = ff_agg_rd_shifted )



  for (file_nb in files2rescale){
    for (marker in "Ce142Di"){ #colnames(ff_agg_rd@exprs)[c(3,17,28:62,71)]){
      selection_file <- which(ff_agg_rd@exprs[,"File"]==file_nb)
      ff_agg_rd_rescaled@exprs[selection_file, marker] <-
        (((ff_agg_rd@exprs[selection_file, marker] - min_rescale[marker])/
            (range_rescale[marker])) * range_ref[marker] ) + min_ref[marker]
    }
  }

  layout(matrix(1:2, nrow = 2))
  subsample <- sample(1:sum(exprs(ff_agg_rd)[,"File"]  %in% files2rescale), 10000)
  plot(exprs(ff_agg_rd)[which(exprs(ff_agg_rd)[,"File"]  %in% files2rescale), c("File_scattered", get_channels(ff_agg_rd, "CD11a"))],
       pch = ".")
  points(files2rescale, max_2rescale[, "Ce142Di"], col = "red", pch = 19)
  points(files2rescale, min_2rescale[, "Ce142Di"], col = "red", pch = 19)
  abline(h = min_rescale["Ce142Di"], col = "red")
  abline(h = max_rescale["Ce142Di"], col = "red")

  plot(exprs(ff_agg_rd_rescaled)[which(exprs(ff_agg_rd_rescaled)[,"File"]  %in% files2rescale), c("File_scattered", get_channels(ff_agg_rd, "CD11a"))],
       pch = ".")
  points(files2rescale, max_2rescale[, "Ce142Di"], col = "red", pch = 19)
  points(files2rescale, min_2rescale[, "Ce142Di"], col = "red", pch = 19)
  abline(h = min_rescale["Ce142Di"], col = "red")
  abline(h = max_rescale["Ce142Di"], col = "red")
}


plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_nat_date_rescaled.png",
                       x_names = samp_rd$Id.Cryostem.R ,ff_agg = ff_agg_rd_rescaled )

save(ff_agg_rd_rescaled,
     file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_date_rescaled.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_date_rescaled.RData")




############################################
##########   RESCALE  PER  DAY   ###########
############################################

ref_files <- which(samp_rd$DATEOFCYTOFEXPERIMENT== names(table(samp_rd$DATEOFCYTOFEXPERIMENT))[1])

min_refs <- plyr::ldply(ref_files, function(ref_id){
  min_ref <- apply(ff_agg_rd@exprs[which(exprs(ff_agg_rd)[,"File"] == ref_id),c(3,17,28:62,71)],2,
                   function(x) quantile(x, 0.01))
})
min_ref <- apply(min_refs, 2, median)

max_refs <- plyr::ldply(ref_files, function(ref_id){
  max_ref <- apply(ff_agg_rd@exprs[which(exprs(ff_agg_rd)[,"File"] == ref_id),c(3,17,28:62,71)],2,
                   function(x) quantile(x, 0.99))
})
max_ref <- apply(max_refs, 2, median)

ff_agg_rd_rescaled <- ff_agg_rd
ff_agg_rd_rescaled_orig <- ff_agg_rd_rescaled

for (day in names(table(samp_rd$DATEOFCYTOFEXPERIMENT))[-1]){
  print(paste0("Rescaling day ", day))
  files2rescale <- which(samp_rd$DATEOFCYTOFEXPERIMENT == day)
  min_2rescale <- plyr::ldply(files2rescale, function(file_id){
    min_ref <- apply(ff_agg_rd@exprs[which(exprs(ff_agg_rd)[,"File"] == file_id),c(3,17,28:62,71)],2,
                     function(x) quantile(x, 0.01))
  })
  min_rescale <- apply(min_2rescale, 2, median)

  max_2rescale<- plyr::ldply(files2rescale, function(file_id){
    max_ref <- apply(ff_agg_rd@exprs[which(exprs(ff_agg_rd)[,"File"] == file_id),c(3,17,28:62,71)],2,
                     function(x) quantile(x, 0.99))
  })
  max_rescale <- apply(max_2rescale, 2, median)

  range_rescale <- max_rescale - min_rescale
  range_ref <- max_ref - min_ref



  for (file_nb in files2rescale){
    for (marker in "Ce142Di"){ #colnames(ff_agg_rd@exprs)[c(3,17,28:62,71)]){
      selection_file <- which(ff_agg_rd@exprs[,"File"]==file_nb)
      ff_agg_rd_rescaled@exprs[selection_file, marker] <-
        (((ff_agg_rd@exprs[selection_file, marker] - min_rescale[marker])/
           (range_rescale[marker])) * range_ref[marker] ) + min_ref[marker]
    }
  }

  layout(matrix(1:2, nrow = 2))
  subsample <- sample(1:sum(exprs(ff_agg_rd)[,"File"]  %in% files2rescale), 10000)
  plot(exprs(ff_agg_rd)[which(exprs(ff_agg_rd)[,"File"]  %in% files2rescale), c("File_scattered", get_channels(ff_agg_rd, "CD11a"))],
       pch = ".")
  points(files2rescale, max_2rescale[, "Ce142Di"], col = "red", pch = 19)
  points(files2rescale, min_2rescale[, "Ce142Di"], col = "red", pch = 19)
  abline(h = min_rescale["Ce142Di"], col = "red")
  abline(h = max_rescale["Ce142Di"], col = "red")

  plot(exprs(ff_agg_rd_rescaled)[which(exprs(ff_agg_rd_rescaled)[,"File"]  %in% files2rescale), c("File_scattered", get_channels(ff_agg_rd, "CD11a"))],
       pch = ".")
  points(files2rescale, max_2rescale[, "Ce142Di"], col = "red", pch = 19)
  points(files2rescale, min_2rescale[, "Ce142Di"], col = "red", pch = 19)
  abline(h = min_rescale["Ce142Di"], col = "red")
  abline(h = max_rescale["Ce142Di"], col = "red")
}


plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_nat_date_rescaled.png",
                       x_names = samp_rd$Id.Cryostem.R ,ff_agg = ff_agg_rd_rescaled )

save(ff_agg_rd_rescaled,
     file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_date_rescaled.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_date_rescaled.RData")



###################
##### FlowSOM #####
###################

### Extract markersToPlot
markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
pheno_marks <- pheno_marks[-which(pheno_marks %in% c("CD45", "CD19"))] # all cells are already pregated on CD45+
markersToPlot <-sapply(pheno_marks, function(m){
  names(prettyMarkerNames[which(prettyMarkerNames==m)])})
names(pheno_marks) <- markersToPlot
colsToUse<-markersToPlot

#### FlowSOM ####
seed <- 1
fsom_rd <- FlowSOM(ff_agg_rd,
                    colsToUse = colsToUse,
                    scale = FALSE,
                    xdim = 15, ydim = 15, # larger grid because ++ markers
                    nClus = 50,
                    seed = seed)
PlotStars(UpdateNodeSize(fsom_rd$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD20",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))],
          view = "MST")

hospital <- samp_rd %>%
  select(HOSPITAL) %>%
  mutate (File = 1:nrow(samp_rd))

ff_agg_info <- as.data.frame(exprs(ff_agg_rd)) %>%
  left_join(hospital, by = "File")

PlotPies(UpdateNodeSize(fsom_rd$FlowSOM, maxNodeSize = 8, reset = TRUE),
         cellTypes = as.factor(ff_agg_info$HOSPITAL))

save(fsom_rd,
     file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/fsom_rd.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/fsom_rd.RData")


## generate percentages :

pctgs_rd <- generate_pctgs(
  recip_names = rd_names,
  fsom = fsom_rd,
  pdf_name = "Plot_Stars_140rd.pdf",
  fcs_dir = fcs_dir,
  output_dir = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/"
)

save(pctgs_rd, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/pctgs_rd.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/pctgs_rd.RData")
tsne <- Rtsne::Rtsne(pctgs_rd)
plot(tsne$Y, col = as.factor(samp_rd$HOSPITAL))

pctgs_meta_rd <- t(apply(pctgs_rd, 1, function(x){tapply(x, fsom_rd$metaclustering, sum)}))

hc <- hclust(dist(pctgs_rd))
library(dendextend)
dend <- as.dendrogram(hc)
labels_colors(dend) <- rainbow(13)[sort_levels_values(
  as.numeric(as.factor(samp_rd$DATEOFCYTOFEXPERIMENT))[order.dendrogram(dend)]
)]
dend <- set(dend, "labels_cex", 0.7)

short_names <- names(rd_names)
long_names <- as.character(rd_names)
names(short_names) <- long_names
reord_names <- short_names[names(labels_colors(dend))]

plot(dend,
     main = "Clustering of the datasets",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = names(table(reord_names)), fill = names(table(labels_colors(dend))),cex=0.75)

