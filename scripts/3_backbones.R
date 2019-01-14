### Idea: generate 3 different flowSOM backbones and compare them
# 1st backbone: generated only on donors
# 2nd backbone: generated on donors + normal recips
# 3rd backbone: on all complete cases (where we have donors and recip couples)

suppressPackageStartupMessages({
  library("openxlsx")
  library("FlowSOM")
  library("tidyverse")
  library("magrittr")
  library("flowCore")
  library("flowWorkspace")
  library("ggraph")
  library("igraph")
})
library(BioGVHD)
options("scipen"=100)

################################################
#####  identify matching recip and donors  #####
################################################

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
rd_names <- fcs_names[-which(names(fcs_names)%in%c("12R","18R","12D","D1071","D369"))]

## import metadat info :
dataPath <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx"
samp_rd <- import_patient_info(data_synthesis_file = dataPath,
                               patient_names = rd_names,
                               patient_type = "donor")

rd_couples <- samp_rd$COUPLENUMBER[which(duplicated(samp_rd$COUPLENUMBER)==TRUE)]
samp_rd <- samp_rd[which(samp_rd$COUPLENUMBER %in% rd_couples),]

rd_names <- rd_names[which(names(rd_names) %in% samp_rd$Id.Cryostem.R)]






###############################################
#######     backbone on donors only     #######
###############################################

ff_agg_d <- fcs_to_agg(fcs_dir= fcs_dir,
                        fcs_names= donor_names,
                        seed = 1,
                        cTotal = 10000*length(donor_names),
                        output_name = "aggregate_d.fcs")

prettyMarkerNames <- ff_agg_d@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_d@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_d)

# Values of 2 first donors were too high compared to the others -> rescale
files2rescale <- which(names(donor_names) %in% c("D1073", "D1502"))
ref_file <- which(names(donor_names) %in% c("D2031"))

min_ref <- apply(ff_agg_d@exprs[which(ff_agg_d@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_d@exprs[which(ff_agg_d@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.999))

for (file_nb in files2rescale){
  for (marker in colnames(ff_agg_d@exprs)[c(3,17,28:62,71)]){
    ff_agg_d@exprs[which(ff_agg_d@exprs[,"File"]==file_nb), marker] <-
      scales::rescale(ff_agg_d@exprs[which(ff_agg_d@exprs[,"File"]==file_nb), marker],
                      to = c(min_ref[marker], max_ref[marker]))
  }
}
# and plot results after scaling:
plot_aggregate_markers(patient_names = donor_names, samp_patients=samp_donor, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_rescaled_date_34donors.png",
                       ff_agg = ff_agg_d )

### markersToPlot
markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
pheno_marks <- pheno_marks[-which(pheno_marks=="CD45")] # all cells are already pregated on CD45+
markersToPlot <-sapply(pheno_marks, function(m){
  names(prettyMarkerNames[which(prettyMarkerNames==m)])})
names(pheno_marks) <- markersToPlot
colsToUse<-markersToPlot

#### FlowSOM ####
seed <- 1
fsom_donors <- FlowSOM(ff_agg_d,
                colsToUse = colsToUse,
                scale = FALSE,
                xdim = 15, ydim = 15, # larger grid because ++ markers
                nClus = 30,
                seed = seed)
PlotStars(UpdateNodeSize(fsom_donors$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))])

cl_30 <- fsom_donors$metaclustering

# redo metaclustering for nclus = 40 and nclus = 50
cl_40 <- as.factor(metaClustering_consensus(fsom_donors$FlowSOM$map$codes,
                                         40, seed = seed))
cl_50 <- as.factor(metaClustering_consensus(fsom_donors$FlowSOM$map$codes,
                                         50, seed = seed))
save(fsom_donors, file=("outputs/data/cyto/3_backbones/bb_donors.RData"))

##########################################################################
#### generate plots for Laetitia with identify_fsom_cell_population() ####
##########################################################################

## table of associations btw markers and populations ------
pop_mark <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_3:04:18/CYTOF_David_Michonneau_excel/Populations and markers_filtered.xlsx",
                      check.names=F) %>%
  dplyr::filter(!is.na(Markers)) %>%
  select(-c(CD45, GranzymeB)) # rm CD45: they are all positive (pre-gated by laetitia)
# rm granzymeB: functional marker, not phenotypic

# extract info of which markers are expressed in which cell types:
cellTypes <- markers_of_cellTypes(marker_table = pop_mark)

# identify cell type for each fsom cluster
fsom_d <- fsom_donors
#fsom_d$FlowSOM$prettyColnames <- gsub("<[A-Za-z0-9]*> ","", fsom_d$FlowSOM$prettyColnames)
identify_fsom_cellPopulations(fsom = fsom_d,
                              prettyMarkerNames, cellTypes,
                              pdf_name = "identify_clusters_donors_grid_view.pdf",
                              view="grid")


#sapply(names(cellTypes[[cellType]]), function(x) get_channels(ff_agg_d, x))
# instead of all of the identify_fsom_cellPopulations function?

PlotStars(UpdateNodeSize(fsom_donors$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))],
          view = "grid")
PlotLabels(UpdateNodeSize(fsom_donors$FlowSOM, maxNodeSize = 8, reset = TRUE),
           labels=cl_50, view = "grid")
PlotNumbers(UpdateNodeSize(fsom_donors$FlowSOM, maxNodeSize = 8, reset = TRUE),
            fontSize = .5)









###############################################
####   backbone on donors + normal recip   ####
###############################################

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
rd_names <- fcs_names[-which(names(fcs_names)%in%c("12R","18R","12D","D1071","D369"))]

apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395")
rdn_names <- rd_names[-(which(names(rd_names)%in%apart))]

## import metadat info :
dataPath <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx"
samp_rdn <- import_patient_info(data_synthesis_file = dataPath,
                               patient_names = rdn_names,
                               patient_type = "donor")

rdn_couples <- samp_rdn$COUPLENUMBER[which(duplicated(samp_rdn$COUPLENUMBER)==TRUE)]
samp_rdn <- samp_rdn[which(samp_rdn$COUPLENUMBER %in% rdn_couples),]
rdn_names <- rdn_names[which(names(rdn_names) %in% samp_rdn$Id.Cryostem.R)]

######################################################
######## aggregate flowframes of rdn (normal) ########
######################################################

ff_agg_rdn <- fcs_to_agg(fcs_dir= fcs_dir,
                       fcs_names= rdn_names,
                       seed = 1,
                       cTotal = 10000*length(rdn_names),
                       output_name = "aggregate_d.fcs")

prettyMarkerNames <- ff_agg_rdn@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_rdn@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_rdn)

# Values of 2 first donors were too high compared to the others -> rescale
files2rescale <- which(names(rdn_names) %in% c("D1073", "D1502"))
ref_file <- which(names(rdn_names) %in% c("D2031"))

min_ref <- apply(ff_agg_rdn@exprs[which(ff_agg_rdn@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_rdn@exprs[which(ff_agg_rdn@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.999))

for (file_nb in files2rescale){
  for (marker in colnames(ff_agg_rdn@exprs)[c(3,17,28:62,71)]){
    ff_agg_rdn@exprs[which(ff_agg_rdn@exprs[,"File"]==file_nb), marker] <-
      scales::rescale(ff_agg_rdn@exprs[which(ff_agg_rdn@exprs[,"File"]==file_nb), marker],
                      to = c(min_ref[marker], max_ref[marker]))
  }
}
# and plot results after scaling:
plot_aggregate_markers(patient_names = rdn_names, samp_patients=samp_rdn, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_rescaled_date_54rdn.png",
                       ff_agg = ff_agg_rdn )

save(ff_agg_rdn, file = "ff_agg_rdn.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/d&rnorm/ff_agg_rdn.RData")

### markersToPlot
markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
pheno_marks <- pheno_marks[-which(pheno_marks=="CD45")] # all cells are already pregated on CD45+
markersToPlot <-sapply(pheno_marks, function(m){
  names(prettyMarkerNames[which(prettyMarkerNames==m)])})
names(pheno_marks) <- markersToPlot
colsToUse<-markersToPlot

#### FlowSOM ####
seed <- 1
fsom_rdn <- FlowSOM(ff_agg_rdn,
                       colsToUse = colsToUse,
                       scale = FALSE,
                       xdim = 15, ydim = 15, # larger grid because ++ markers
                       nClus = 30,
                       seed = seed)
PlotStars(UpdateNodeSize(fsom_rdn$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))],
          view = "grid")

PlotNumbers(UpdateNodeSize(fsom_rdn$FlowSOM, maxNodeSize = 8, reset = TRUE),
            fontSize = .5, view = "grid")

save(fsom_rdn, file = "fsom_rdn.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/d&rnorm/fsom_rdn.RData")

# identify cell type for each fsom cluster
identify_fsom_cellPopulations(fsom = fsom_rdn,
                              prettyMarkerNames, cellTypes,
                              pdf_name = "identify_clusters_rdn_tree_view.pdf",
                              view="MST")



###############################################
####   backbone on D&R tolerant 1&2 only   ####
###############################################

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
rd_names <- fcs_names[-which(names(fcs_names)%in%c("12R","18R","12D","D1071","D369"))]

## import metadat info :
dataPath <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx"
samp_rd <- import_patient_info(data_synthesis_file = dataPath,
                                patient_names = rd_names,
                                patient_type = "donor")

samp_rdtol <- samp_rd[-which(samp_rd$GROUP=="non_tolerant"),]
rdtol_couples <- samp_rdtol$COUPLENUMBER[which(duplicated(samp_rdtol$COUPLENUMBER)==TRUE)]
samp_rdtol <- samp_rdtol[which(samp_rdtol$COUPLENUMBER %in% rdtol_couples),]
rdtol_names <- rd_names[which(names(rd_names) %in% samp_rdtol$Id.Cryostem.R)]

######################################################
######## aggregate flowframes of rd tolerants ########
######################################################

ff_agg_rdtol <- fcs_to_agg(fcs_dir= fcs_dir,
                         fcs_names= rdtol_names,
                         seed = 1,
                         cTotal = 10000*length(rdtol_names),
                         output_name = "aggregate_rdtol.fcs")

prettyMarkerNames <- ff_agg_rdtol@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_rdtol@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_rdtol)

# Values of 2 first donors were too high compared to the others -> rescale
files2rescale <- which(names(rdtol_names) %in% c("D1073", "D1502"))
ref_file <- which(names(rdtol_names) %in% c("D2031"))

min_ref <- apply(ff_agg_rdtol@exprs[which(ff_agg_rdtol@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_rdtol@exprs[which(ff_agg_rdtol@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.999))

for (file_nb in files2rescale){
  for (marker in colnames(ff_agg_rdtol@exprs)[c(3,17,28:62,71)]){
    ff_agg_rdtol@exprs[which(ff_agg_rdtol@exprs[,"File"]==file_nb), marker] <-
      scales::rescale(ff_agg_rdtol@exprs[which(ff_agg_rdtol@exprs[,"File"]==file_nb), marker],
                      to = c(min_ref[marker], max_ref[marker]))
  }
}
# and plot results after scaling:
plot_aggregate_markers(patient_names = rdtol_names, samp_patients=samp_rdtol, color_by = "DATEOFCYTOFEXPERIMENT",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_rescaled_date_26rdtol.png",
                       ff_agg = ff_agg_rdtol )

save(ff_agg_rdtol, file = "ff_agg_rdtol.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_4_D&Rtol/ff_agg_rdtol.RData")

### markersToPlot
markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]
pheno_marks <- pheno_marks[-which(pheno_marks=="CD45")] # all cells are already pregated on CD45+
markersToPlot <-sapply(pheno_marks, function(m){
  names(prettyMarkerNames[which(prettyMarkerNames==m)])})
names(pheno_marks) <- markersToPlot
colsToUse<-markersToPlot

#### FlowSOM ####
seed <- 1
fsom_rdtol <- FlowSOM(ff_agg_rdtol,
                    colsToUse = colsToUse,
                    scale = FALSE,
                    xdim = 15, ydim = 15, # larger grid because ++ markers
                    nClus = 30,
                    seed = seed)
PlotStars(UpdateNodeSize(fsom_rdtol$FlowSOM, maxNodeSize = 8, reset = TRUE),
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                                                                           "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
                                                                           "foxP3","CD24","CXCR5"))],
          view = "grid")

PlotNumbers(UpdateNodeSize(fsom_rdtol$FlowSOM, maxNodeSize = 8, reset = TRUE),
            fontSize = .5, view = "MST")

save(fsom_rdtol, file = "fsom_rdtol.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_4_D&Rtol/fsom_rdtol.RData")

# identify cell type for each fsom cluster
identify_fsom_cellPopulations(fsom = fsom_rdtol,
                              prettyMarkerNames, cellTypes,
                              pdf_name = "identify_clusters_rdtol_grid_view.pdf",
                              view="grid")



