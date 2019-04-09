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
  library("scales")
  library("dendextend")
  library("cowplot")
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







##############################################
####    backbone on donors + all recip    ####
##############################################

fcs_dir <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
fcs_names[["R1044"]] <- "20170530_R1044_01_normalized_livecellswithoutbeads__livecells__ADN__time__Ungated_____ICOScleaned__Ungated_"
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

#####################################################
######## aggregate flowframes of all R and D ########
#####################################################

# ff_agg_rd <- fcs_to_agg(fcs_dir= fcs_dir,
#                          fcs_names= rd_names,
#                          seed = 1,
#                          cTotal = 10000*length(rd_names),
#                          output_name = "aggregate_rd.fcs")
#
# prettyMarkerNames <- ff_agg_rd@parameters@data[,"desc"] #change names of markers in flowSOM
# prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
# prettyMarkerNames[is.na(prettyMarkerNames)] <-
#   ff_agg_rd@parameters@data[,"name"][is.na(prettyMarkerNames)]
# names(prettyMarkerNames) <- colnames(ff_agg_rd)

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/ff_agg_rd.RData")

# Values of 2 first donors were too high compared to the others -> rescale
files2rescale <- which(names(rd_names) %in% c("D1073", "D1502"))
ref_file <- which(names(rd_names) %in% c("D2031"))

min_ref <- apply(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==ref_file),c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.999))

# for (file_nb in files2rescale){
#   for (marker in colnames(ff_agg_rd@exprs)[c(3,17,28:62,71)]){
#     ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==file_nb), marker] <-
#       scales::rescale(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==file_nb), marker],
#                       to = c(min_ref[marker], max_ref[marker]))
#   }
# }
# # and plot results after scaling:
# plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "DATEOFCYTOFEXPERIMENT",
#                        prettyMarkerNames, pheno_marks, png_name= "Aggregate_rescaled_date_68rd.png",
#                        x_names = rd_names, ff_agg = ff_agg_rd )
#
# save(ff_agg_rd, file = "ff_agg_rd.RData")

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
# seed <- 1
# fsom_rd <- FlowSOM(ff_agg_rd,
#                     colsToUse = colsToUse,
#                     scale = FALSE,
#                     xdim = 15, ydim = 15, # larger grid because ++ markers
#                     nClus = 30,
#                     seed = seed)
# PlotStars(UpdateNodeSize(fsom_rd$FlowSOM, maxNodeSize = 8, reset = TRUE),
#           markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD11a","CD16","CD127","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
#                                                                            "CD38","CD161","CCR7","CD27","CCR4","CCR5","CD5","CXCR3","Fas",
#                                                                            "foxP3","CD24","CXCR5"))],
#           view = "MST")
#
# PlotNumbers(UpdateNodeSize(fsom_rd$FlowSOM, maxNodeSize = 8, reset = TRUE),
#             fontSize = .5, view = "grid")

#save(fsom_rd, file = "fsom_rd.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/fsom_rd.RData")

# identify cell type for each fsom cluster
# identify_fsom_cellPopulations(fsom = fsom_rd,
#                               prettyMarkerNames, cellTypes,
#                               pdf_name = "identify_clusters_rd_tree_view.pdf",
#                               view="MST")

# pctgs_rd <- generate_pctgs(
#   recip_names = rd_names,
#   fsom = fsom_rd,
#   pdf_name = "Plot_Stars_68rd.pdf",
#   fcs_dir = fcs_dir,
#   output_dir = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/",
#   min_ref = min_ref,
#   max_ref = max_ref,
#   files2rescale = c("D1073", "D1502")
# )
#pctgs <- t(apply(counts, 1, function(x){x/sum(x)}))
rownames(pctgs_rd) <- names(rownames(pctgs_rd))
save(pctgs_rd, file = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_rd.RData")
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_rd.RData")

# manual annotation with Laetitia's labels:
#bb_rd <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_14:01:19/Pop ID Backbone R&Dall.xlsx")
#bb_rd_ordered <- bb_rd[order(bb_rd$Cluster),]
bb_rd <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/CYTOF/PopID_Backbone_R&Dall_Final_1.04.19.xlsx")
bb_rd_ordered <- bb_rd[order(bb_rd$Cluster.FS),]
fsom_rd$metaclustering <- bb_rd_ordered[,3]

pctgs_meta_rd <- t(apply(pctgs_rd, 1, function(x){tapply(x, fsom_rd$metaclustering, sum)}))
save(pctgs_meta_rd, file="/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_1_04.RData")
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_1_04_with_metaclust_labels.RData")

# with Laetitia's names :
meta_rd_labels <- seq_along(table(bb_rd_ordered[,3]))
new_labels <- lapply(seq_along(meta_rd_labels), function(i){
  meta_rd_labels[i] <- bb_rd[which(bb_rd[,3] == i),2][1]
})
colnames(pctgs_meta_rd) <- new_labels

#colnames(pctgs_meta_rd) <- names(table(bb_rd_ordered$Population.Correction.with.FlowJo))
save(pctgs_meta_rd, file="/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_1_04_with_metaclust_labels.RData")
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_18_03_with_metaclust_labels.RData")
write.xlsx(pctgs_meta_rd, file = "~/Desktop/table_pctgs_old_metaclusters_DR.xlsx", row.names=T)

load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_with_metaclust_labels.RData")
head(pctgs_meta_rd[,1:5])
write.xlsx(pctgs_meta_rd, file = "~/Desktop/table_pctgs_18_03_metaclusters.xlsx", row.names=T)

#######################################
########  FUNCTIONAL MARKERS  #########
#######################################

load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/fsom_meta_rd_old.RData")
#wsp_file <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/Threshold_41BB_ICOScleanedFCS.wsp"
wsp_file <- "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/Threshold_functionalmarkers_ICOScleaned_14032019_LagT .wsp"

wsp <- openWorkspace(wsp_file)
gates <- parseWorkspace(wsp, "All Samples", sampNloc = "sampleNode")

plot(gates)
leaf_nodes <- c("41BB+", "CD24+", "CD25+", "CD38+", "CTLA4+", "Granzyme B+",
                "HLADR+", "ICOS+", "IL10+", "Lag3+", "OX40+", "PD1+", "Tim3+")
leaf_nodes <- getNodes(gates)[-1]

gates_matrix <- lapply(gates, function(x){
  getIndiceMat(x, paste(leaf_nodes, collapse = "|"))
})

gates_manual <- lapply(gates_matrix, function(x){
  FlowSOMworkshop::manual_vector(x, leaf_nodes)
})
names(gates_manual) <- gsub("_[0-9]*$", "", names(gates_manual))
names(gates_manual)[1] <- "20170530_D2031_01_livecellswithoutbeads__livecells__ADN__time__Ungated____.fcs"

#file <- names(gates_manual)[6]
result <- list()
for(file in names(gates_manual)){
  ff <- read.FCS(file.path(fcs_dir, file))
  ff <- flowCore::transform(ff,
                            flowCore::transformList(colnames(ff)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
  fsom_tmp <- NewData(fsom_meta_rd$FlowSOM, ff)

  cluster_assignment <- table(GetClusters(fsom_tmp), gates_manual[[file]])
  cluster_labels <- colnames(cluster_assignment)[-1][apply(cluster_assignment[,-1], 1, which.max)]

  celltype_colors <- grDevices::colorRampPalette(c("white", "#00007F",
                                                   "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))(19)
  names(celltype_colors) <- levels(gates_manual[[file]])
  file_red <- gsub("^[0-9]*_([^_]*)_.*", "\\1", file)
  pdf(file = paste0("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/plots/cyto/3_backbones/backbone_2_D&Rall/fsom_1_04_funct_markers/patient_", file_red,".pdf"))
  PlotPies(fsom_tmp, gates_manual[[file]],
           main = paste0(file_red,"_", samp_rd[file_red, "GROUP"], "_", samp_rd[file_red, "DATEOFCYTOFEXPERIMENT"]))
  dev.off()
  result[[file]] <- list(cluster_labels = cluster_labels,
                         cluster_assignment = cluster_assignment)
}

save(result, file = "outputs/data/cyto/3_backbones/backbone_2_D&Rall/ratios_funct_meta_rd_1_04.RData")

###################################################
### re-arranging the functional marker matrices ###
###################################################

load("outputs/data/cyto/3_backbones/backbone_2_D&Rall/ratios_funct_meta_rd_1_04.RData")

tables_res <- lapply(2:14, function(marker){
  names(marker) <- colnames(result[[1]]$cluster_assignment)[marker]
})

nClus = 38

for(pat in seq_along(result)){
  for (marker in 2:14){
    if(length(result[[pat]]$cluster_assignment[,marker]) != nClus){
      patient_vector <- rep(0, nClus)
      names(patient_vector) <- 1:nClus
      patient_vector[names(result[[pat]]$cluster_assignment[,marker])] <- result[[pat]]$cluster_assignment[,marker]

      tables_res[[marker-1]] <- rbind(tables_res[[marker-1]], patient_vector)
    } else {
      tables_res[[marker-1]] <- rbind(tables_res[[marker-1]], result[[pat]]$cluster_assignment[,marker])
    }
  }
}

names(tables_res) <- colnames(result[[1]]$cluster_assignment)[-1]
tables_res <- lapply(tables_res, function(mat){
  mat2 <- apply(mat[-1,], 1, function(row_mat){
    row_mat <- as.numeric(row_mat)
    row_mat <- row_mat/(sum(row_mat))
  })
  rownames(mat2) <- paste0(as.character(mat[1,1]), "_", 1:nClus)
  t(mat2)
})

save(tables_res, file = "outputs/data/cyto/3_backbones/backbone_2_D&Rall/funct_mark_res_1_04.RData")

funct_big_table <- do.call(cbind, tables_res)

load("outputs/data/cyto/3_backbones/backbone_2_D&Rall/funct_mark_res_18_03.RData")



# fsom with metaclusters instead of clusters:
set.seed(1)
fsom_meta_rd <- fsom2fsom_meta(fsom = fsom_rd, colsToUse = markersToPlot,
                               pctgs_patients = pctgs_rd, plot_size = "equal_size")

PlotStars(fsom_meta_rd$FlowSOM,
          markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))])
PlotNumbers(fsom_meta_rd$FlowSOM)

save(fsom_meta_rd, file = "/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/fsom_meta_rd_1_04.RData")
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/fsom_meta_rd_old.RData")
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/fsom_meta_rd.RData")

meta_rd_labels <- seq_len(nrow(fsom_meta_rd$FlowSOM$map$codes))
new_labels <- lapply(seq_along(meta_rd_labels), function(i){
  meta_rd_labels[i] <- bb_rd[which(bb_rd$Metacluster.15.03==i),2][1]
})

PlotLabels(fsom_meta_rd$FlowSOM, labels = new_labels, fontSize = .8)

group = samp_rd[rownames(pctgs_meta_rd), "GROUP"]
PlotPies(fsom_meta_rd$FlowSOM, cellTypes = group[ff_agg_rd@exprs[,"File"]])

orig_pctgs <- table(group[ff_agg_rd@exprs[,"File"]])/sum(table(group[ff_agg_rd@exprs[,"File"]]))
df1 <- df <- data.frame(
  group = names(orig_pctgs),
  value = as.numeric(orig_pctgs)
)

plot_pie(df1, c("white", "cyan", "red"))

patient_type = rep(0, length(rownames(samp_rd)))
patient_type[grep("R", rownames(samp_rd))] <- "Recipient"
patient_type[grep("D", rownames(samp_rd))] <- "Donor"
PlotPies(fsom_meta_rd$FlowSOM, cellTypes = patient_type[ff_agg_rd@exprs[,"File"]])

apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395")
strange_normal <- rep("Normal", length(rownames(samp_rd)))
strange_normal[which(rownames(samp_rd) %in% apart)] <- "Strange"
PlotPies(fsom_meta_rd$FlowSOM, cellTypes = strange_normal[ff_agg_rd@exprs[,"File"]])

PlotPies(fsom_meta_rd$FlowSOM, cellTypes = group[ff_agg_rd@exprs[,"File"]])

orig_pctgs_2 <- table(strange_normal[ff_agg_rd@exprs[,"File"]])/sum(table(strange_normal[ff_agg_rd@exprs[,"File"]]))
df1 <- df <- data.frame(
  group = names(orig_pctgs_2),
  value = as.numeric(orig_pctgs_2)
)

plot_pie(df1, c("white", "red"))


### Plot pctgs' PCA:
pca <- prcomp(pctgs_meta_rd)
plot(pca$x)

distances = c()

for (i in names(table(samp_rd$COUPLENUMBER))){
  print(i)
  group_status <- samp_rd$GROUP[which(samp_rd$COUPLENUMBER==i)]
  pt_coord <- pca$x[rownames(samp_rd)[which(samp_rd$COUPLENUMBER==i)],c(1:2)]
  if(group_status[1]=="non_tolerant"){
    lines(pt_coord,col="red")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  } else if(group_status[1]=="primary_tolerant"){
    lines(pt_coord,col="green")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  } else {
    lines(pt_coord,col="blue")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  }
  distances = c(distances, sqrt((pt_coord[2,1]-pt_coord[1,1])^2 + (pt_coord[2,2]-pt_coord[1,2])^2))
}

list_couples <- names(table(samp_rd$COUPLENUMBER))
names(distances) <- paste0("couple_",samp_rd$COUPLENUMBER[which(samp_rd$COUPLENUMBER %in% list_couples)][1:34])
dis_colors <- c("red","green","blue")[samp_rd$GROUP[which(samp_rd$COUPLENUMBER %in% list_couples)][1:34]]

age <- samp_rd$DOB
age <- lapply(strsplit(as.character(samp_rd$DOB), "-"), function(day){
  day[1]
})
samp_rd <- samp_rd %>% mutate(AGE = round(as.numeric(unlist(age)), digits = -1))

dis_colors <- c("blue","green","yellow","orange","red")[as.factor(samp_rd$AGE[which(samp_rd$COUPLENUMBER %in% list_couples)][1:34])]
order_dis <- order(distances)

plot(distances[order_dis], col = dis_colors[order_dis], pch=19,
     main = "Distance between donor and recipient from same couple",
     ylab = "Distance D/R", xaxt = "n", xlab = "")
axis(1, at=1:34, labels=FALSE)
text(x=1:34, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels=names(distances[order_dis]), srt=45, adj=1, xpd=TRUE, cex = .7)
legend("topleft", c("non_tolerant","primary_tolerant","secondary_tolerant"),
       col = c("red","green","blue"), pch = 19)

legend("topleft", names(table(samp_rd$AGE)),
       col = c("blue","green","yellow","orange","red","pink"), pch=19)


##########
### plotting a heatmap to see how some metaclusters are changing in the D and R:
### (this figure is not very easy to analyse, probably useless)

pheatmap::pheatmap(pctgs_meta_rd[,"B naive cells", "CD4 Naive cells",
                                 "CD4 Treg", "CD8 naives", "CD8 Tc1 (TEMRA)"])

pheatmap::pheatmap(pctgs_meta_rd[,c(2, 6, 18, 21, 24)])

mat2plot <- pctgs_meta_rd[,c(2, 6, 18, 21, 24)]
mat_info <- merge.data.frame(mat2plot, samp_rd, by = "row.names" ) %>%
  group_by(COUPLENUMBER)

values <- mat_info %>%
  summarise_if(is.numeric, ~abs(.[1]-.[2]))

couple_gr <- rep(0, nrow(values))
non_tol <- samp_rd$COUPLENUMBER[which(samp_rd$GROUP=="non_tolerant")]
couple_gr[which(values$COUPLENUMBER %in% non_tol)] <- "non_tol"
tol_1 <- samp_rd$COUPLENUMBER[which(samp_rd$GROUP=="primary_tolerant")]
couple_gr[which(values$COUPLENUMBER %in% tol_1)] <- "primary_tol"
tol_2 <- samp_rd$COUPLENUMBER[which(samp_rd$GROUP=="secondary_tolerant")]
couple_gr[which(values$COUPLENUMBER %in% tol_2)] <- "secondary_tol"

val_ordered <- values %>% mutate(couple_gr = couple_gr) %>%
  arrange(couple_gr) %>%
  column_to_rownames("COUPLENUMBER")



pheatmap::pheatmap(val_ordered[,1:5], cluster_rows = F, annotation_row = val_ordered[,c(6,6)])


## correlations between donors and recips :
pctgs_cor <- merge(data.frame(pctgs_meta_rd), samp_rd, by="row.names")
pctgs_corD <- pctgs_cor[grep(pattern = "D", pctgs_cor$Row.names),] %>%
  arrange(COUPLENUMBER)
pctgs_corR <- pctgs_cor[grep(pattern = "R", pctgs_cor$Row.names),]%>%
  arrange(COUPLENUMBER)

correlations <- sapply(seq_len(nrow(pctgs_corD)), function(i){
  cor(as.numeric(pctgs_corD[i, 2:37]),
      as.numeric(pctgs_corR[i, 2:37]))
})

correlations <- 1-correlations # to compare it to the distances

names(correlations) <- paste0("couple_",pctgs_corD$COUPLENUMBER)
cor_colors <- c("red","green","blue")[pctgs_corD$GROUP]
order_cor <- order(correlations)
plot(correlations[order_cor], col = cor_colors[order_cor], pch=19,
     main = "Correlation between donor and recipient from same couple",
     ylab = "1 - Correlation D/R", xaxt = "n", xlab ="")
axis(1, at=1:34, labels=FALSE)
text(x=1:34, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels=names(correlations[order_cor]), srt=45, adj=1, xpd=TRUE, cex = .7)
legend("topleft", c("non_tolerant","primary_tolerant","secondary_tolerant"),
       col = c("red","green","blue"), pch = 19)




#Distribution of patients in the 36 metaclusters?
png("~/Documents/VIB/Projects/Integrative_Paris/Integrative/distrib_patients.png")
barplot(pctgs_cor[,2], col = c("red","blue","green")[pctgs_cor$GROUP])

png("~/Documents/VIB/Projects/Integrative_Paris/Integrative/distrib_patients.png",
    width = 9000,
    height = 4500)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:40, nrow=8, byrow = TRUE))
for(metacluster in 2:37){
  print(paste0("Plotting metacluster ",metacluster-1))
  barplot(pctgs_cor[,metacluster], col = c("red","green","blue")[pctgs_cor$GROUP],
          main = new_labels[metacluster-1])
}
plot.new()
legend("center",legend=levels(pctgs_cor$GROUP), col=c("red","blue","green"),
       pch=19, cex=3)
dev.off()





### paired t.tests between D and R per group
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_with_metaclust_labels.RData")
pctg_data <- merge(pctgs_meta_rd, samp_rd, by="row.names")
rown <- pctg_data$Row.names
pctg_data <- pctg_data[,c("COUPLENUMBER", "GROUP", colnames(pctgs_meta_rd))]#[,c(41, 39, 2:37)]
colnames(pctg_data)[c(1,2)] <- c("couplenb", "group")
rownames(pctg_data) <- rown

pctg_D <- pctg_data[grep(pattern = "D", rownames(pctg_data)),]
pctg_D <- pctg_D[order(pctg_D$couplenb),]

pctg_R <- pctg_data[grep(pattern = "R", rownames(pctg_data)),]
pctg_R <- pctg_R[order(pctg_R$couplenb),]

ttests <- function(mat_D, mat_R, group){
  mat_D_gr <- mat_D[which(mat_D$group==group),]
  mat_R_gr <- mat_R[which(mat_R$group==group),]

  ttests_gr <- sapply(3:ncol(mat_D), function(pctgp){
    paired_ttest <- t.test(as.numeric(mat_D_gr[, pctgp]),
                           as.numeric(mat_R_gr[, pctgp]), paired = T)
    pval <- paired_ttest$p.value
    names(pval) <- colnames(pctg_D)[pctgp]
    pval
  })
  return(ttests_gr)
}

ttests_NT <- ttests(pctg_D, pctg_R, "non_tolerant")
# correct for multiple tests
ttests_NT <- p.adjust(ttests_NT, "BH")

DE_NT <- ttests_NT[which(ttests_NT <0.05)]
table_NT <- data.frame(as.numeric(sort(DE_NT)))
rownames(table_NT) <- names(sort(DE_NT))
colnames(table_NT) <- "pvalue"
table_NT
length(DE_NT)
# 18

ttests_1T <- ttests(pctg_D, pctg_R, "primary_tolerant")
# correct for multiple tests
ttests_1T <- p.adjust(ttests_1T, "BH")

DE_1T <- ttests_1T[which(ttests_1T <0.05)]
table_1T <- data.frame(as.numeric(sort(DE_1T)))
rownames(table_1T) <- names(sort(DE_1T))
colnames(table_1T) <- "pvalue"
table_1T
length(DE_1T)
# 0

ttests_2T <- ttests(pctg_D, pctg_R, "secondary_tolerant")
# correct for multiple tests
ttests_2T <- p.adjust(ttests_2T, "BH")

DE_2T <- ttests_2T[which(ttests_2T <0.05)]
table_2T <- data.frame(as.numeric(sort(DE_2T)))
rownames(table_2T) <- names(sort(DE_2T))
colnames(table_2T) <- "pvalue"
table_2T
length(DE_2T)
# 0


### clustering on the R/D couples
big_mat <- merge(pctgs_meta_rd, samp_rd, by = "row.names")
big_D <- big_mat[grep("D", big_mat$Row.names),] %>%
  arrange(COUPLENUMBER)
big_R <- big_mat[grep("R", big_mat$Row.names),]%>%
  arrange(COUPLENUMBER)

pctg_couples <- matrix(rep(0, 108*34), nrow = 34)
for (i in 1:nrow(big_D)){
  d_vec <- big_D[i, 2:37]
  r_vec <- big_R[i, 2:37]
  both <- rbind(d_vec, r_vec)
  means_both <- apply(both,2,mean)
  pctg_couples[i,] <- as.numeric(c((r_vec-d_vec), means_both, abs(r_vec-d_vec)))
}
rownames(pctg_couples) <- paste0("couple_",big_D$COUPLENUMBER)

hc <- hclust(dist(pctg_couples))
plot(hc, col = big_D$GROUP)


dend <- as.dendrogram(hc)
grop<-as.numeric(big_D$GROUP)
labels_colors(dend) <-
  c("red","blue","green")[sort_levels_values(
    as.numeric(grop)[order.dendrogram(dend)]
  )]
dend <- set(dend, "labels_cex", 0.7)

plot(dend,
     main = "Clustering on metaclusters only",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = c("non_tolerant","primary_tol","secondary_tol"), fill = c("red","green","blue"),cex=0.75)


apart<- c("R690","R830","R219","R598","R2798","R836","R2589","03R","R419","R395")
which(big_R$Row.names%in% apart)


grop<-rep(1, 34)
grop[which(big_R$Row.names%in% apart)] <- 2
labels_colors(dend) <-
  c("black", "red")[sort_levels_values(
    as.numeric(grop)[order.dendrogram(dend)]
  )]
dend <- set(dend, "labels_cex", 0.7)

plot(dend,
     main = "Clustering on metaclusters only",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = c("normal","CMV_strange"), fill = c("black","red"),cex=0.75)


### RF on the DR couples:
table1 <- data.frame(cbind(pctg_D$group ,pctg_couples))
table1[,1] <- as.factor(table1[,1])
colnames(table1) <- c("group", paste0("substraction_",1:36), paste0("mean_",1:36),
                      paste0("abs_",1:36))
rf_r<-randomForest(group~., table1, ntree=5000, mtry=30)
rf_r
plot(rf_r, main = "Random Forest on the D/R couple CYTOF info")









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



