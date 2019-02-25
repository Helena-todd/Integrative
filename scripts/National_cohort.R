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

samp_rd_all <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/National_cohort/CYTOF/Data synthesis national cohort Cryostem 22222019.xlsx")
samp_rd_all$DATEOFCYTOFEXPERIMENT <- as.Date(samp_rd_all$DATEOFCYTOFEXPERIMENT, format = "%d.%m.%Y")

## Filter out data :
samp_rd <- samp_rd_all %>%
  dplyr::filter(DELAY_SAMPLE >= 148) %>%
  dplyr::filter(!Couple.number %in% c(32, 71)) %>%
  arrange(HOSPITAL,DATEOFCYTOFEXPERIMENT)

rd_names <- fcs_names[samp_rd$Id.Cryostem.R] # rearrange per hospital for the plots


## plot aggregate :

ff_agg_rd <- fcs_to_agg(fcs_dir= fcs_dir,
                       fcs_names= rd_names,
                       seed = 1,
                       cTotal = 10000*length(rd_names),
                       output_name = "aggregate_rd_national.fcs")

prettyMarkerNames <- ff_agg_rd@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_rd@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_rd)

x_names <- paste0(samp_rd$Id.Cryostem.R, "_", samp_rd$DATEOFCYTOFEXPERIMENT)
plot_aggregate_markers(patient_names = rd_names, samp_patients=samp_rd, color_by = "HOSPITAL",
                       prettyMarkerNames, pheno_marks, png_name= "Aggregate_nat_144rd.png",
                       x_names = x_names ,ff_agg = ff_agg_rd )

save(ff_agg_rd,
     file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto_national/ff_agg_rd_nat_144.RData")


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
tsne <- Rtsne::Rtsne(pctgs_rd)
plot(tsne$Y, col = as.factor(samp_rd$HOSPITAL))

hc <- hclust(dist(pctgs_rd))
library(dendextend)
dend <- as.dendrogram(hc)
labels_colors(dend) <- rainbow(13)[sort_levels_values(
  as.numeric(as.factor(samp_rd$HOSPITAL))[order.dendrogram(dend)]
)]
dend <- set(dend, "labels_cex", 0.7)

plot(dend,
     main = "Clustering of the datasets",
     horiz =  F,  nodePar = list(cex = .007))
legend("topleft", legend = c("non_tolerant","primary_tol","secondary_tol"), fill = c("red","green","blue"),cex=0.75)

