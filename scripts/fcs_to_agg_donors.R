library(FlowSOM)
library(flowCore)
library(openxlsx)
library(dplyr)

## only for donors ##

samples <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018.xlsx",
                           check.names = FALSE) %>%
  mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
         GROUP = tolower(GROUP)) %>%
  dplyr::filter(!is.na(FCSNAME))

fcs_dir<-"~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir,
                        pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)

sample_recip<- samples[grep("R",samples[,1]),]
sample_recip<- sample_recip[-which(sample_recip$Id.Cryostem.R%in%c("12R","18R")),]

donor_names<-fcs_names[grep("D",names(fcs_names))]
#donor_names<- donor_names[-which(names(donor_names)%in%c("12R","18R"))]
#donor_names<- donor_names[myord]


## aggregate donor fcs files ----------------------------------

setwd(dir = fcs_dir)
seed <- 1
set.seed(seed)
ff_agg_donor <- AggregateFlowFrames(donor_names,
                                    cTotal=10000*37, writeOutput = TRUE,
                                    outputFile = "aggregate_donor.fcs")
ff_agg_donor <- transform(ff_agg_donor,
                          transformList(colnames(ff_agg_donor)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))

prettyMarkerNames <- ff_agg_donor@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <-
  ff_agg_donor@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg_donor)

sample_donor<- samples[grep("D",samples[,1]),]
#sample_donor<- sample_donor[-which(sample_donor$Id.Cryostem.R%in%c("12R","18R")),]

## plot expression of pheno markers in patients (QC):
markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF corrigé 07.2018.xlsx",
                     check.names = FALSE) ## excel file with info about markers
pheno_marks<-markers[which(markers[,4]==1),1]
samp_donor<- sample_donor %>% arrange(DATEOFCYTOFEXPERIMENT) # reorder files by date of experiment

markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
png(file.path(paste0("Aggregate_group_donor_all_marks2.png")),
    width = 9000,
    height = 4500)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:35, nrow=7, byrow = TRUE))
gr_levels<-levels(as.factor(samp_donor$DATEOFCYTOFEXPERIMENT))
for(marker in markersToPlot){
  print(paste0("Plotting ",prettyMarkerNames[marker]," for the aggregated file."))
  plot(0,#exprs(ff_agg[,c("File_scattered",marker)]),
       type="n",
       xlab = "Files",
       ylab = prettyMarkerNames[marker],
       ylim = c(0, max(10, ff_agg_donor@exprs[,marker])),
       xlim= c(0,length(donor_names)+1),
       xaxt="n")
  abline(v=seq_along(donor_names), col="lightgrey")
  axis(side=1, at=seq_along(donor_names), labels=names(donor_names), las=2,cex=0.7)
  # points(myff[,c("File_scattered",marker)],
  #        pch=".",
  #        col = scales::hue_pal()(3)[factor(samp_donor$GROUP[myff[,"File"]],levels=gr_levels)])
  #
  points(ff_agg_donor@exprs[,c("File_scattered",marker)],
         pch=".",
         col = scales::hue_pal()(length(gr_levels))[factor(samp_donor$DATEOFCYTOFEXPERIMENT[ff_agg_donor@exprs[,"File"]],levels=gr_levels)])
}
plot.new()
legend("center",legend=gr_levels, col=scales::hue_pal()(7),
       pch=19, cex=3)
dev.off()

save(ff_agg_donor, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/ff_agg_donor")
save(donor_names, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/donor_names.RData")
save(sample_donor, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/donors/sample_donor.RData")
