# Libraries and functions -----------------------------------------------------

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("flowCore"))
#suppressPackageStartupMessages(library("flowAI"))
suppressPackageStartupMessages(library("flowWorkspace"))
options("scipen"=100)
suppressPackageStartupMessages(library("MetaCyto"))

# QC between samples ----------------------------------------------------------

## excel file with info about samples
samples <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018.xlsx", 
                     check.names = FALSE) %>% 
  mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
         GROUP = tolower(GROUP)) %>% 
  dplyr::filter(!is.na(FCSNAME))
rownames(samples)<- samples$Id.Cryostem.R

# samples<-samples[-which(is.na(samples$FCSNAME)),]

markers <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/PANEL CYTOF .xlsx", 
                     check.names = FALSE) ## excel file with info about markers

marker_pops<-read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_3:04:18/CYTOF_David_Michonneau_excel/Populations and markers.xlsx", 
                       check.names = FALSE) ## excel file with info about markers

seed <- 1
set.seed(seed)
fcs_dir<-"~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/fcs/"
fcs_names <- list.files(fcs_dir, 
                        pattern="^2.*fcs$")#,
#recursive = TRUE)
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)
# fcs_names<-fcs_names[samples$Id.Cryostem.R]
samples <- samples[names(fcs_names), ]

## aggregate fcs files:
setwd(dir = fcs_dir)
set.seed(seed)
ff_agg <- AggregateFlowFrames(fcs_names,  
                              cTotal=10000*88, writeOutput = TRUE,
                              outputFile = "aggregate.fcs")

## transform aggregated file:
ff_agg <- transform(ff_agg,
                    transformList(colnames(ff_agg)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))

prettyMarkerNames <- ff_agg@parameters@data[,"desc"] #change names of markers in flowSOM
prettyMarkerNames <- gsub(".*_", "", prettyMarkerNames)
prettyMarkerNames[is.na(prettyMarkerNames)] <- 
  ff_agg@parameters@data[,"name"][is.na(prettyMarkerNames)]
names(prettyMarkerNames) <- colnames(ff_agg)

## find only phenotypic markers
pheno_marks<-markers[which(markers[,4]==1),1]

## plot expression of pheno markers in patients (QC):
#markersToPlot <- colnames(ff_agg)[c(1,4,7:18)]
markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
png(file.path(paste0("Aggregate_time.png")),
    width = 9000,
    height = 4500)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:30, nrow=5, byrow = TRUE))
gr_levels<-levels(as.factor(samples$DATEOFCYTOFEXPERIMENT))
for(marker in markersToPlot){
  print(paste0("Plotting ",prettyMarkerNames[marker]," for the aggregated file."))
  plot(0,#exprs(ff_agg[,c("File_scattered",marker)]),
       type="n",
       xlab = "Files",
       ylab = prettyMarkerNames[marker],
       ylim = c(0, max(10, exprs(ff_agg[,marker]))),
       xlim= c(0,length(fcs_names)+1),
       xaxt="n")
  abline(v=seq_along(fcs_names), col="lightgrey")
  axis(side=1, at=seq_along(fcs_names), labels=names(fcs_names), las=2,cex=0.7)
  points(exprs(ff_agg[,c("File_scattered",marker)]),
         pch=".",
         col = scales::hue_pal()(10)[factor(samples$DATEOFCYTOFEXPERIMENT[exprs(ff_agg)[,"File"]],levels=gr_levels)])
  #col = scales::hue_pal()(3)[factor(samples$GROUP[exprs(ff_agg)[,"File"]],levels=gr_levels)])
}
plot.new()
legend("center",legend=gr_levels, col=scales::hue_pal()(10),
       pch=19, cex=3)
dev.off()

### 2 first processed patients (D1073 and D1502) have a higher expression for all markers
## I normalise them so that they fit in their batch from date 30.05.17

# ff_goal = 7 "standard" files from day 30-05-2017
ff_goal<- exprs(ff_agg)[which(exprs(ff_agg)[,72]%in%c(3:9)),]
ff_BE<-ff_agg

## normalise with mean sd ---------------------
mean_goal <- apply(ff_goal[,markersToPlot],2,mean)
sd_goal <- apply(ff_goal[,markersToPlot],2,sd)

for (marker in markersToPlot){
  ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker] <- scale(ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker]) * sd_goal[marker] + mean_goal[marker]
}

## normalise with quantiles 0,001 and 0,999 ---------------------
min_goal <- apply(ff_goal[,markersToPlot],2,function(x) quantile(x,0.001))
max_goal <- apply(ff_goal[,markersToPlot],2,function(x) quantile(x,0.999))

for (marker in markersToPlot){
  ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker] <- (((ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker] - min(ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker])) / 
                                                                 (max(ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker]) - min(ff_BE@exprs[which(exprs(ff_agg)[,72]%in%c(1,2)),marker]))) * (max_goal[marker] - min_goal[marker]) ) + min_goal[marker]
}

png(file.path(paste0("testtt.png")),
    width = 9000,
    height = 4500)
par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
layout(matrix(1:30, nrow=5, byrow = TRUE))
gr_levels<-levels(as.factor(samples$DATEOFCYTOFEXPERIMENT))
for(marker in markersToPlot){
  print(paste0("Plotting ",prettyMarkerNames[marker]," for the aggregated file."))
  plot(0,#exprs(ff_agg[,c("File_scattered",marker)]),
       type="n",
       xlab = "Files",
       ylab = prettyMarkerNames[marker],
       ylim = c(0, max(10, exprs(ff_agg[,marker]))),
       xlim= c(0,length(fcs_names)+1),
       xaxt="n")
  abline(v=seq_along(fcs_names), col="lightgrey")
  axis(side=1, at=seq_along(fcs_names), labels=names(fcs_names), las=2,cex=0.7)
  points(ff_BE@exprs[,c("File_scattered",marker)],
         pch=".",
         col = scales::hue_pal()(10)[factor(samples$DATEOFCYTOFEXPERIMENT[exprs(ff_agg)[,"File"]],levels=gr_levels)])
  #col = scales::hue_pal()(3)[factor(samples$GROUP[exprs(ff_agg)[,"File"]],levels=gr_levels)])
}
plot.new()
legend("center",legend=gr_levels, col=scales::hue_pal()(10),
       pch=19, cex=3)
dev.off()