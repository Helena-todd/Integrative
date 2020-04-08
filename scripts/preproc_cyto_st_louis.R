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


#####################################################
#####  identify matching recipients and donors  #####
#####################################################

fcs_dir <- "~/Documents/Projects/Integrative_Paris/fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)

# remove samples that were generated on 2017/10/25 and 2018/01/05 because all
# cells are CD19+ in these samples
rd_names <- fcs_names[-which(names(fcs_names)%in%c("12R","18R","12D","D1071","D369"))]

# select only the files for which we have the donor and the recipient:
rd_names <- c("D1073","D1502","D2031","D2497","R1044","R1503","R2032","R2498",
              "D1302","D2013","D2901","D391","D835","R2014","R2798","R346",
              "R836","R997","D1949","D1964","D217","D2670","D2944","D562",
              "D708","R1950","R1965","R212","R2671","R2945","R563","R709",
              "D1008","D218","D2793","D597","D689","D829","R1009","R219",
              "R2794","R598","R690","R830","D1169","D1641","D2361","D2850",
              "D296","D766","R1152","R1329","R2322","R2851","R297","R610",
              "R767","03D","03R","07D","07R","08D","08R","09D","09R","R2257",
              "D2256","D609")


####################################################################
############  aggregate and arcsinh transform samples  #############
####################################################################

set.seed(1)
ff_agg <- FlowSOM::AggregateFlowFrames(paste0(fcs_dir, "/", rd_names),
                                       cTotal= 10000*length(rd_names),
                                       writeOutput = TRUE,
                                       outputFile = paste0(fcs_dir, "/",
                                                           "aggregate_rd.fcs"))
ff_agg_rd <- flowCore::transform(ff_agg,
             flowCore::transformList(colnames(ff_agg@exprs)[c(3,17,28:62,71)],
                                     flowCore::arcsinhTransform(b=1/5, a=0, c=0)))


######################################################
########   Rescale values of 2 first donors   ########
######################################################

# 2 first files were very stretched compared to all others
# rescale them using the quantiles of a file from the same date as reference
files2rescale <- which(names(rd_names) %in% c("D1073", "D1502"))
ref_file <- which(names(rd_names) %in% c("D2031"))

min_ref <- apply(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==ref_file),
                                 c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.001))
max_ref <- apply(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==ref_file),
                                 c(3,17,28:62,71)],2,
                 function(x) quantile(x, 0.999))

for (file_nb in files2rescale){
  for (marker in colnames(ff_agg_rd@exprs)[c(3,17,28:62,71)]){
    ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==file_nb), marker] <-
      scales::rescale(ff_agg_rd@exprs[which(ff_agg_rd@exprs[,"File"]==file_nb), marker],
                      to = c(min_ref[marker], max_ref[marker]))
  }
}

