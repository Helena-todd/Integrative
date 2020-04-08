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

fcs_dir <- "~/Documents/Projects/Integrative_Paris/National_fcs/"
fcs_names <- list.files(fcs_dir, pattern="^2.*fcs$")
names(fcs_names) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", fcs_names)

# select only the files for which we have the donor and the recipient:
rd_names <- c("D4147","R4147","D4316","R4316","D585","D4491","R4491","R2824",
              "R585","D1840","R1840","D3504","R3504","D3574","R3574","D557",
              "R557","D4000","R4000","D5366","R5366","D4047","R4047","D3594",
              "R3594","D3975","R3975","D3715","R3715","D3920","R3920","D2854",
              "R2854","D3636","R3636","D4617","R4617","D3209","D6312","R6312",
              "D4602","R4602","D3653","R3653","D2915","R3209","D3421","R3421",
              "D908","R908","D4045","R4045","R3030","R4596","D2896","R2896",
              "D639","R639","D4112","R4112","D3979","R985","D2048","R2048",
              "D3484","R3484","D4525","R4525","D4001","R4001","D5404","R5404",
              "R3979","D1626","R1626","D4412","R4412","D1665","R1665","D2655",
              "R2655","D5340","R5340","D1256","D3245","D170","R170","D1306",
              "R1306","D3298","R3298","D216","R216","D3096","R3096","D379",
              "R379","D5555","R5555","R1256","R3245","D1539","R1539","D3741",
              "R3741","D3917","R3917","D3939","R3939","D5242","R5242","D5280",
              "R5280","D4465","D1710","R1710","R5554","D2291","R2291","D3116",
              "R3116","D3534","R3534","D3294","R3294","D3110","R3110","D985",
              "D2824","R2915","D3824","R3824","D2986","R2986","D3747","R3747",
              "D3030","D5554","R4465","D4596")

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




