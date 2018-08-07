
#' generate_counts
#'
#' Generates a matrix of counts and pdf fsom files of all cells mapped to fsom for each patient
#'
#' @param recip_names vector containing the names of recipients' fcs files
#' @param fsom FlowSOM object
#' @param pdf_name choose a name for the exported pdf which will contain one fsom map per patient
#' @param fcs_dir path to the directory containing fcs files, corresponding to recip_names
#'
#' @return a matrix and pdf
#' @export
#'
#' @examples bla
generate_counts <- function(recip_names, fsom, pdf_name, fcs_dir){
  counts <- matrix(0,
                   length(recip_names),
                   ncol = fsom$FlowSOM$map$nNodes,
                   dimnames = list(recip_names,
                                   as.character(1:fsom$FlowSOM$map$nNodes)))
  i <- 1
  pdf(file = pdf_name)
  for (i in seq_along(recip_names)){
    file<-recip_names[[i]]
    message(file)
    ff <- read.FCS(file.path(fcs_dir,file))
    ff <- transform(ff,
                    transformList(colnames(fsom$FlowSOM$data)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
    fsom_tmp <- NewData(fsom$FlowSOM, ff)
    name<-names(recip_names)[i]
    PlotStars(fsom_tmp,main = name)

    #PlotPies(UpdateNodeSize(fsom_tmp, reset= TRUE, maxNodeSize = 8),
    #          backgroundValues = fsom$metaclustering,
    #          main = file)

    t <- table(fsom_tmp$map$mapping[,1])
    counts[file, names(t)] <- t
  }
  dev.off()
  return(counts)
}
