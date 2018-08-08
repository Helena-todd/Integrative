
#' generate_pctgs
#'
#' Generates a matrix of pctgs and pdf fsom files of pctgs of each patient's cells mapped to fsom clusters
#'
#' @param recip_names vector containing the names of recipients' fcs files
#' @param fsom FlowSOM object
#' @param pdf_name choose a name for the exported pdf which will contain one fsom map per patient
#' @param fcs_dir path to the directory containing fcs files, corresponding to recip_names
#'
#' @return a pctgs matrix and pdf containing one map per patient
#' @export
#'
#' @examples
#' pctgs <- generate_pctgs(recip_names, fsom, pdf_name = "my_pdf.pdf", fcs_dir)
generate_pctgs <- function(recip_names, fsom, pdf_name, fcs_dir){
  pctgs <- matrix(0,
                   length(recip_names),
                   ncol = fsom$FlowSOM$map$nNodes,
                   dimnames = list(recip_names,
                                   as.character(1:fsom$FlowSOM$map$nNodes)))
  #i <- 1
  pdf(file = pdf_name)
  for (i in seq_along(recip_names)){
    file<-recip_names[[i]]
    message(file)
    ff <- read.FCS(file.path(fcs_dir,file))
    ff <- flowCore::transform(ff,
                    transformList(colnames(fsom$FlowSOM$data)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
    fsom_tmp <- NewData(fsom$FlowSOM, ff)
    name<-names(recip_names)[i]
    #PlotStars(fsom_tmp,main = name)

    counts <- table(GetClusters(fsom_tmp))
    #pctgs <- rep(0, fsom_tmp$map$nNodes)
    #names(pctgs) <- as.character(seq_len(fsom_tmp$map$nNodes))
    pctgs[file,names(counts)] <- counts / sum(counts)
    fsom_tmp$MST$size <- sqrt(pctgs[file,] * 500)
    PlotStars(fsom_tmp,
              markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))],
              main = name)

    #PlotPies(UpdateNodeSize(fsom_tmp, reset= TRUE, maxNodeSize = 8),
    #          backgroundValues = fsom$metaclustering,
    #          main = file)

    #t <- table(fsom_tmp$map$mapping[,1])
    #counts[file, names(t)] <- t
  }
  dev.off()
  return(pctgs)
}
