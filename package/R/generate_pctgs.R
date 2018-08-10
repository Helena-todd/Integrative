
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
generate_pctgs <- function(recip_names, fsom, pdf_name, fcs_dir, output_dir){
  pctgs <- matrix(
    0,
    length(recip_names),
    ncol = fsom$FlowSOM$map$nNodes,
    dimnames = list(
      recip_names,
      as.character(1:fsom$FlowSOM$map$nNodes)))
  #i <- 1
  pdf(file = pdf_name)
  for (i in seq_along(recip_names)){
    file <- recip_names[[i]]
    file_in<-file.path(fcs_dir, file)
    file_out<-file.path(output_dir, file)
    message(file)
    ff <- flowCore::read.FCS(file_in)
    from <- flowCore::colnames(ff)[c(3,17,28:62,71)]
    tlist <- flowCore::transformList(
      from = from,
      tfun = flowCore::arcsinhTransform(b=1/5, a=0, c=0),
      to = from
    )
    ff_t <- flowCore::transform(
      ff,
      tlist
    )
    fsom_tmp <- FlowSOM::NewData(fsom$FlowSOM, ff_t)
    name<-names(recip_names)[i]
    #PlotStars(fsom_tmp,main = name)

    counts <- table(FlowSOM::GetClusters(fsom_tmp))
    #pctgs <- rep(0, fsom_tmp$map$nNodes)
    #names(pctgs) <- as.character(seq_len(fsom_tmp$map$nNodes))
    pctgs[file,names(counts)] <- counts / sum(counts)
    fsom_tmp$MST$size <- sqrt(pctgs[file,] * 500)
    FlowSOM::PlotStars(fsom_tmp,
              markers = names(prettyMarkerNames)[which(prettyMarkerNames%in% c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))],
              main = name)

    #PlotPies(UpdateNodeSize(fsom_tmp, reset= TRUE, maxNodeSize = 8),
    #          backgroundValues = fsom$metaclustering,
    #          main = file)

    #t <- table(fsom_tmp$map$mapping[,1])
    #counts[file, names(t)] <- t

    m <- matrix(0, nrow = nrow(ff), ncol = 3, dimnames = list(NULL,
                                                              c("FlowSOM-clusters",
                                                                "FlowSOM-metaclusters",
                                                                "FlowSOM-metaclusters-jittered")))
    m[, "FlowSOM-clusters"] <- FlowSOM::GetClusters(fsom_tmp)
    m[, "FlowSOM-metaclusters"] <- FlowSOM::GetMetaclusters(fsom_tmp, fsom$metaclustering)
    m[, "FlowSOM-metaclusters-jittered"] <- as.numeric(FlowSOM::GetMetaclusters(fsom_tmp, fsom$metaclustering)) +
      rnorm(nrow(ff), sd = 0.1)

    ff_updated <- flowCore::flowFrame(exprs = cbind(flowCore::exprs(ff), m))
    ff_updated@parameters@data[,"desc"] <- c(ff@parameters@data[,"desc"],
                                             "FlowSOM-clusters",
                                             "FlowSOM-metaclusters",
                                             "FlowSOM-metaclusters-jittered")
    ff_updated@description <- ff@description
    #
    # ff_updated <- cbind2(ff, m)

    flowCore::write.FCS(ff_updated, file_out)
  }
  dev.off()
  return(pctgs)
}
