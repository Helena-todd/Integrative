
#' generate_donor_pctgs
#'
#' Generates a matrix of pctgs and pdf fsom files of pctgs of each patient's cells mapped to fsom clusters
#'
#' @param donor_names vector containing the names of donors' fcs files
#' @param fsom FlowSOM object
#' @param pdf_name choose a name for the exported pdf which will contain one fsom map per patient
#' @param fcs_dir path to the directory containing fcs files, corresponding to donor_names
#' @param output_dir
#'
#' @return a pctgs matrix and pdf containing one map per patient
#' @export
#'
#' @examples
#' pctgs <- generate_donor_pctgs(donor_names, fsom, pdf_name = "my_pdf.pdf", fcs_dir)
generate_donor_pctgs <- function(donor_names, fsom, pdf_name, fcs_dir, output_dir){
  pctgs <- matrix(
    0,
    length(donor_names),
    ncol = fsom$FlowSOM$map$nNodes,
    dimnames = list(
      donor_names,
      as.character(1:fsom$FlowSOM$map$nNodes)))
  #i <- 1
  pdf(file = pdf_name)

  normal_names <- donor_names[-which(names(donor_names) %in% c("D1073", "D1502", "D2031"))]
  ref_names <- donor_names[c("D2031")]
  toscale_names <- donor_names[c("D1073", "D1502")]

  file <- ref_names[[1]]
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

  min_ref <- apply(ff_t@exprs[,c(3,17,28:62,71)],2,
                   function(x) quantile(x, 0.001))
  max_ref <- apply(ff_t@exprs[,c(3,17,28:62,71)],2,
                   function(x) quantile(x, 0.999))

  fsom_tmp <- FlowSOM::NewData(fsom$FlowSOM, ff_t)
  name<-names(ref_names)[1]
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


  for (j in seq_along(toscale_names)){
    file <- toscale_names[[j]]
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

    for (marker in colnames(ff_t@exprs)[c(3,17,28:62,71)]){
      ff_t@exprs[, marker] <-
        scales::rescale(ff_t@exprs[, marker],
                        to = c(min_ref[marker], max_ref[marker]))
    }

    fsom_tmp <- FlowSOM::NewData(fsom$FlowSOM, ff_t)
    name<-names(toscale_names)[j]
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


  for (i in seq_along(normal_names)){
    file <- normal_names[[i]]
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
    name<-names(normal_names)[i]
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
