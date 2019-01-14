
#' generate_meta_MFIs
#'
#' Generates a matrix of MFIs of each patient's cells mapped to fsom metaclusters
#'
#' @param recip_names vector containing the names of recipients' fcs files
#' @param fsom FlowSOM object
#' @param cols_to_use functional markers for which to extract MFIs per metacluster per patient
#'
#' @return a list of MFIs matrices for each column to use
#' @export
#'
#' @examples
#' pctgs <- generate_pctgs(recip_names, fsom, pdf_name = "my_pdf.pdf", fcs_dir)
generate_meta_MFIs <- function(recip_names, fsom, cols_to_use){
  MFIs <- lapply(cols_to_use, function(x){
    matrix(
      NA,
      nrow = length(recip_names),
      ncol = max(fsom$metaclustering),
      dimnames = list(
        recip_names,
        paste0("MC", 1:max(fsom$metaclustering))))
  })
  names(MFIs) <- cols_to_use

  #i <- 1
  for (i in seq_along(recip_names)){
    file <- recip_names[[i]]
    file_in<-file.path(fcs_dir, file)
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

    mfis_tmp <- FlowSOM::MetaclusterMFIs(list(FlowSOM = fsom_tmp,
                                              metaclustering = factor(fsom$metaclustering)))
    for(col in cols_to_use){
      MFIs[[col]][file, ] <- mfis_tmp[, col]
    }

  }
  return(MFIs)
}

