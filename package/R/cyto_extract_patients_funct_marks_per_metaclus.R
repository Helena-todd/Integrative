
#' extract_funct_markers
#'
#' Generates a list of matrices (one matrix per metacluster), containing the
#' expression of the functional markers in the cells from that cluster, from
#' all patients
#'
#' @param recip_names vector containing the names of recipients' fcs files
#' @param fsom FlowSOM object
#' @param cols_to_use functional markers for which to extract MFIs per metacluster per patient
#' @param min_ref vector of minimal values to use for rescaling
#' @param max_ref vector of maximal values to use for rescaling
#' @param files2rescale vector of names of the files that need to be rescaled
#' @param samp_patients table providing information on the patients (group, batch, ...)
#'
#' @return a list of MFIs matrices for each column to use
#' @export
#'
#' @examples
#' funct_matrices <- extract_funct_markers(recip_names, fsom, fcs_dir, cols_to_use,
#'                                         min_ref, max_ref, files2rescale = c("D1073"))
#'
extract_funct_markers <- function(recip_names, fsom, cols_to_use,
                               min_ref, max_ref, files2rescale,
                               samp_patients){

  samp_pat <- samp_patients[names(recip_names),]

  metaclust_values <- list()

  # read in a patient's file and transform it accurately:
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

    if(names(recip_names[i]) %in% files2rescale){
      for (marker in colnames(exprs(ff_t))[c(3,17,28:62,71)]){
        exprs(ff_t)[, marker] <-
          scales::rescale(exprs(ff_t)[, marker],
                          to = c(min_ref[marker], max_ref[marker]))
      }
    }

    fsom_tmp <- FlowSOM::NewData(fsom$FlowSOM, ff_t)
    date_cytof <- samp_pat$DATEOFCYTOFEXPERIMENT[i]
    group_patient <- samp_pat$GROUP[i]

    for (j in seq_along(levels(factor(fsom$metaclustering)))){
      metaclust_cells <- subset(fsom_tmp$data,
                                fsom$metaclustering[fsom_tmp$map$mapping[, 1]] == j)
      # if(is.null(dim(metaclust_values[j]))){
      if(length(unlist(metaclust_values[j]))==0){
        metaclust_values[[j]] <- data.frame(metaclust_cells[,which(colnames(metaclust_cells)%in%names(funct_marks))]) %>%
          mutate(patients = rep(i, nrow(metaclust_cells)),
                 batch = rep(date_cytof, nrow(metaclust_cells)),
                 group = rep(group_patient, nrow(metaclust_cells)))
      } else if(nrow(metaclust_cells)>1){
        metaclust_values_tmp <- data.frame(metaclust_cells[,which(colnames(metaclust_cells)%in%names(funct_marks))]) %>%
          mutate(patients = rep(i, nrow(metaclust_cells)),
                 batch = rep(date_cytof, nrow(metaclust_cells)),
                 group = rep(group_patient, nrow(metaclust_cells)))
        metaclust_values[[j]] <- rbind(metaclust_values[[j]], metaclust_values_tmp)
      } else if(nrow(metaclust_cells)==1){
        metaclust_values_tmp <- as.data.frame(t(metaclust_cells[,which(colnames(metaclust_cells)%in%names(funct_marks))])) %>%
          mutate( patients = i, batch = date_cytof, group = group_patient)
        metaclust_values[[j]] <- rbind(metaclust_values[[j]], metaclust_values_tmp)
      }
    }


  }
  return(metaclust_values)
}

