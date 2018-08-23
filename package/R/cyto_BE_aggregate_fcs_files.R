#' BE_aggregate_fcs_files
#'
#' @param patient_names vector containing the fcs names of the patients, and their short names as names
#' @param fsom FlowSOM map to match the fcs files to
#' @param fcs_dir Directory of the FCS files
#' @param markers markers to keep in the final table
#' @param metadata_patients metadata information about the patients. Should contain at least one column "DATEOFCYTOFEXPERIMENT", and have the short names of the patients as rownames.
#'
#' @return a data table containing the expression of the markers of interest in the cells of all patients, with three extracolumns : cluster_id, patient_id and day_id.
#' @export
#'
#' @examples
#'
#' aggreg_table <- BE_aggregate_fcs_files(patient_names = recip_names, fsom = fsom_recip,
#' fcs_dir = fcs_dir, markers = functional_marks,
#' metadata_patients = samp_recip)
BE_aggregate_fcs_files <- function(patient_names, fsom, fcs_dir, markers, metadata_patients){

  aggreg_table <- purrr::map_df(seq_along(patient_names), function(i){
    file <- patient_names[[i]]
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

    data.table::data.table(
      ff_t@exprs[, which(colnames(ff_t@exprs)%in%names(markers))],
      file_id = rep(names(patient_names)[i], nrow(ff_t)),
      cluster_id = fsom_tmp$map$mapping[,1],
      day_id = rep(metadata_patients[names(patient_names)[i],"DATEOFCYTOFEXPERIMENT"],nrow(ff_t))
    )
  })
  return(aggreg_table)
}
