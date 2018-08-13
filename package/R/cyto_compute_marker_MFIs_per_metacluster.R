# compute, per metacluster, median expression of funct markers in patients:
#' Title
#'
#' @param ff_agg aggregated cells from each patient, as returned by AggregateFlowFrames FlowSOM function
#' @param pctgs_meta Matrix of n(patients) * p(metaclusters) containing the cell percentages from patients in fsom metaclusters
#' @param fsom FlowSOM object
#' @param functional_marks list of markers to use to compute median MFIs in fsom metaclusters
#'
#' @return description of output
#' @export
compute_marker_MFIs_per_metacluster <- function( ff_agg, pctgs_meta, fsom, functional_marks){
  patient_IDs <- as.factor(ff_agg@exprs[,72])
  pctgs_and_MFIs <- pctgs_meta
  metacluster_ids <- FlowSOM::GetMetaclusters(fsom)
  markersToCluster <- names(functional_marks)
  for(i in 1:length(table(fsom$metaclustering))){
    patient_MFIs <- apply(ff_agg@exprs[which(metacluster_ids==i),markersToCluster], 2, function(markerValues){
      tapply(markerValues, patient_IDs[which(metacluster_ids==i)] , median)
    })
    colnames(patient_MFIs) <- as.character(prettyMarkerNames[which(names(prettyMarkerNames)%in%colnames(patient_MFIs))])
    colnames(patient_MFIs) <- paste0("meta",i,"_",colnames(patient_MFIs))
    pctgs_and_MFIs<-cbind(pctgs_and_MFIs, patient_MFIs)
  }

  pctgs_and_MFIs <- as.data.frame(pctgs_and_MFIs)
  pctgs_and_MFIs[is.na(pctgs_and_MFIs)] <- 0
  return(pctgs_and_MFIs)
}
