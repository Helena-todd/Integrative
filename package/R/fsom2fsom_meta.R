#' fsom2fsom_meta
#'
#' This fuction reconstructs a FlowSOM object structure with metaclusters instead of clusters. The resulting FlowSOM object can be plotted with functions from the FlowSOM package such as PlotStars, PlotNumbers...
#'
#' @param fsom FlowSOM object
#' @param colsToUse names of the fluorochromes of the markers to use
#' @param pctgs_patients pctgs of cells of each patient in fsom clusters, matrix of n(patients) * p(clusters)
#' @param plot_size Parameter to set the size of the nodes in fsom map. Can be set to "equal_size", or the node size can be proportional to the number of cells from one patient. In this case, plot_size should be a numeric, corresponding to the desired patient.
#'
#' @return New FlowSOM object, with metaclusters instead of clusters.
#' @export
#'
#' @examples
#'fsom_meta_recip <- fsom2fsom_meta(fsom = fsom_recip, colsToUse = markersToPlot,
#'pctgs_patients = pctgs_recip, plot_size = "equal_size")
#'PlotNumbers(fsom_meta_recip$FlowSOM)
#'
#'
fsom2fsom_meta <- function(fsom, colsToUse, pctgs_patients, plot_size){
  metacluster_MFIs <- t(apply(fsom$FlowSOM$data[, colsToUse], 2, function(x){
    tapply(x, fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]], median)
  }))

  metacluster_pctgs <- t(apply(pctgs_patients, 1, function(x){
    tapply(x, fsom$metaclustering, sum)
  }))

  fsom_meta_recip <- fsom
  fsom_meta_recip$FlowSOM$map$codes <- t(metacluster_MFIs)
  fsom_meta_recip$FlowSOM$map$medianValues <- t(metacluster_MFIs)
  fsom_meta_recip$FlowSOM$map$colsUsed <- seq_along(rownames(metacluster_MFIs))
  fsom_meta_recip$FlowSOM <- BuildMST(fsom_meta_recip$FlowSOM, silent = T)

  if((class(plot_size) == "numeric")||(class(plot_size) == "integer")){
    # Size depending on 1 file
    fsom_meta_recip$FlowSOM$MST$size <- sqrt(metacluster_pctgs[plot_size,]*500)
  } else if (plot_size == "equal_size"){
    # Size equal for all nodes
    fsom_meta_recip$FlowSOM$MST$size <- rep(15, ncol(metacluster_MFIs))
  }
  #fsom_meta_recip$FlowSOM$prettyColnames <- fsom_meta_recip$FlowSOM$prettyColnames[fsom$FlowSOM$map$colsUsed]
  return(fsom_meta_recip)
}
