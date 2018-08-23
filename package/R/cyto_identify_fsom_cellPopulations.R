#' identify_fsom_cellPopulations
#'
#' Matches cell types to the clusters of a FlowSOM object.
#'
#' @param fsom A FlowSOM object
#' @param prettyMarkerNames vector containing the biological names of the markers and the associated fluorochromes
#' @param cellTypes A list of the markers expressed in each cell type, as returned by the "markers_of_cellTypes" function
#' @param pdf_name Name of the pdf to be generated, which will contain one map per cell type, with the clusters corresponding to the cell type highlighted in red.
#' @param view either "grid" or "MST", given as a parameter to FlowSOM::PlotStars
#'
#' @return Cell types associated to the clusters in the FlowSOM object
#' @export
#'
#' @examples
#'
#' celllabels <- identify_fsom_cellPopulations(fsom = fsom_meta_recip, prettyMarkerNames, cellTypes,
#' pdf_name = "pops_filt_marks_grid.pdf", view="MST")
identify_fsom_cellPopulations <- function(fsom, prettyMarkerNames, cellTypes, pdf_name, view){
  labels <- rep("Unknown", fsom$FlowSOM$map$nNodes)
  if(class(colnames(fsom$FlowSOM$map$medianValues))=="character"){
    colnames(fsom$FlowSOM$map$medianValues)<-prettyMarkerNames[colnames(fsom$FlowSOM$map$medianValues)]
  } else {colnames(fsom$FlowSOM$map$medianValues) <- colnames(fsom$FlowSOM$map$medianValues)}

  names(fsom$FlowSOM$prettyColnames) <- fsom$FlowSOM$prettyColnames
  pdf(pdf_name)
  for(cellType in names(cellTypes)){
    query_res <- QueryStarPlot(UpdateNodeSize(fsom$FlowSOM, reset = TRUE, maxNodeSize = 12),
                               query = cellTypes[[cellType]], main = cellType, view=view)
    labels[query_res$selected] <- cellType
  }
  dev.off()
  return(labels)
}
