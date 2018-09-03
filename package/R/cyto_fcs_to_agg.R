#' Title
#'
#' @param fcs_dir param
#' @param fcs_names param
#' @param seed param
#' @param cTotal param
#' @param output_name param
#'
#' @return aggregated flowframe as returned by the FlowSOM::AggregateFlowFrames function
#' @export
#'
#' @examples bla
fcs_to_agg <- function(fcs_dir, fcs_names, seed, cTotal, output_name){
  setwd(dir = fcs_dir)
  set.seed(seed)
  ff_agg <- FlowSOM::AggregateFlowFrames(fcs_names,
                                cTotal= cTotal, writeOutput = TRUE,
                                outputFile = output_name)
  ff_agg <- flowCore::transform(ff_agg,
                      flowCore::transformList(colnames(ff_agg@exprs)[c(3,17,28:62,71)], arcsinhTransform(b=1/5, a=0, c=0)))
  return(ff_agg)
}
