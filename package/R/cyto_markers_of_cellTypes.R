#' markers_of_cellTypes
#'
#'Extracts the necessary information from a table describing which markers are expressed or not in each cell population.
#'Returns a list of all populations, with the markers that allow to identify them.
#'
#' @param marker_table A dataframe containing the populations in rows and the markers in columns. Values are either "+","-" or "NA".
#'
#' @return A list of markers describing each cell population.
#' @export
#'
#' @examples
#'
#' cellTypes <- markers_of_cellTypes(marker_table = pop_mark)
markers_of_cellTypes <- function(marker_table){
  # replace + and - by high and low in all the table
  dam2 = reshape2::dcast(
    dplyr::mutate(
      reshape2::melt(marker_table,id.var="Markers"),
      value=plyr::mapvalues(
        value, c("+","-","High","Low"),c("high","low","high","low"))
    ),Markers~variable)

  ## define phenotype of cell types with full table (no NAs)
  # cellTypes <- lapply(seq_along(dam2[,1]), function(i){
  #   x <- dam2[i,]
  #   unlist(x[-1])
  # })
  # names(cellTypes) <- dam2[,1]
  #
  # cellTypes <- lapply(cellTypes, sort)

  # define phenotype of cell types with sparse table (NAs ++)
  cellTypes <- lapply(seq_along(dam2[,2]), function(i){
    dam2[i,which(!is.na(dam2[i,]))[-1]]
  })
  names(cellTypes) <- dam2[,1]
  return(cellTypes)
}
