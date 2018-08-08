#' ggplot_analysis_results
#'
#' @param analysis Type of analysis to be performed. For now, either "PCA" or "tSNE" can be run.
#' @param perplexity_nb If tSNE analysis, choose the perplexity parameter. Default is 8.
#' @param data_matrix Matrix on which to run the PCA or tSNE analysis, with patients in rows and parameters in columns.
#' @param metadata Dataframe containing information about the patients, with patients in rows and variables in columns.
#' @param col_by Metadata column to use as a factor to color points in analysis.
#' @param shape_by Metadata column to use as a factor to display the shape of points in analysis.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' ggplot_analysis_results("tSNE", perplexity = 8, data_matrix = pctgs_recip, metadata = samp_recip,
#' col_by = "DATEOFCYTOFEXPERIMENT", shape_by = "GROUP")
#'
#'

ggplot_analysis_results <- function(analysis, perplexity_nb=8, data_matrix, metadata, col_by, shape_by){
  if (analysis=="PCA"){
    analyss <- prcomp(data_matrix)
    df <- metadata %>%
      slice(match(rownames(data_matrix), metadata$Id.Cryostem.R)) %>%
      dplyr::mutate(axis_1 = analyss$x[,1], axis_2 = analyss$x[,2])
  } else if (analysis == "tSNE"){
    set.seed(1)
    analyss <- Rtsne::Rtsne(data_matrix, perplexity = perplexity_nb)
    df <- metadata %>%
      slice(match(rownames(data_matrix), metadata$Id.Cryostem.R)) %>%
      dplyr::mutate(axis_1 = analyss$Y[,1], axis_2 = analyss$Y[,2])
  }

  ggplot(df) +
    geom_point(aes(x = axis_1, y = axis_2,
                   col = as.factor(!!rlang::sym(col_by)),
                   shape = !!rlang::sym(shape_by)), size = 4) +
    geom_text(aes(x = axis_1, y = axis_2, label= Id.Cryostem.R)) +
    theme_minimal() +
    theme(text = element_text(size = 12))
}
