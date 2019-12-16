#' plot_correlations
#'
#' @param ft_df The dataframe with samples in rows, features + group in columns
#' @param compar The comparison: either "TNT" or "PS"
#'
#' @return A plot
#' @export
#'
plot_correlations <- function(ft_df, compar){
  correlations <- ft_df %>%
    dplyr::select(-group) %>%
    as.matrix %>%
    cor %>%
    as.data.frame %>%
    rownames_to_column(var = 'var1') %>%
    gather(var2, value, -var1) %>%
    dplyr::filter(var1 < var2)

  correlated_features <- correlations %>%
    dplyr::filter(value > 0.7)
  gr <- graph_from_data_frame(correlated_features, directed = FALSE,
                              vertices=colnames(ft_df)[-which(colnames(ft_df)=="group")])
  nodes <- get.vertex.attribute(gr)

  if(compar == "TNT"){
    color_nodes <- rep("#FF000080", length(nodes$name))

    gr_medians <- ft_df %>%
      group_by(group) %>%
      summarize_if(is.numeric, median) %>%
      as.data.frame()
    rownames(gr_medians) <- gr_medians$group

    medians <- gr_medians[,nodes$name]
    for (gene in seq_along(colnames(medians))){
      if (medians["tolerant", gene] == max(medians[, gene])){
        color_nodes[gene] <- "#0000FF80"
      }
    }
  } else if(compar == "PS"){
    color_nodes <- rep("#0000FF80", length(nodes$name))

    gr_medians <- ft_df %>%
      group_by(group) %>%
      summarize_if(is.numeric, median) %>%
      as.data.frame()
    rownames(gr_medians) <- gr_medians$group

    medians <- gr_medians[,nodes$name]
    for (gene in seq_along(colnames(medians))){
      if (medians["primary_tolerant", gene] == max(medians[, gene])){
        color_nodes[gene] <- "#00FF0080"
      }
    }
  }

  set.seed(1)
  plot(gr, vertex.size = 7, vertex.label.cex = .5, edge.width = (E(gr)$value - .7) * 20,
       vertex.color = color_nodes, layout = layout_nicely)
}

