#' tree_func
#'
#' @param final_model result of a random forest
#'
#' @return plots a tree
#' @export
#'
#' @examples
#' tree_func(rf_result)
tree_func <- function(final_model){#,
                      #tree_num) {

  # get tree by index
  tree <- randomForest::getTree(final_model,
                                #k = tree_num,
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))

  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))

  # convert to graph and delete the last node that we don't want to plot
  graph <- igraph::graph_from_data_frame(graph_frame) %>%
    igraph::delete_vertices("0")

  # set node labels
  igraph::V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  igraph::V(graph)$leaf_label <- as.character(tree$prediction)
  igraph::V(graph)$split <- as.character(round(tree$`split point`, digits = 2))

  # plot
  plot <- ggraph::ggraph(graph, 'dendrogram') +
    theme_bw() +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point() +
    ggraph::geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    ggraph::geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    ggraph::geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE,
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))

  print(plot)
}
#We can now plot, e.g. the tree with the smalles number of nodes:

#tree_num <- which(model_rf$finalModel$forest$ndbigtree == min(model_rf$finalModel$forest$ndbigtree))

#tree_func(final_model = model_rf$finalModel, tree_num)

