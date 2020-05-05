#' plot_correlations
#'
#' @param ft_df The dataframe with samples in rows, features + group in columns
#' @param compar The comparison: either "TNT" or "PS"
#'
#' @return A plot and the graph
#' @export
#'
plot_correlations <- function(ft_df,
                              ft_file_other_cohort,
                              compar,
                              pval_threshold){
  correlations <- ft_df %>%
    dplyr::select(-group) %>%
    as.matrix %>%
    cor %>%
    as.data.frame %>%
    rownames_to_column(var = 'var1') %>%
    gather(var2, value, -var1) %>%
    dplyr::filter(var1 < var2)

  # cor.tests to define which correlations are significant
  cor_sel <- correlations[which(correlations[,3] > 0),]
  pvals <- rep(1, nrow(cor_sel))
  for(i in 1: nrow(cor_sel)){
    var1 <- ft_df[,cor_sel[i,1]]
    var2 <- ft_df[,cor_sel[i,2]]
    pvals[i] <- cor.test(var1, var2, method = "spearman", exact = FALSE)$p.value
  }

  cor_sel <- cor_sel[which(pvals < pval_threshold),]

  if(grepl(".RDS", ft_file_other_cohort)){
    other <- readRDS(ft_file_other_cohort) %>%
      dplyr::select(colnames(ft_df))
  } else if(grepl(".RData", ft_file_other_cohort)){
    tmp_name <- load(ft_file_other_cohort)
    other <- eval(parse(text = tmp_name)) %>%
      dplyr::select(colnames(ft_df))
  }
  #
  # load(ft_file_other_cohort)
  # other <- rna_recip_cryo_TNT_90 %>%
  #   dplyr::select(colnames(ft_df))

  correlations_other <- other %>%
    dplyr::select(-group) %>%
    as.matrix %>%
    cor %>%
    as.data.frame %>%
    rownames_to_column(var = 'var1') %>%
    gather(var2, value, -var1) %>%
    dplyr::filter(var1 < var2)

  # cor.tests to define which correlations are significant
  cor_sel_other <- correlations_other[which(correlations_other[,3] > 0),]
  pvals_other <- rep(1, nrow(cor_sel_other))
  for(i in 1: nrow(cor_sel_other)){
    var1 <- other[,cor_sel_other[i,1]]
    var2 <- other[,cor_sel_other[i,2]]
    pvals_other[i] <- cor.test(var1, var2, method = "spearman", exact = FALSE)$p.value
  }

  cor_sel_other <- cor_sel_other[which(pvals_other < pval_threshold),]

  # correlated_features <- correlations %>%
  #   dplyr::filter(value > 0.7)
  # correlated_other <- correlations_other %>%
  #   dplyr::filter(value > 0.7)

  # correlated_both <- intersect(correlated_features[,1:2], correlated_other[,1:2])
  correlated_both <- dplyr::intersect(cor_sel[,1:2], cor_sel_other[,1:2])
  cor <- rep(0, nrow(correlated_both))
  for(i in 1:nrow(correlated_both)){
    cor1 <- cor_sel[which((cor_sel[,1] == correlated_both[i,1]) &
                                        (cor_sel[,2] == correlated_both[i,2])), 3]
    cor2 <- cor_sel_other[which((cor_sel_other[,1] == correlated_both[i,1]) &
                                        (cor_sel_other[,2] == correlated_both[i,2])), 3]
    cor[i] <- mean(cor1, cor2)
  }
  correlated_both$value <- cor

  gr <- graph_from_data_frame(correlated_both, directed = FALSE,
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
  plot(gr, vertex.size = 5, vertex.label.cex = .6, edge.width = (E(gr)$value - round(min(E(gr)$value), digits = 2)) * 10,
       vertex.color = color_nodes, layout = layout_nicely)
  legend(1, 1, legend=c(round(min(E(gr)$value), digits = 3), round(max(E(gr)$value), digits = 3)),
         col=c("grey", "grey"),
         lwd=c((min(E(gr)$value) - round(min(E(gr)$value), digits = 2)) * 10,
               (max(E(gr)$value) - round(min(E(gr)$value), digits = 2)) * 10),
         cex=0.8,
         title="Edge weight", text.font=4, bg='white')
  return(gr)
}

