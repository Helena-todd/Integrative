#' visualise_features_in_boxplots
#'
#' @param list_common_features The list containing the common features and the local and other cohort dataframes, as returned by the function find_common_features
#'
#' @return Boxplots
#' @export

visualise_features_in_boxplots <- function(list_common_features){
  for(i in list_common_features$common_features){
    cryo_tmp <- common_features$local_df %>%
      select(i, group) %>% magrittr::set_colnames(c("Var", "group"))

    louis_tmp <- common_features$other_df %>%
      select(i, group) %>% magrittr::set_colnames(c("Var", "group"))

    ymin <- min(c(cryo_tmp$Var, louis_tmp$Var))
    ymax <- max(c(cryo_tmp$Var, louis_tmp$Var))

    g1<- ggplot(cryo_tmp) + geom_boxplot(aes(x = group, y = Var, col = group)) +
      geom_jitter(aes(x = group, y = Var, col = group)) + theme_minimal() +
      coord_cartesian(ylim = c(ymin, ymax)) +
      ggtitle(paste0(i," in Cryostem"))+
      theme(plot.title = element_text(size = 10, face = "bold"))

    g2<- ggplot(louis_tmp) + geom_boxplot(aes(x = group, y = Var, col = group)) +
      geom_jitter(aes(x = group, y = Var, col = group)) + theme_minimal() +
      coord_cartesian(ylim = c(ymin, ymax)) +
      ggtitle(paste0(i," in St Louis")) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    gmix <- patchwork::wrap_plots(g1, g2)
    print(gmix)
  }
}
