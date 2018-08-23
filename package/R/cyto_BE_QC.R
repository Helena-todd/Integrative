#BE_QC <- function( aggreg_table, aggreg_rescaled_table, png_name_1, png_names_2){
#' BE_QC
#'
#' @param aggreg_table Table containing all cells from patients of interest, as returned by the function BE_aggregate_to_table
#' @param png_name_1 Name of the output file
#'
#' @return A png file
#' @export
#'
#' @examples
#'
#' aggreg_table <- BE_aggregate_fcs_files(patient_names = recip_names[1:2], fsom = fsom_recip,
#' fcs_dir = fcs_dir, markers = functional_marks,
#' metadata_patients = samp_recip)
#' BE_QC(aggreg_table, "my_output_file.png")
BE_QC <- function( aggreg_table, png_name_1){
  ix <- aggreg_table %>%
    mutate(i = row_number()) %>%
    group_by(file_id) %>%
    sample_n(5000) %>%
    ungroup() %>%
    pull(i) # find indices to use both on aggreg_table and aggreg_rescaled_table

  sampled_aggreg_table <-
    aggreg_table[ix,] %>%
    arrange(day_id) %>%
    mutate(file_id = factor(file_id, levels = unique(file_id)))

  grDevices::png(file.path(paste0(png_name_1)),
      width = 5000,
      height = 2500)

  sat_gathered <- sampled_aggreg_table %>%
    gather(marker, value, -file_id, -cluster_id, -day_id) %>%
    mutate(marker_name = prettyMarkerNames[marker])

  g <- ggplot(sat_gathered) + # %>% dplyr::filter(marker_name == "IL10")) +
    geom_violin(aes(file_id, value, colour = factor(day_id))) +
    labs(x = "", y = "") +
    facet_wrap(~ marker_name)


  #ggsave()

  # ggplots <- lapply(seq_along(functional_marks), function(i){
  #   marker <- functional_marks[i]
  #   mark <- names(prettyMarkerNames[which(prettyMarkerNames==marker)])
  #   ggplot(sampled_aggreg_table) +
  #     geom_point(aes_string("file_id", mark, colour = "factor(day)")) +
  #     labs(x = "", y = marker)
  #   day <- as.factor(sampled_aggreg_table$day_id)
  #   ggplot2::qplot(sampled_aggreg_table$file_id,
  #                  sampled_aggreg_table[,mark],
  #                  col = day,
  #                  xlab = "", ylab = marker#, geom = "jitter"
  #                  )
  #})
  #print(patchwork::wrap_plots(ggplots, ncol = 4))
  print(g)
  dev.off()
}


