plot_lines_dr <- function(samp_rd_tmp, selected_group, color, featur){
  samp_rd_gr <- samp_rd_tmp %>%
    rownames_to_column("rn") %>%
    dplyr::filter(GROUP==selected_group) %>%
    column_to_rownames("rn")
  mat_pf_gr <- mat_pheno_funct[rownames(samp_rd_gr),]
  colnames(mat_pf_gr) <- make.names(colnames(mat_pf_gr))

  d_y_coords <- mat_pf_gr[grep("D", rownames(mat_pf_gr)), featur]
  r_y_coords <- mat_pf_gr[grep("R", rownames(mat_pf_gr)), featur]

  data <-bind_rows(
    tibble(
      y = d_y_coords,
      x = 1,
      group = 1:length(d_y_coords)
    ),
    tibble(
      y = r_y_coords,
      x = 2,
      group = 1:length(d_y_coords)
    )
  )

  ggplot(data) + geom_line(aes(x, y, group = group), col = color) +
    ggtitle(featur)
}


