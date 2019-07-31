#' plot_density_dr
#'
#' A function to plot both densities of donors and recipients in the 3 groups,
#' to help compare differences in expression
#'
#' @param samp_rd_tmp the clinical data table, with rownames = rownames of mat_pheno_funct
#' @param selected_group the group of choice ("non_tolerant", "primary_tolerant" or "secondary_tolerant")
#' @param featur the feture to plot, must be in colnames of mat_pheno_funct
#' @param mat_pheno_funct the expression data matrix, must have same rownames as samp_rd_tmp
#'
#' @return a ggplot
#' @export
#'
plot_density_dr <- function(samp_rd_tmp, selected_group, featur, mat_pheno_funct){
  samp_rd_gr <- samp_rd_tmp %>%
    rownames_to_column("rn") %>%
    dplyr::filter(GROUP==selected_group) %>%
    column_to_rownames("rn")
  mat_pf_gr <- mat_pheno_funct[rownames(samp_rd_gr),]
  colnames(mat_pf_gr) <- make.names(colnames(mat_pf_gr))

  d_y_coords <- mat_pf_gr[grep("D", rownames(mat_pf_gr)), featur]
  r_y_coords <- mat_pf_gr[-grep("D", rownames(mat_pf_gr)), featur]

  data <-bind_rows(
    tibble(
      y = d_y_coords,
      x = "donor",
      group = 1:length(d_y_coords)
    ),
    tibble(
      y = r_y_coords,
      x = "recipient",
      group = 1:length(d_y_coords)
    )
  )%>%
    mutate(x = as.factor(x))

  ggplot(data, aes(x=y, fill=x)) +
    geom_density(alpha = 0.4)+
    scale_fill_manual(values=c("#000099", "#009900"))+
    ggtitle(paste0(featur, "_", selected_group))
}


