#' Title
#'
#' @param aggreg_table Table containing all cells from patients of interest, as returned by the function BE_aggregate_to_table
#' @param fsom FlowSOM object
#'
#' @return A dataframe of n(patients)*[p(markers)*m(metaclusters)], containing the mean expression of each functional marker in each metacluster per patient
#' @export
#'
#' @examples
#'
#'aggreg_table <- BE_aggregate_fcs_files(patient_names = recip_names[1:2], fsom = fsom_recip,
#' fcs_dir = fcs_dir, markers = functional_marks,
#' metadata_patients = samp_recip)
#'
#'mean_MFI_df <- BE_compute_MFI_per_metaclust(aggreg_table, fsom)
BE_compute_MFI_per_metaclust <- function(aggreg_table, fsom){
  metaclst <- fsom$metaclustering

  MFI_table <- purrr::map_dfc(seq_along(table(metaclst)), function(i){
    metaclst_cells <- which(aggreg_table$cluster_id %in% which(metaclst==i))
    if (length(metaclst_cells) != 0 ){
      meta_table <- aggreg_table[metaclst_cells,] %>%
        group_by(file_id) %>%
        summarise_at(colnames(aggreg_table)[1:10], mean) %>%
        as.data.frame() %>%
        column_to_rownames("file_id")
      colnames(meta_table) <- paste0("meta",i,"_",colnames(meta_table))
      if(nrow(meta_table) != length(table(aggreg_table$file_id))){
        pat_names <- as.character(names(table(aggreg_table$file_id)))
        extra <- pat_names[which(!pat_names%in%rownames(meta_table))]
        extra_mat <- matrix(rep(0, length(extra)*10), nrow = length(extra)) %>%
          as.data.frame()
        colnames(extra_mat) <- colnames(meta_table)
        rownames(extra_mat) <- extra
        meta_table <- rbind(meta_table, extra_mat)
        meta_table <- meta_table[names(table(aggreg_table$file_id)),]
      }
    } else {
      meta_table <- as.data.frame(rep(0, length(table(aggreg_table$file_id))))
      colnames(meta_table) <- paste0("meta",i)
      rownames(meta_table) <- names(table(aggreg_table$file_id))
    }
    meta_table
  })
  rownames(MFI_table) <- names(table(aggreg_table$file_id))
  MFI_table
}
