#' Select_features
#'
#' Select features that had a value above a cetain threshold in the permutaion
#' values results.
#'
#' @param perm_vals The permutation distributions results
#' @param norm_data The normalised data with patients in rows and features in columns
#' @param alone_or_demo 1: the permutation values for genes alone, 2: the permutation values for models incorporating the genes + age + gender
#' @param threshold The threshold to use
#' @param file_path The file_path
#' @param file_name The file_name (typically: perm_val_TNT_95_thresh.xlsx)
#'
#' @return A vector containing the selected genes
#' @export
#'
select_features <- function(perm_vals,
                            norm_data,
                            alone_or_demo = 1,
                            threshold,
                            file_path,
                            file_name){
  qt <- unlist(map(perm_vals, alone_or_demo))
  selected_ft_qt <- colnames(norm_data)[(which(qt>threshold))]
  selected_ft_qt <- selected_ft_qt[!is.na(selected_ft_qt)]

  write.xlsx(selected_ft_qt,
             file = paste0(file_path, file_name))
  return(selected_ft_qt)
}
