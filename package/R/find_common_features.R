#' Find_common_features
#'
#' @param selected_features The selected features in this cohort
#' @param ndata The normalised data of this cohort (with patients in rows)
#' @param other_features_file The path to the excel file containing the selected features of the other cohort
#' @param other_norm_data_file The path to the RDS file containing the normalised data of the other cohort (with patients in rows)
#' @param file_path Where new results should be saved
#' @param file_name The file name
#'
#' @return A list containing 1) the common features between the 2 cohorts
#' 2) the dataframe of the common features in the local cohort
#' 3) the dataframe of the common features in the other cohort
#' @export

find_common_features <- function(
  selected_features,
  ndata,
  other_features_file,
  other_norm_data_file,
  file_path,
  file_name){

  sel_other <- read.xlsx(other_features_file, colNames = FALSE)
  sel_common <- which(selected_features %in% sel_other$X1)
  common_features <- selected_features[sel_common]
  print(paste0(length(sel_common), " features were commonly found in both cohorts"))

  local_df <- ndata[,c("group", common_features)]
  ndata_other <- readRDS(other_norm_data_file)
  other_df <- ndata_other[,c("group", common_features)]

  saveRDS(local_df, file = paste0(file_path, file_name))

  return(lst(common_features, local_df, other_df))
}
