#' find_common_behaviour_features
#'
#' Function to find features that had a similar group behaviour in the two cohorts
#'
#' @param list_common_features The list containing the common features and the local and other cohort dataframes, as returned by the function find_common_features
#' @param file_path The path to the folder where all results should be stored
#' @param common_ft_file_name  the excel file name to save common features
#' @param common_beh_file_name the excel file name to save common features with the same behaviour
#'
#' @return The features that had the same behaviour in both cohorts
#' @export

find_common_behaviour_features <- function(list_common_features,
                                           file_path,
                                           common_ft_file_name,
                                           common_beh_file_name){
  dir_local <- list_common_features$local_df %>%
    group_by(group) %>%
    summarise_all(median)

  overexp <- dir_local %>%
    select(-group) %>%
    apply(2, function(column){
      dir_local$group[which.max(as.matrix(column))]
    })

  local_tags <- as.character(overexp)

  dir_other <- list_common_features$other_df %>%
    group_by(group) %>%
    summarise_all(median)

  overexp <- dir_other %>%
    select(-group) %>%
    apply(2, function(column){
      dir_other$group[which.max(as.matrix(column))]
    })

  other_tags <- as.character(overexp)

  common_beh <- list_common_features$common_features[which(local_tags == other_tags)]

  df1 <- data.frame(features = list_common_features$common_features,
                    dir = local_tags)
  write.xlsx(df1,
             file = paste0(file_path, common_ft_file_name))

  df2 <- data.frame(features = common_beh,
                    dir = local_tags[which(local_tags == other_tags)])

  write.xlsx(df2,
             file = paste0(file_path, common_beh_file_name))

  print(paste0(length(common_beh),
               " features out of the ", length(list_common_features$common_features),
               " common selected features have the same behaviour in both cohorts."))
  return(common_beh)

}
