#' Title
#'
#' @param df_orig Original dataframe on which rf should be run. Should contain a "group" column, with "primary_tolerant", "secondary_tolerant" or "non_tolerant" factor levels
#' @param mtry Number of parameters to be tested by RF at every split
#' @param ntree Number of trees in the RF
#'
#' @return RF results on classification of tol vs non-tol patients
#' @export

rf_tol1_tol2 <- function(df_orig, mtry = 3, ntree = 10000){
  which(df_orig$group=="non_tolerant")
  df_tol <- df_orig[-which(df_orig$group=="non_tolerant"),]
  df_tol$group <- factor(as.character(df_tol$group))

  rf_res <- randomForest::randomForest(group~.,df_tol, mtry = mtry, ntree = ntree,
                                         importance = T, proximity = T)
  return(rf_res)
}
