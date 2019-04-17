#' rf_tol_non_tol
#'
#' @param df_orig Original dataframe on which rf should be run. Should contain a "group" column, with "primary_tolerant", "secondary_tolerant" or "non_tolerant" factor levels
#' @param mtry Number of parameters to be tested by RF at every split
#' @param ntree Number of trees in the RF
#' @param package_2use The random forest implementation to use (either "randomForest" or "ranger")
#'
#' @return RF results on classification of tol vs non-tol patients
#' @export
#'
rf_tol_non_tol <- function(df_orig, mtry = 3, ntree = 10000, package_2use = "randomForest"){
  gr_tmp <- as.character(df_orig$group)
  gr_tmp[which(gr_tmp %in% c("primary_tolerant", "secondary_tolerant"))] <- "tolerant"
  df_tmp <- df_orig
  df_tmp$group <- factor(gr_tmp)
  if(package_2use == "randomForest"){
    rf_res <- randomForest::randomForest(group~., df_tmp, mtry = mtry, ntree = ntree,
                                         importance = T, proximity = T)
  } else {
    rf_res <- ranger::ranger(group~., df_tmp, mtry = mtry, num.trees = ntree,
                             importance = "impurity")
  }

  return(rf_res)
}
