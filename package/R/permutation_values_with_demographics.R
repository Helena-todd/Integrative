#' Permutations value with demographics
#'
#' This function takes as input a dataset and a feature id and returns the quantile value of the prediction score of the feature compared to the permutation distribution
#'
#' @param data A normalised dataset (dataframe with patients in rows and features in columns)
#' @param nb_perm The number of permutations to generate the permutation distribution
#' @param ind_feature the feature index, a column number
#' @param orig_AUC the original AUC, computed only on the demographics
#'
#' @return the quantile value of the feature in the permutation distribution
#' @export

perm_val_demo <- function (data, nb_perm, ind_feature, orig_AUC){
  print(paste0("Computing feature ", ind_feature))
  quantInv <- function(distr, value) ecdf(distr)(value)
  data_tmp <- data
  res_perm <- list(AUC_demo = c(), AUC = c())
  for(i in 1:(nb_perm-1)){
    indices <- sample(1:nrow(data), size = nrow(data), replace = FALSE) # was TRUE before but sometimes generated bugs when it resultd in only 2 classes instaed of 3
    data_tmp[,"group"] <- data[indices, "group"]
    # if(length(levels(as.factor(as.character(data_tmp[,"group"]))))!=3){
    #   browser()
    # }
    #
    # data_tmp$group <- as.factor(as.character(data_tmp$group))
    fit1 <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data_tmp, trace=F)
    fit2 <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])) + gender + age, data=data_tmp,trace=F)
    res_tmp1 <- pROC::multiclass.roc(data_tmp$group,fitted(fit1))$auc[[1]]
    res_tmp2 <- pROC::multiclass.roc(data_tmp$group,fitted(fit2))$auc[[1]]
    diff_demo <- res_tmp2 - orig_AUC
    res_perm[[1]] <- c(res_perm[[1]], res_tmp1)
    res_perm[[2]] <- c(res_perm[[2]], diff_demo)
  }
  fit <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data,trace=F)
  data_tmp[,"group"] <- data[, "group"]
  fit_demo <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])) + gender + age, data=data_tmp,trace=F)
  x <- pROC::multiclass.roc(data$group,fitted(fit))$auc[[1]]
  x_demo <- pROC::multiclass.roc(data$group,fitted(fit_demo))$auc[[1]]
  diff_demo <- x_demo - orig_AUC
  qt <- quantInv(res_perm[[1]], x)
  qt_demo <- quantInv(res_perm[[2]], diff_demo)

  return(tibble::lst(quantile = qt, quantile_demo = qt_demo))
}


