#' Permutations value log regression with demographics
#'
#' This function takes as input a dataset and a feature id and returns
#' 1) the quantile value of the prediction score of the feature compared to
#' the permutation distribution
#' 2) the quantile value when taking the age and gender demographics into account
#'
#' @param data A normalised dataset (dataframe with patients in rows and features in columns)
#' @param nb_perm The number of permutations to generate the permutation distribution
#' @param ind_feature the feature index, a column number
#'
#' @return the quantile value of the feature in the permutation distribution
#' @export

perm_val_LR_demo <- function (data, nb_perm, ind_feature){
  print(paste0("Computing feature ", ind_feature))
  quantInv <- function(distr, value) ecdf(distr)(value)
  data_tmp <- data
  res_perm <- list(AUC_demo = c(), AUC = c())
  for(i in 1:(nb_perm-1)){
    indices <- sample(1:nrow(data), size = nrow(data), replace = FALSE) # was TRUE before but sometimes generated bugs when it resultd in only 2 classes instaed of 3
    data_tmp[,"group"] <- data[indices, "group"]

    fit <- glm(group ~ gender_comp + age_recip, data=data_tmp,trace=F, family = "binomial")
    orig_AUC <- as.numeric(pROC::roc(data_tmp$group,fitted(fit))$auc)

    fit1 <- glm(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data_tmp, trace=F, family = "binomial")
    fit2 <- glm(group ~ eval(parse(text=colnames(data)[ind_feature])) + gender_comp + age_recip, data=data_tmp,trace=F, family = "binomial")
    res_tmp1 <- as.numeric(pROC::roc(data_tmp$group,fitted(fit1))$auc)
    res_tmp2 <- as.numeric(pROC::roc(data_tmp$group,fitted(fit2))$auc)
    diff_demo <- res_tmp2 - orig_AUC
    res_perm[[1]] <- c(res_perm[[1]], res_tmp1)
    res_perm[[2]] <- c(res_perm[[2]], diff_demo)
  }
  data_tmp[,"group"] <- data[, "group"]
  fit <- glm(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data,trace=F, family = "binomial")
  fit_demo <- glm(group ~ eval(parse(text=colnames(data)[ind_feature])) + gender_comp + age_recip, data=data_tmp,trace=F, family = "binomial")
  x <- as.numeric(pROC::roc(data$group,fitted(fit))$auc)
  x_demo <- as.numeric(pROC::roc(data$group,fitted(fit_demo))$auc)

  fit_ag <- glm(group ~ gender_comp + age_recip, data=data_tmp,trace=F, family = "binomial")
  orig_AUC <- as.numeric(pROC::roc(data_tmp$group,fitted(fit_ag))$auc)

  diff_demo <- x_demo - orig_AUC
  qt <- quantInv(res_perm[[1]], x)
  qt_demo <- quantInv(res_perm[[2]], diff_demo)

  return(tibble::lst(quantile = qt, quantile_demo = qt_demo))
}


