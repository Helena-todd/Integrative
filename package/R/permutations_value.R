#' Permutations value
#'
#' This function takes as input a dataset and a feature id and returns the quantile value of the prediction score of the feature compared to the permutation distribution
#'
#' @param data A normalised dataset (dataframe with patients in rows and features in columns)
#' @param nb_perm The number of permutations to generate the permutation distribution
#' @param ind_feature the feature index, a column number
#'
#' @return the quantile value of the feature in the permutation distribution
#' @export

perm_val <- function (data, nb_perm, ind_feature){
  print(paste0("Computing feature ", ind_feature))
  quantInv <- function(distr, value) ecdf(distr)(value)
  data_tmp <- data
  res_perm <- c()
  for(i in 1:(nb_perm-1)){
    indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data_tmp[,ind_feature] <- data[indices, ind_feature]
    fit <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data_tmp,trace=F)
    res_tmp <- pROC::multiclass.roc(data$group,fitted(fit))$auc[[1]]
    res_perm <- c(res_perm, res_tmp)
  }
  fit <- nnet::multinom(group ~ eval(parse(text=colnames(data)[ind_feature])), data=data,trace=F)
  x <- pROC::multiclass.roc(data$group,fitted(fit))$auc[[1]]
  quantInv(res_perm, x)
}


