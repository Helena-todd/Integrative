#' metabo_ plot correlation between 2 cohorts
#'
#' @param model_stat1 The model stat of the cryostem cohort
#' @param norm_data1 The norm_data of the cryostem cohort
#' @param model_stat2 The model stat of the st louis cohort
#' @param norm_data2 The norm_data of the st louis cohort
#' @param stat_of_interest 2 = global AUC, 3 = tol 1 vs non, 4 = tol 2 vs non, 5 = tol 1 vs tol 2
#' @param meta_metabo the information on the pathways and subpathways
#' @param perm_distribution Boolean defining if the stats to be compared between the two cohorts are AUCs or permutation quantiles
#'
#' @return plots
#' @export

metabo_cor_2_cohorts <- function (model_stat1, norm_data1,
                                  model_stat2, norm_data2,
                                  stat_of_interest,
                                  perm_distribution = FALSE,
                                  meta_metabo){
  load(model_stat1)
  load(norm_data1)
  if(perm_distribution == FALSE){
    g_aucs_r_cryo <- purrr::map(model_stats, stat_of_interest)
  } else {g_aucs_r_cryo <- perm_vals}
  names(g_aucs_r_cryo) <- colnames(norm_data)[-1]
  load(model_stat2)
  load(norm_data2)
  if(perm_distribution == FALSE){
    g_aucs_r_louis <- purrr::map(model_stats, stat_of_interest)
  } else {g_aucs_r_louis <- perm_vals}
  names(g_aucs_r_louis) <- colnames(norm_data)[-1]

  print(paste0("Number of metabolites in common between the two cohorts = ",length(which(names(g_aucs_r_cryo) %in% names(g_aucs_r_louis)))))
  g_r_cryo <- g_aucs_r_cryo[which(names(g_aucs_r_cryo) %in% names(g_aucs_r_louis))]
  g_r_louis <- g_aucs_r_louis[names(g_r_cryo)]

  plot(g_r_louis, g_r_cryo, #xlim = c(0.4,0.8), ylim = c(0.4,0.8),
       xlab = "global AUCs St Louis",
       ylab = "global AUCs Cryostem",
       main = "Metabolite importance in the 2 cohorts")
  abline(a = 0, b=1, col = "red", lty = 2)

  meta_metabo_filt <- meta_metabo[,make.names(names(g_r_cryo))]

  for (i in 1:length(table((as.factor(meta_metabo_filt[2,]))))){
    if(perm_distribution==FALSE){
      plot(g_r_louis, g_r_cryo,
           col = c("black", "red")[as.factor(meta_metabo_filt[2,]==names(table(as.factor(meta_metabo_filt[2,])))[i])],
           pch = 19,
           xlab = "global AUCs St Louis",
           ylab = "global AUCs Cryostem",
           main = paste0(names(table(as.factor(meta_metabo_filt[2,])))[i], " importance in the 2 cohorts"))
      abline(a = 0, b=1, col = "red", lty = 2)
    } else {
      plot(g_r_louis, g_r_cryo,
           col = c("black", "red")[as.factor(meta_metabo_filt[2,]==names(table(as.factor(meta_metabo_filt[2,])))[i])],
           pch = 19,
           xlab = "global AUC permutation distribution quantiles St Louis",
           ylab = "global AUC permutation distribution quantiles Cryostem",
           main = paste0(names(table(as.factor(meta_metabo_filt[2,])))[i], " importance in the 2 cohorts"))
      abline(a = 0, b=1, col = "red", lty = 2)
    }
  }
}
