#' metabo_ plot correlation between 2 cohorts
#'
#' @param model_stat1 The model stat of the cryostem cohort
#' @param norm_data1 The norm_data of the cryostem cohort
#' @param model_stat2 The model stat of the st louis cohort
#' @param norm_data2 The norm_data of the st louis cohort
#' @param stat_of_interest 2 = global AUC, 3 = tol 1 vs non, 4 = tol 2 vs non, 5 = tol 1 vs tol 2
#' @param meta_metabo the information on the pathways and subpathways
#'
#' @return plots
#' @export

metabo_cor_2_cohorts_LR <- function (model_stat1, norm_data1,
                                  model_stat2, norm_data2,
                                  stat_of_interest,
                                  meta_metabo){
  load(model_stat1)
  load(norm_data1)
  qt_cryo <- purrr::map(perm_vals, stat_of_interest)
  names(qt_cryo) <- colnames(norm_data)[-which(colnames(norm_data) %in% c("group", "gender_comp", "age_recip"))]
  load(model_stat2)
  load(norm_data2)
  qt_louis <- purrr::map(perm_vals, stat_of_interest)
  names(qt_louis) <- colnames(norm_data)[-which(colnames(norm_data) %in% c("group", "gender_comp", "age_recip"))]

  print(paste0("Number of metabolites in common between the two cohorts = ",length(which(names(qt_cryo) %in% names(qt_louis)))))
  g_r_cryo <- qt_cryo[which(names(qt_cryo) %in% names(qt_louis))]
  g_r_louis <- qt_louis[names(g_r_cryo)]

  plot(g_r_louis, g_r_cryo, #xlim = c(0.4,0.8), ylim = c(0.4,0.8),
       xlab = "quantiles St Louis",
       ylab = "quantiles Cryostem",
       main = "Metabolite importance in the 2 cohorts")
  abline(a = 0, b=1, col = "red", lty = 2)

  colnames(meta_metabo) <- colnames(norm_data)[-which(colnames(norm_data) %in% c("group", "gender_comp", "age_recip"))]
  meta_metabo_filt <- meta_metabo[,make.names(names(g_r_cryo))]

  for (i in 1:length(table((as.factor(meta_metabo_filt[2,]))))){
    plot(g_r_louis, g_r_cryo,
         col = c("black", "red")[as.factor(meta_metabo_filt[2,]==names(table(as.factor(meta_metabo_filt[2,])))[i])],
         pch = 19,
         xlab = "quantiles St Louis",
         ylab = "quantiles Cryostem",
         main = paste0(names(table(as.factor(meta_metabo_filt[2,])))[i], " importance in the 2 cohorts"))
    abline(a = 0, b=1, col = "red", lty = 2)
  }
}
