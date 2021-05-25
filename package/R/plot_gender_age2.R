#' plots_gender_age
#'
#' Generates forest plots and other informative plots to better understand the
#' role of a feature in tolerance with respect to age and gender
#'
#' @param norm_data A dataframe containing patients in rows, features in columns, and the columns "group", "gender_comp" and "age_recip"
#' @param features A vector of the features of interest
#' @param pdf_name The pdf name
#'
#' @return a PDF
#' @export
#'
plots_gender_age2 <- function(norm_data,
                             features,
                             pdf_name,
                             excel_name,
                             compar){
  pdf(pdf_name, width = 14, height = 10)
  toplot <- norm_data[,c("gender_comp", "age_recip", "group",
                         features)]
  toplot$age_recip <- round(toplot$age_recip/365)

  coeffiscients <- matrix(rep(0, 8 * length(features)),
                          nrow = length(features))

  for(ft in seq_along(features)){
    coeffs <- lapply(new_models[[ft]], function(mod){
      sum1 <- summary(mod)
      beta <- sum1$coefficients["Var", 1]
      se <- sum1$coefficients["Var", 2]
      pvalue <- sum1$coefficients["Var", 4]
      data.frame(beta = beta, se = se, pvalue = pvalue)
    })
    df_coeff <- do.call(rbind, coeffs)
    name <- c("Var", "Var + gender", "Var + age", "Var + gender + age")
    df_c <- cbind(name, df_coeff)

    coeffiscients[ft,] <- c(df_c$beta, df_c$pvalue)

    p1 <- ggforestplot::forestplot(
      df = df_c,
      estimate = beta,
      pvalue = pvalue,
      psignif = 0.05,
      logodds = TRUE,
      xlab = "1-SD increment in BMI\nper 1-SD increment in biomarker concentration",
      title = paste0("Associations of the feature to tolerance")
    )
    i<-ft
    data_tmp_gender <- data.frame(Var = toplot[, i+3],gender = toplot[, 1], group = toplot[,3])
    data_tmp_age <- data.frame(Var = toplot[, i+3],age_recip = toplot[, 2], group = toplot[,3])

    if(compar == "TNT"){
      g1 <- ggplot(data_tmp_gender) + geom_boxplot(aes(x = gender, y = Var)) +
        geom_jitter(aes(x = gender, y = Var, col = group), shape=16, position=position_jitter(0.2)) +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
      g2 <- ggplot(data_tmp_age, aes(x = age_recip, y = Var)) +
        geom_point(aes(x = age_recip, y = Var, col = as.factor(group))) + geom_smooth(method = "loess") +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
      g0 <- ggplot(data_tmp_age, aes(x = as.factor(group), y = Var)) +
        geom_boxplot(aes(x = as.factor(group), y = Var, col = as.factor(group))) +
        geom_jitter(aes(x = group, y = Var, col = group), shape=16, position=position_jitter(0.2)) +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
    } else if (compar == "PTST"){
      g1 <- ggplot(data_tmp_gender) + geom_boxplot(aes(x = gender, y = Var)) +
        geom_jitter(aes(x = gender, y = Var, col = group), shape=16, position=position_jitter(0.2)) +
        scale_color_manual(values = c("green", "blue")) +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
      g2 <- ggplot(data_tmp_age, aes(x = age_recip, y = Var)) +
        geom_point(aes(x = age_recip, y = Var, col = as.factor(group))) + geom_smooth(method = "loess") +
        scale_color_manual(values = c("green", "blue")) +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
      g0 <- ggplot(data_tmp_age, aes(x = as.factor(group), y = Var)) +
        geom_boxplot(aes(x = as.factor(group), y = Var, col = as.factor(group))) +
        scale_color_manual(values = c("green", "blue")) +
        geom_jitter(aes(x = group, y = Var, col = group), shape=16, position=position_jitter(0.2)) +
        theme_minimal() + ggtitle(colnames(toplot)[i+3])
    }

    grid.arrange(p1, g0, g1, g2, nrow = 2)
  }
  dev.off()
  rownames(coeffiscients) <- features
  colnames(coeffiscients) <- c("beta_var", "beta_var_gender", "beta_var_age", "beta_var_gender_age",
                               "pval_var", "pval_var_gender", "pval_var_age", "pval_var_gender_age")
  write.xlsx(coeffiscients, file = excel_name, rowNames = TRUE)
}

