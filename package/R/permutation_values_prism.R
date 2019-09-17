perm_vals_PRISM <- function(norm_data, PRISM_name){
  handle3 <- qsub_lapply(
    X = seq_along(colnames(norm_data)[-which(colnames(norm_data) %in% c("group", "age_recip", "gender_comp"))]),
    qsub_environment = c("norm_data", "perm_val_LR_demo"),
    qsub_packages = c("tidyverse", "BioGVHD", "magrittr", "purrr", "dplyr", "nnet", "pROC"),
    qsub_config = override_qsub_config(
      memory = "15G",
      max_wall_time = "24:00:00",
      name = "BioGVHD",

      wait = FALSE,
      remove_tmp_folder = FALSE,
      stop_on_error = FALSE
    ),

    FUN = function(feat) {
      time0 <- Sys.time()
      tryCatch({
        PVD <- perm_val_LR_demo(data = norm_data, nb_perm = 1000, ind_feature = feat+1)
      }, error = function(e) {
        tibble(feature=feat, error = e$message) #dataset=dta$id,
      })
      time1 <- Sys.time()
      time_diff <- as.numeric(difftime(time1, time0, units = "secs"))
      #list(method_output = PVD, time_diff = time_diff)
      PVD
    })

  saveRDS(handle3, "handle_louis_donors_TNT.rds")
}


# library(qsub)
#
# readRDS("handle2.rds")
# output2<-qsub_retrieve(handle2)
# save(output2, file="dataset1_5D.RData")
# save(clus_res, file="dataset1_5D_umap.RData")
#
# map(output2, 1)
