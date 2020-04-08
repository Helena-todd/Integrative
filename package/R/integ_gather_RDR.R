#' Gather Recipients and D-R
#'
#' This function takes as imput: information on the patients under the form of
#' a couple_info dataframe, as well as the path to the matrices of the selected
#' features for recipients, donors and couples.
#' optionally, it can also take the lists of genes that had the same behaviour
#' in the two cohorts and select only these features.
#' It returns a big dataframe, regrouping the selected features of recipients,
#' donors and couples, as well as some meta info.
#'
#' @param data_couples The path to the matrix with couples in rows and selected feaatures in columns
#' @param beh_couples The path to the list of couple features that had the same behaviour in the 2 cohorts
#' @param data_recip The path to the matrix with recip in rows and selected feaatures in columns
#' @param beh_recip The path to the list of recip features that had the same behaviour in the 2 cohorts
#' @param couple_info The dataframe containing all donors and recipients in rows and meta information on them
#' @param couple_info_2keep The columns in couple info that you wish to keep
#'
#' @return a dataframe containing the groupes features of D, R, D-R and some meta info
#' @export

gather_RDR <- function(data_couples, beh_couples = NULL,
                        data_recip, beh_recip = NULL,
                        couple_info, couple_info_2keep){
  load(data_couples)
  if(!is.null(beh_couples)){
    # select only the features that had the same behaviour in the two cohorts
    sel_beh <- read.xlsx(beh_couples)
    pctgs_sel_ft_couples_90 <- pctgs_sel_ft_couples_90[,c("group", "COUPLENUMBER", sel_beh$sel_beh)]
  }
  # add "couple_" in front of the feature names
  colnames2keep <- which(colnames(pctgs_sel_ft_couples_90) %in% c("group", "COUPLENUMBER"))
  colnames(pctgs_sel_ft_couples_90)[-colnames2keep] <-
    gsub("X", "couple_", colnames(pctgs_sel_ft_couples_90)[-colnames2keep])
  # combine feature info with couple info about gender_comp...
  pctgs_couples <- couple_info %>%
    rownames_to_column("Id.Cryostem.R") %>%
    dplyr::select(couple_info_2keep) %>%
    dplyr::filter(grepl("R", rownames(couple_info))) %>%
    inner_join(pctgs_sel_ft_couples_90, by = "COUPLENUMBER")

  # recip info:
  pctgs_sel_ft_recip_90 <- readRDS(data_recip)
  if(!is.null(beh_recip)){
    # select only the features that had the same behaviour in the two cohorts
    sel_beh <- read.xlsx(beh_recip)
    pctgs_sel_ft_recip_90 <- pctgs_sel_ft_recip_90[,c("group", sel_beh$features)]
  }
  colnames2keep <- which(colnames(pctgs_sel_ft_recip_90) %in% c("group"))
  colnames(pctgs_sel_ft_recip_90)[-colnames2keep] <-
    gsub("X", "R_", colnames(pctgs_sel_ft_recip_90)[-colnames2keep])
  pctgs_recip <- couple_info %>%
    rownames_to_column("Id.Cryostem.R") %>%
    dplyr::select("Id.Cryostem.R", "COUPLENUMBER") %>%
    right_join(pctgs_sel_ft_recip_90 %>% rownames_to_column("Id.Cryostem.R"), by = "Id.Cryostem.R")


  # combine :
  pctgs_all <- pctgs_couples %>%
    left_join(pctgs_recip, by = c("Id.Cryostem.R", "COUPLENUMBER", "group"))
}
