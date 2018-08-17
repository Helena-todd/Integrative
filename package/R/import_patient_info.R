#' import_patient_info
#'
#' Imports an excel table containing all the metadata concerning patients
#' Columns conataining NAs are removed
#' Columns containing characters are turned into factors
#' Columns containing dates are imported as such
#'
#' @param data_synthesis_file The excel file name
#' @param patient_names The FCS names of the patients for whom information should be kept
#'
#' @return A data frame corresponding to the excel file, with only the patients of interest and no NAs
#' @export
#'
#' @examples
#'
#' mytable <- import_patient_info(data_synthesis_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx",
#' patient_names = recip_names)

import_patient_info <- function(data_synthesis_file, patient_names){
  samples <- read.xlsx(data_synthesis_file,
            check.names = FALSE) %>%
    dplyr::filter(!is.na(FCSNAME))
  rownames(samples)<- samples$Id.Cryostem.R

  samp_patients <- samples[which(rownames(samples) %in% names(patient_names)),] %>%
    select_if(~!any(is.na(.))) %>%
    mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
           DOB = as.Date(DOB, "%d.%m.%Y"),
           DOG = as.Date(DOG, "%d.%m.%Y"),
           DATEOFSAMPLE = as.Date(DATEOFSAMPLE, "%d.%m.%Y"),
           lastfollowup = as.Date(lastfollowup, "%d.%m.%Y"),
           GROUP = tolower(GROUP)) %>%
    mutate_if(is.character, as.factor)
  rownames(samp_patients)<- samp_patients$Id.Cryostem.R
  return(samp_patients)
}
