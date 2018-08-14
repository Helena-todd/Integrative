#' extract_info_metabo
#'
#' @param metadata_file dataframe containing patients in rows and some information in columns
#' @param metabo_file excel file containing the metabolomics data
#'
#' @return A list, containing patient metadata, the metabolites data and the pathways and subpathways of these metabolites.
#' @export
#'
#' @examples
#'
#' merged <- merge_metabo_metadata(metadata_file = samp_recip,
#'  metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis_filtered.xlsx")

extract_info_metabo <- function(metadata_file, metabo_file){
  patient_meta <- metadata_file[which(!is.na(metadata_file$METABONAME)),]
  rownames(patient_meta) <- patient_meta$METABONAME

  ## strange numbers after row 390, I only read 80 first rows (corresponding to patients)
  metabo<- read.xlsx(metabo_file, rows = c(1:82), colNames = F, rowNames = T)
  meta_metabo <- metabo[1:3,]
  data_metabolites <- metabo[which(rownames(metabo)%in%rownames(patient_meta)),]
  colnames(data_metabolites) <- metabo[which(rownames(metabo)=="metabolite"),]
  return(list(subset_meta=patient_meta, data_metabolites= data_metabolites, meta_metabo=meta_metabo))
}
