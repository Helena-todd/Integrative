#' metabo_preprocess
#'
#' @param patient_type A charcter string, either "donors", "recipients" or "donors and recipients"
#' @param data_metabolites The metabo data frame, with patients in rows and metabolites in columns
#' @param meta_metabo A dataframe containing information on the metabolites, as returned by the extract_info_metabo function
#' @param data_metabolites_national The metabo data frame of the unused cohort, to extract common metabolites
#' @param rd_meta The dataframe containing clinical information, as returned by the extract_info_metabo function
#' @param names_patients A character string containing the names of the patients to analyse
#' @param pdf_name Name of the figure containing the metabolites with too many NAs
#' @param pdf_variance_name Name of the figure containing the metabolite variance, and the cutoff applied to filter lowly variable metabolites
#' @param save_results If True, the folder where the results should be saved has to be provided to the "save_res_repository"
#' @param save_res_repository The folder where the results should be saved
#'
#' @return a list of objects : the log data, normalised data, metabo information, and big_matrix containing the metabo + clinical info
#' @export
#'
#' @example
#' metabo_preprocess <- function("donors and recipients", data_metabolites, meta_metabo,
#' data_metabolites_national, rd_meta, names_patients,
#' pdf_name = "rm_metabo.pdf", pdf_variance_name = "variance_cutoff.pdf", save_results = T, save_res_repository = "Desktop/")
metabo_preprocess <- function(patient_type, data_metabolites, meta_metabo,
                              data_metabolites_national, rd_meta, names_patients,
                              pdf_name, pdf_variance_name, save_results, save_res_repository){
  ##### filtering #####

  data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics
  national_metabolites <- data_metabo_national$X1[-c(1,2)]
  # remove metabolites that are not present in both the national and the St_Louis cohort :
  data_metabolites <- data_metabolites[,which(colnames(data_metabolites) %in% national_metabolites)]

  ##### generate figure similar to David's #####

  if(patient_type == "donors and recipients"){
    big_mat <- data_metabolites %>%
      mutate_all(as.numeric) %>%
      mutate(METABONAME = names_patients) %>%
      left_join(rd_meta[,c("GROUP", "COUPLENUMBER", "METABONAME")], by = "METABONAME") %>%
      mutate(status = ifelse(test = grepl("D", METABONAME), yes = "D", no = "R"))
    na_pctgs <- big_mat %>%
      group_by(GROUP, status) %>%
      summarize_all(funs(sum(is.na(.)) / length(.)))
    melted <- melt(na_pctgs, id = 1:2)

    to_rm <- lapply(names(table(melted$variable)), function(metabolite){
      values <- melted$value[which(melted$variable == metabolite)]
      all(values > 0.5)
    })

    names2rm <- as.character(names(table(melted$variable))[which(to_rm == T)])
    plot2rm <- melted[which(melted$variable %in% names2rm),]
    table_pctg_na <- plot2rm %>% arrange(status) %>%
      mutate(variable = paste0(variable, GROUP))

    pdf(pdf_name)
    ggplot(data = table_pctg_na,
           mapping = aes(x = variable, fill = GROUP,
                         y = ifelse(test = status == "D",
                                    yes = -value, no = value))) +
      geom_bar(stat = "identity") +
      scale_y_continuous(labels = abs, limits = max(table_pctg_na$value) * c(-1,1)) +
      labs(y = "pctg of NA") +
      coord_flip()
    dev.off()

  } else {
    big_mat <- data_metabolites[names_patients,] %>%
      mutate_all(as.numeric) %>%
      mutate(METABONAME = names_patients) %>%
      left_join(rd_meta[,c("GROUP", "COUPLENUMBER", "METABONAME")], by = "METABONAME")
    na_pctgs <- big_mat %>%
      group_by(GROUP) %>%
      summarize_all(funs(sum(is.na(.)) / length(.)))
    melted <- melt(na_pctgs, id = 1:2)

    to_rm <- lapply(names(table(melted$variable)), function(metabolite){
      values <- melted$value[which(melted$variable == metabolite)]
      all(values > 0.5)
    })
    names2rm <- as.character(names(table(melted$variable))[which(to_rm == T)])
    plot2rm <- melted[which(melted$variable %in% names2rm),]
    table_pctg_na <- plot2rm %>%
      mutate(variable = paste0(variable, GROUP))
    pdf(pdf_name)
    g <- ggplot(data = table_pctg_na,
           mapping = aes(x = variable, fill = GROUP,
                         y = value)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(labels = abs, limits = max(table_pctg_na$value) * c(0,1)) +
      labs(y = "pctg of NA") +
      coord_flip()
    print(g)
    dev.off()
  }

  ##### new filtering : > 50% na in all groups, donors and recipients

  high_na <- which(to_rm==T)
  data_metabo <- data_metabolites[names_patients, , -high_na]

  for(i in 1:ncol(data_metabo)){ #replace remaining NA by 1/2 min column value + noise
    x <- data_metabo[,i]
    to_replace <- data_metabo[is.na(x),i]
    with_noise <- jitter(rep(0.5*(min(as.numeric(x[-which(is.na(x))]))), length(to_replace)))
    data_metabo[is.na(x),i] <- with_noise
  }

  # ignore metabolites with too small variance:
  sds <- apply(data_metabo,2,sd)
  pdf(file = pdf_variance_name)
  plot(sort(sds), type = "l", ylim = c(0,10^7))
  abline(h=quantile(sds, 0.3), col="red")
  dev.off()
  data_metabo <- data_metabo[,which(sds>=quantile(sds, 0.3))]
  rnames <- rownames(data_metabo)
  data_metabo <- apply(data_metabo,2,as.numeric)
  rownames(data_metabo) <- rnames

  # logtransform and normalise:
  mat2use <- merge.data.frame(as.data.frame(rd_meta[,2], row.names = rownames(rd_meta)),
                              data_metabo, by = "row.names") %>%
    tibble::column_to_rownames(var="Row.names")

  logdata <- LogTransform(mat2use)
  normdata <- Normalise(logdata$output, method = "median")
  colnames(normdata$output) <- colnames(logdata$output)
  norm_data <- normdata$output

  # merge dataframes:
  big_mat <- merge.data.frame(normdata$output, rd_meta, by = "row.names") %>%
    tibble::column_to_rownames("Row.names")

  logdata <- logdata$output
  meta_metabo <- meta_metabo[,which(as.character(meta_metabo[3,])%in%colnames(norm_data))] # only selected metabolites

  if (save_results == TRUE){
    save(logdata, file=paste0(save_res_repository,"logdata.RData"))
    save(norm_data, file=paste0(save_res_repository, "norm_data.RData"))
    save(meta_metabo, file=paste0(save_res_repository, "meta_metabo.RData"))
    save(big_mat, file=paste0(save_res_repository, "big_mat.RData"))
  }
  return(tibble::lst(logdata = logdata, norm_data = norm_data,
                     meta_metabo = meta_metabo, big_mat = big_mat))
}
