#' export_DEstats_and_graphinfo
#'
#' @param gr A graph as returned by the function plot_correlations
#' @param common_features The features that are common to both cohorts
#' @param common_beh_features The features that have the same behaviour in both cohorts
#' @param norm_data The normalised data
#' @param comparison The comparison to perform: either "TNT" (non tolerant vs tolerant) or "PTST" (primary vs secondary tolerant)
#' @param folder_path The path to the folder where files should be exported
#' @param common_ft_name The name of the file that will contain information on the common features
#' @param common_beh_ft_name The name of the file that will contain information on the common beh features
#'
#' @return excel files
#' @export

export_DEstats_and_graphinfo <- function(gr,
                                         common_features,
                                         common_beh_features,
                                         norm_data,
                                         comparison,
                                         folder_path,
                                         common_ft_name,
                                         common_beh_ft_name){
  # Identify connected components in the graph:
  gr_components <- decompose(gr, mode = "strong", max.comps = NA,
                             min.vertices = 0)
  comp_sizes <- sapply(gr_components, diameter)

  # Generate a new annotated excel table
  annot_genes <- as.data.frame(common_features)

  # Compute the FC and pvalues:
  TS <- tolower(norm_data$group)

  TS <- factor(TS, levels=unique(TS))
  design <- model.matrix(~0+TS)
  colnames(design) <- levels(TS)
  #design
  fit <- lmFit(t(norm_data[,common_features]), design)

  #### Create contrast matrix
  if(comparison == "TNT"){
    cont.matrix <- makeContrasts(group1=non_tolerant-tolerant,
                                 group2=tolerant-non_tolerant,
                                 levels=design)
  } else if(comparison == "PTST"){
    cont.matrix <- makeContrasts(group1=primary_tolerant-secondary_tolerant,
                                 group2=secondary_tolerant-primary_tolerant,
                                 levels=design)
  } else {
    print("Unknown comparison, please enter another comparison type")
  }

  #cont.matrix
  fit2 = contrasts.fit(fit, cont.matrix)

  #### Moderate T-test with 'shrunken' se
  fit.eb <- eBayes(fit2)

  getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
    topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
    genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
    genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
    ##Sort genes on logFC
    genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
    genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
    genes_de_sorted<-rbind(genes_up, genes_down)

    return(genes_de_sorted)
  }

  allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
  DEgenesGroup1<-getDEgenes(allGenesGroup1,1,0)

  DEgenesGroup1 <- DEgenesGroup1[common_features,]
  annot_genes <- cbind(annot_genes, DEgenesGroup1)

  # add the information on correlated genes to the table

  for( i in which(comp_sizes != 0)){
    vertex_names <- vertex_attr(gr_components[[i]], index = V(gr_components[[i]]))[[1]]
    v_names <- c(vertex_names, rep(0, nrow(annot_genes) - length(vertex_names)))
    annot_genes <- cbind(annot_genes, v_names)
  }

  write.xlsx(annot_genes, file = paste0(folder_path, common_ft_name))

  DegenesGroup1_beh <- annot_genes[common_beh_features, 2:7]
  write.xlsx(DegenesGroup1_beh, file = paste0(folder_path, common_beh_ft_name), row.names = TRUE)
}


