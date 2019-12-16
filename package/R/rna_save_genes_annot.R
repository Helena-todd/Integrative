#' Save_genes_annot
#' 
#' Function to save the names of the selected genes, as well as a data frame 
#' containing the FC and pvalues for these genes, and the graph informations
#' (ie how the genes were correlated)
#'
#' @param norm_data The normalised data matrix
#' @param path_obj Where to save the vector of selected genes
#' @param selected_genes The selected genes to save
#' @param compar_design Either "TNT" or "PS"
#' @param graph The graph from which correlation information needs to be saved
#' @param path_annot_obj  Where to save the annotation object
#'
#' @return Two excel objects: the vector of selected genes and the dataframe of
#' the annotated selected genes.
#' @export 
#'
 

save_genes_annot <- function(norm_data,
                             path_obj,
                             selected_genes,
                             compar_design,
                             graph = NULL,
                             path_annot_obj
                             ){
  sel_genes <- norm_data[,c("group", selected_genes)]

  write_rds(sel_genes, path = path_obj)
  
  # Generate a new annotated excel table
  annot_genes <- as.data.frame(selected_genes)
  
  if(compar_design == "TNT"){
    patient_groups <- tolower(colData$GROUP)
    patient_groups[which(patient_groups != "non_tolerant")] <- "tolerant"
    TS <- patient_groups[c(donor_id)]
  } else if (compar_design == "PS"){
    TS <- colData[rownames(norm_data),] %>% 
      dplyr::select(GROUP) %>% 
      mutate(GROUP = tolower(GROUP))
    TS <- TS$GROUP
  }
  
  TS <- factor(TS, levels=unique(TS))
  design <- model.matrix(~0+TS)
  colnames(design) <- levels(TS)
  #design
  fit <- lmFit(t(norm_data[,selected_genes]), design)
  
  if(compar_design == "TNT"){
    cont.matrix <- makeContrasts(group1=non_tolerant-tolerant,
                                 group2=tolerant-non_tolerant,
                                 levels=design) 
  } else {
    cont.matrix <- makeContrasts(group1=primary_tolerant-secondary_tolerant,
                                 group2=secondary_tolerant-primary_tolerant,
                                 levels=design)
  }
  
  
  #cont.matrix
  fit2 = contrasts.fit(fit, cont.matrix)
  
  #### Moderate T-test with 'shrunken' se
  fit.eb <- eBayes(fit2)
  
  
  allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
  DEgenesGroup1<-getDEgenes(allGenesGroup1,1,0)
  
  # length(which(selected_genes_qt[sel_common] %in% rownames(DEgenesGroup1)))
  
  # add them to the table
  
  colnames(annot_genes) <- "V1"
  DEgenesGroup1 <- DEgenesGroup1[as.character(annot_genes$V1),] 
  annot_genes <- cbind(annot_genes, DEgenesGroup1)
  
  if(is.igraph(graph) == TRUE){
    # Identify connected components in the graph:
    gr_components <- decompose(graph, mode = "strong", max.comps = NA,
                               min.vertices = 0)
    comp_sizes <- sapply(gr_components, diameter)
    
    for( i in which(comp_sizes != 0)){
      vertex_names <- vertex_attr(gr_components[[i]], index = V(gr_components[[i]]))[[1]]
      v_names <- c(vertex_names, rep(0, nrow(annot_genes) - length(vertex_names)))
      annot_genes <- cbind(annot_genes, v_names)
    }
  }
  write.xlsx(annot_genes, file = path_annot_obj)
}
