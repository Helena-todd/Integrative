#' find_subpatways_driving_PC
#'
#' Identifies which which sub-pathways drive a pca axes of interest
#'
#' @param meta_metabo A dataframe containing the metabolites in columns, and 3 rows: the sub-pathway, the superpathway, and the metabolite name
#' @param pca The result of a PCA, as returned by the ade4::dudi.pca function
#' @param PC A numeric, corresponding to the principal component to be investigated
#'
#' @return A list containing: 1) the "ordered_subpathways", a vector containing all subpathways ordered by importance along the PC of interect
#' and 2) "princ_axes", a dataframe containing both the meta_metabo data and the importance of the different metabolites along each PC of the pca
#' @export
#'
#' @examples
#' pca <- ade4::dudi.pca(metabo_data)
#' result_pc1 <- find_subpatways_driving_PC(meta_metabo = meta_metabo, pca = pca, PC = 1)
#
find_subpatways_driving_PC <- function(meta_metabo, pca, PC){
  metadat <- t(meta_metabo)
  rownames(metadat) <- metadat[,3]
  princ_axes <- merge.data.frame(metadat, pca$c1, by = "row.names")
  sub_pathways <- names(table(princ_axes[,2]))
  sub_path_ranks <- lapply(seq_along(sub_pathways), function(i){
    ranks <- which(princ_axes[order(princ_axes[,4+PC]),2]==sub_pathways[i])
    mean_rank <- sum(ranks)/length(ranks)
    names(mean_rank) <- sub_pathways[i]
    mean_rank
  }) # returns mean rank of each sub_pathway
  ordered_subpathways <- unlist(sub_path_ranks[order(unlist(sub_path_ranks))])
  return(list(ordered_subpathways = ordered_subpathways, princ_axes = princ_axes))
}
