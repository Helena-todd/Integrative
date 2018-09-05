#' compute_pca_distance_scores
#'
#' @param pca_metaclust pca on fsom metaclusters, result of the prcomp function
#' @param samp_rd metadata on donors and recipients
#'
#' @return distances between donors and recipients on the pca's 2 first principal components. A score of one corresponds to the smallest distance ( = maximum resemblance), a score of 0.5 corresponds to the biggest distance.
#' @export
#'
#' @examples
#'
#' compute_pca_distance_scores(pca_metaclust, samp_rd)
#'
compute_pca_distance_scores <- function(pca_metaclust, samp_rd){
  pca_distance_scores <- rep(0, 2)

  for ( i in unique(samp_rd$COUPLENUMBER)){
    couple_names <- rownames(samp_rd[which(samp_rd$COUPLENUMBER==i),])
    pc1_diff <- pca_metaclust$x[couple_names[1],1] - pca_metaclust$x[couple_names[2],1]
    pc2_diff <- pca_metaclust$x[couple_names[1],2] - pca_metaclust$x[couple_names[2],2]

    tot_score <- c(abs(pc1_diff), abs(pc2_diff))
    pca_distance_scores <- rbind( pca_distance_scores, tot_score)
  }

  pca_distance_scores <- pca_distance_scores[-1,]
  pca_distance_scores <- cbind( unique(samp_rd$COUPLENUMBER), pca_distance_scores)
  colnames(pca_distance_scores) <- c("couple_number", "pc1_diff", "pc2_diff")
  rownames(pca_distance_scores) <- paste0("couple_", unique(samp_rd$COUPLENUMBER))

  # re-scale : min distance couple gets a score of 1, max distance couple gets a score of 0.5
  pca_distance_scores[,2] <- scales::rescale(pca_distance_scores[,2],
                                             to = c(0.5, 1))
  pca_distance_scores[,2] <- (1-pca_distance_scores[,2]) + 0.5

  pca_distance_scores[,3] <- scales::rescale(pca_distance_scores[,3],
                                             to = c(0.5, 1))
  pca_distance_scores[,3] <- (1-pca_distance_scores[,3]) + 0.5
}


