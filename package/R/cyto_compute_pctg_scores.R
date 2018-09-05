#' compute_pctg_scores
#'
#' @param pctgs_meta pca on fsom metaclusters, result of the prcomp function
#' @param samp_rd metadata on donors and recipients
#'
#' @return distances between donors and recipients on the pca's 2 first principal components. A score of one corresponds to the smallest distance ( = maximum resemblance), a score of 0.5 corresponds to the biggest distance.
#' @export
#'
#' @examples
#'
#' compute_pctg_scores(pctgs_meta, samp_rd)
#'
compute_pctg_scores <- function(pctgs_meta, samp_rd){
  pctg_scores <- rep(0, ncol(pctgs_meta))

  for ( i in unique(samp_rd$COUPLENUMBER)){
    couple_names <- rownames(samp_rd[which(samp_rd$COUPLENUMBER==i),])
    couple_pctgs <- pctgs_meta[couple_names,]
    pctg_diff <- abs(couple_pctgs[1,]-couple_pctgs[2,])

    pctg_scores <- rbind( pctg_scores, pctg_diff)
  }

  pctg_scores <- pctg_scores[-1,]
  pctg_scores <- cbind( unique(samp_rd$COUPLENUMBER), pctg_scores)
  rownames(pctg_scores) <- paste0("couple_", unique(samp_rd$COUPLENUMBER))

  # re-scale : min difference between couples gets a score of 1, max difference between couples gets a score of 0.5
  rescaled <- apply(pctg_scores[,2:ncol(pctg_scores)],2, function(x){
    scales::rescale(x, to = c(0.5, 1))
    })
  final_score <- apply(rescaled[,1:ncol(rescaled)],2, function(x){
    (1-x) + 0.5
  })
  return(final_score)
}


