#' compute_compatibility_scores
#'
#' @param samp_rd metadata table containing the information on donors and recipients
#'
#' @return the differences between donors and recipients
#' @export
#'
#' @examples
#'
#' compatibility_scores <- compute_compatibility_score(samp_rd)
compute_compatibility_score <- function(samp_rd){
  compatibility_scores <- rep(0, 4)

  for ( i in unique(samp_rd$COUPLENUMBER)){
    couple <- samp_rd[which(samp_rd$COUPLENUMBER==i),]
    # gender score : HH & FF > HF > FH
    if (couple$GENDER[2]==couple$GENDER[1]){
      gender_score <- 1
    } else if (couple$GENDER[2] == "M"){
      gender_score <- 0.75
    } else { gender_score <- 0.5 }

    # age score : youngyoung > oy & yo > oo
    med_age <- median(samp_rd$DOB)
    if ((couple$DOB[2]>med_age)&(couple$DOB[1]>med_age)){
      age_score <- 1
    } else if ((couple$DOB[2]<med_age)&(couple$DOB[1]<med_age)){
      age_score <- 0.5
    } else { age_score <- 0.75 }

    # CMV score : -- > ++ > -+ & +-
    if ((couple$CMVStatus[2]==0)&(couple$CMVStatus[1]==0)){
      cmv_score <- 1
    } else if ((couple$CMVStatus[2]==1)&(couple$CMVStatus[1]==1)){
      cmv_score <- 0.75
    } else { cmv_score <- 0.5 }

    # blood_score : AA & BB & OO > 0A & OB > AB & BA (not interested in rhesus)
    if (couple$GROUPE[2] == couple$GROUPE[1]){
      groupe_score <- 1
    } else { groupe_score <- 0.5 }

    tot_score <- c(gender_score, age_score, cmv_score, groupe_score)
    compatibility_scores <- rbind( compatibility_scores, tot_score)
  }

  compatibility_scores <- compatibility_scores[-1,]
  compatibility_scores <- cbind( unique(samp_rd$COUPLENUMBER), compatibility_scores)
  colnames(compatibility_scores) <- c("couple_number", "gender_score", "age_score", "cmv_score", "groupe_score")
  rownames(compatibility_scores) <- paste0("couple_", unique(samp_rd$COUPLENUMBER))
  return(compatibility_scores)
}


