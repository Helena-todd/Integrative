BE_rescale_per_day <- function(aggreg_table){
  days <- names(table(aggreg_table$day_id))
  min_day1 <- apply(aggreg_table[which(aggreg_table$day_id==days[1]),1:10],2,
                    function(x) quantile(x, 0))
  max_day1 <- apply(aggreg_table[which(aggreg_table$day_id==days[1]),1:10],2,
                    function(x) quantile(x, 1))
  aggreg_table <- as.data.frame(aggreg_table)

  for(i in 2:length(days)){
    for (marker in colnames(aggreg_table)[1:10]){
      # aggreg_table[which(aggreg_table$day_id==days[[i]]), ..marker] <- (((aggreg_table[which(aggreg_table$day_id==days[[i]]), ..marker] -
      #                                                                     min(aggreg_table[which(aggreg_table$day_id==days[[i]]), ..marker])) /
      #                                                                    (max(aggreg_table[which(aggreg_table$day_id==days[[i]]), ..marker]) - min(aggreg_table[which(aggreg_table$day_id==days[[i]]), ..marker]))) *
      #                                                                   (max_day1[marker] - min_day1[marker]) ) + min_day1[marker]

      aggreg_table[which(aggreg_table$day_id==days[[i]]), marker] <- scales::rescale(aggreg_table[which(aggreg_table$day_id==days[[i]]), marker],
                                                                                       to = c(min_day1[marker], max_day1[marker]))
    }
  }
  data.table::as.data.table(aggreg_table)
}

