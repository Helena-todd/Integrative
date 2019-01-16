#' plot_pie
#'
#' This fuction plots a pie
#'
#' @param df A dataframe, with 2 columns: "group", containing the names of the groups, and "value", containing the associated percentages
#' @param colors The colors to use for every group in the pie
#'
#' @return Plots the pie
#' @export
#'
#' @examples
#' df <- data.frame(
#' group = c("Male", "Female", "Child"),
#' value = c(25, 25, 50))
#'
#' plot_pie(df, colors = c("blue","green","yellow"))
#'
#'
#'
plot_pie <- function(df, colors){
  bp<- ggplot(df, aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0)

  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )

  pie +  blank_theme +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values= colors)
}
