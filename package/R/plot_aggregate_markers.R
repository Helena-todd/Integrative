#' cyto_plot_aggregate_markers
#'
#' Plots expression of phenotypic markers across patients (QC):
#'
#' @param samp_recip Metadata about the patients
#' @param prettyMarkerNames A vector containing the biological marker names as well as the fluorochrome names
#' @param pheno_marks A vector containing phenotypic markers
#' @param png_name Name of the png file to be exported
#' @param ff_agg Aggregated cells from each patient, as returned by the AggregateFlowFrames FlowSOM function
#'
#' @return A png
#' @export
#'
#' @examples bla
plot_aggregate_markers <- function(samp_recip, arrange_by, prettyMarkerNames, pheno_marks, png_name, ff_agg ){
  samp_recip<- samp_recip %>% arrange(!!rlang::sym(arrange_by)) # reorder files by date of experiment
  markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
  png(file.path(paste0(png_name)),
      width = 9000,
      height = 4500)
  par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
  layout(matrix(1:35, nrow=7, byrow = TRUE))
  selected_column <- samp_recip[[arrange_by]]
  gr_levels<-levels(as.factor(selected_column))
  for(marker in markersToPlot){
    print(paste0("Plotting ",prettyMarkerNames[marker]," for the aggregated file."))
    plot(0,#exprs(ff_agg[,c("File_scattered",marker)]),
         type="n",
         xlab = "Files",
         ylab = prettyMarkerNames[marker],
         ylim = c(0, max(10, ff_agg_recip@exprs[,marker])),
         xlim= c(0,length(recip_names)+1),
         xaxt="n")
    abline(v=seq_along(rownames(samp_recip)), col="lightgrey")
    axis(side=1, at=seq_along(rownames(samp_recip)), labels=samp_recip$Id.Cryostem.R, las=2,cex=0.7)
    points(ff_agg_recip@exprs[,c("File_scattered",marker)],
           pch=".",
           col = scales::hue_pal()(length(gr_levels))[factor(selected_column[ff_agg_recip@exprs[,"File"]],levels=gr_levels)])
  }
  plot.new()
  legend("center",legend=gr_levels, col=scales::hue_pal()(length(gr_levels)),
         pch=19, cex=3)
  dev.off()
}
