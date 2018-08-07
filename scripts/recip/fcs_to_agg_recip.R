## plot expression of pheno markers in patients (QC):
plot_aggregate_markers <- function(samp_recip, prettyMarkerNames, pheno_marks, png_name ){
  samp_recip<- samp_recip %>% dplyr::arrange(DATEOFCYTOFEXPERIMENT) # reorder files by date of experiment
  markersToPlot <- names(prettyMarkerNames)[which(prettyMarkerNames%in%pheno_marks)]
  png(file.path(paste0(png_name)),
      width = 9000,
      height = 4500)
  par(cex.lab = 2.5, mar = c(4.1,5.1,2.1,2.1))
  layout(matrix(1:35, nrow=7, byrow = TRUE))
  gr_levels<-levels(as.factor(samp_recip$DATEOFCYTOFEXPERIMENT))
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
           col = scales::hue_pal()(length(gr_levels))[factor(samp_recip$DATEOFCYTOFEXPERIMENT[ff_agg_recip@exprs[,"File"]],levels=gr_levels)])
  }
  plot.new()
  legend("center",legend=gr_levels, col=scales::hue_pal()(7),
         pch=19, cex=3)
  dev.off()

}
