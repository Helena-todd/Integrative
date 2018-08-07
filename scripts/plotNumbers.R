PlotNumbers <- function (fsom, view = "MST", main = NULL, nodeSize = fsom$MST$size, 
                         backgroundValues = NULL, backgroundColor = function(n) {
                           grDevices::rainbow(n, alpha = 0.3)
                         }, backgroundLim = NULL, backgroundBreaks = NULL) 
{
  switch(view, MST = {
    layout <- fsom$MST$l
    lty <- 1
  }, grid = {
    layout <- as.matrix(fsom$map$grid)
    lty <- 0
  }, tSNE = {
    layout <- fsom$MST$l2
    lty <- 0
  })
  if (!is.null(backgroundValues)) {
    background <- FlowSOM:::computeBackgroundColor(backgroundValues, 
                                         backgroundColor, backgroundLim, backgroundBreaks)
  }
  else {
    background <- NULL
  }
  igraph::plot.igraph(fsom$MST$graph, layout = layout, vertex.size = nodeSize, 
                      vertex.label = seq_len(nrow(fsom$map$codes)), 
                      vertex.label.cex = 0.5,
                      edge.lty = lty, 
                      mark.groups = background$groups, mark.col = background$col[background$values], 
                      mark.border = background$col[background$values])
}
