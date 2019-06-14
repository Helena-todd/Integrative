PlotStars_pctgPos <- function (fsom,
                               pctg_pos,
                               view = "MST",
                               colorPalette = grDevices::colorRampPalette(
                                 c("#00007F",
                                   "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00",
                                   "red", "#7F0000")),
                               starBg = "white",
                               backgroundValues = NULL,
                               backgroundColor = function(n) {
                                 grDevices::rainbow(n, alpha = 0.3)
                               },
                               backgroundLim = NULL, backgroundBreaks = NULL, backgroundSize = NULL,
                               thresholds = NULL, legend = TRUE, query = NULL, main = "")
{
  add.vertex.shape("star", clip = igraph.shape.noclip, plot = FlowSOM:::mystar,
                   parameters = list(vertex.data = NULL, vertex.cP = colorPalette,
                                     vertex.scale = TRUE, vertex.bg = starBg))

  data <- pctg_pos
  scale <- FALSE

  switch(view, MST = {
    layout <- fsom$MST$l
    lty <- 1
  }, grid = {
    layout <- as.matrix(fsom$map$grid)
    lty <- 0
  }, tSNE = {
    layout <- fsom$MST$l2
    lty <- 0
  }, stop("The view should be MST, grid or tSNE. tSNE will only work\n                   if you specified this when building the MST."))
  if (!is.null(backgroundValues)) {
    background <- FlowSOM:::computeBackgroundColor(backgroundValues,
                                         backgroundColor, backgroundLim, backgroundBreaks)
    if (is.null(backgroundSize)) {
      backgroundSize <- fsom$MST$size
      backgroundSize[backgroundSize == 0] <- 3
    }
  }
  else {
    background <- NULL
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(1, 1, 1, 1))
  if (legend) {
    if (!is.null(backgroundValues)) {
      graphics::layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE),
                       widths = c(1, 2), heights = c(1))
    }
    else {
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE),
                       widths = c(1, 2), heights = c(1))
    }
    if (is.null(query)) {
      FlowSOM:::plotStarLegend(colnames(pctg_pos), colorPalette(ncol(data)))
    }
    else {
      FlowSOM:::plotStarQuery(colnames(pctg_pos), values = query ==
                      "high", colorPalette(ncol(data)))
    }
    if (!is.null(backgroundValues)) {
      FlowSOM:::PlotBackgroundLegend(backgroundValues, background)
    }
  }
  igraph::plot.igraph(fsom$MST$g, vertex.shape = "star", vertex.label = NA,
                      vertex.size = fsom$MST$size, vertex.data = data, vertex.cP = colorPalette(ncol(data)),
                      vertex.scale = scale, layout = layout, edge.lty = lty,
                      mark.groups = background$groups, mark.col = background$col[background$values],
                      mark.border = background$col[background$values], mark.expand = backgroundSize,
                      main = main)
  graphics::par(oldpar)
  graphics::layout(1)
}
