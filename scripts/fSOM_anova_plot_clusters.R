multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





#Then the actual code:
  
  
  
fsom_test2 <- fsom
fsom_test2$FlowSOM$map$codes <- t(sapply(1:225, function(i){means[,i]/max(pctgs[,i])}))
fsom_test2$FlowSOM$map$medianValues <- t(sapply(1:225, function(i){means[,i]/max(pctgs[,i])}))

#fsom_test2$FlowSOM$map$codes <- t(apply(means, 2, function(x){x / max(x)}))
#fsom_test2$FlowSOM$map$medianValues <- t(apply(means, 2, function(x){x / max(x)}))
fsom_test2$FlowSOM$prettyColnames <- rownames(means)

pdf("Test_anova.pdf", width = 20, height = 20)
PlotNumbers(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 0), view = "grid")
PlotStars(UpdateNodeSize(fsom$FlowSOM, reset =TRUE), view = "grid")
PlotStars(UpdateNodeSize(fsom_test2$FlowSOM, reset =TRUE),
          1:3, view = "grid",
          backgroundValues = factor(p_v < 0.05, levels = c(FALSE, TRUE)),
          backgroundColor = c("#FFFFFF00", "#FF000055"))
plots <- list()
for(i in 1:225){
  data_tmp <- data.frame(Var = pctgs[, i],
                         Visit = group_res$group)
  plots[[i]] <- ggplot(data_tmp) + 
    geom_boxplot(aes(x = Visit, y = Var), #binwidth = 0.01,
                 col = c("black", "red")[1 + (p_v[i] < 0.05)]) + 
    ggtitle(paste0(i," (p = ",round(p_v[i], 3),")"))+
    theme_minimal()
  
  # ggplot(data_tmp) + 
  #   geom_dotplot(aes(x = Visit, y = Var), binaxis = "y",stackdir="center",
  #                col = c("black", "red")[1 + (p_v[i] < 0.05)]) + 
  #   ggtitle(paste0(i," (p = ",round(p_v[i], 3),")"))+
  #   theme_minimal()
}
multiplot(plotlist = plots, 
          layout = matrix(1:225, nrow = 15, byrow = TRUE)[15:1,])
dev.off()
