load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/fsom_meta.RData")


## recompute MFIs, pctgs, p_v, means and mean_norms for metaclusters --------

metacluster_MFIs <- t(apply(fsom$FlowSOM$data[, colsToUse], 2, function(x){
  tapply(x, fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]], median)
}))

metacluster_pctgs <- t(apply(pctgs, 1, function(x){ 
  tapply(x, fsom$metaclustering, sum) 
}))

p_v_meta <- rep(NA, ncol(metacluster_pctgs))
for (i in seq_len(ncol(metacluster_pctgs))) {
  data_tmp <- data.frame(Var = metacluster_pctgs[, i],Visit = group_res$groups)
#  library(ggplot2)
#  ggplot(data_tmp) + geom_boxplot(aes(x = Visit, y = Var)) + theme_minimal()
  fit <- aov(Var ~ Visit, data = data_tmp)
  p_v_meta[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
}

means_meta <- apply(metacluster_pctgs, 2, function(x){tapply(x, group_res$groups, mean)})
means_norm_meta <- (means - min(means))/(max(means) - min(means))




###########################################################
###############   plot metacluster SOM   ##################
###########################################################

fsom_meta <- fsom
fsom_meta$FlowSOM$map$codes <- t(metacluster_MFIs)
fsom_meta$FlowSOM$map$medianValues <- t(metacluster_MFIs)
colnames(fsom_meta$FlowSOM$map$grid)
fsom_meta$FlowSOM$map$colsUsed <- seq_along(rownames(metacluster_MFIs))
fsom_meta$FlowSOM <- BuildMST(fsom_meta$FlowSOM)
fsom_meta$FlowSOM$map$grid<-as.data.frame(cbind(rep(1:6,5), rep(c(1:5),each=6)))

# Size depending on 1 file
#fsom_meta$FlowSOM$MST$size <- metacluster_pctgs[1,]*100 
# Size equal for all nodes
fsom_meta$FlowSOM$MST$size <- rep(15, ncol(metacluster_MFIs))
fsom_meta$FlowSOM$prettyColnames <- fsom$FlowSOM$prettyColnames[fsom$FlowSOM$map$colsUsed]
names(fsom_meta$FlowSOM$prettyColnames) <- fsom$FlowSOM$prettyColnames[fsom$FlowSOM$map$colsUsed]

colnames(fsom_meta$FlowSOM$map$medianValues) <- prettyMarkerNames[colnames(fsom_meta$FlowSOM$map$medianValues)]
PlotStars(fsom_meta$FlowSOM, 
          markers = c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"),
          view="grid")
PlotStars(fsom_meta$FlowSOM,
          markers = c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))

PlotNumbers(fsom_meta$FlowSOM)
save(fsom_meta, file = "fsom_meta.RData")

# plot size proportional to number of cells:
meta_size <- lapply(seq_along(metaclusters), function(i){
  nrow(metaclusters[[i]])
}) %>% unlist()

fsom_meta$FlowSOM$MST$size <- sqrt(meta_size) / max(sqrt(meta_size)) * 15
PlotStars(fsom_meta$FlowSOM, 
          markers = c("CD11a","CD16","CD27","CD3","CD4","CD45RA","CD8a","HLADR","CD19",
                      "CD38","CD161","CCR5"),
          backgroundValues = my_labels)

# generate fcs files containing one metacluster each:
clusters <- lapply(seq_along(fsom$metaclustering), function(i){
  clust <- fsom$FlowSOM$data[which(fsom$FlowSOM$map$mapping[,1]==i),]
})

metaclusters <- lapply(seq_along(1:length(table(fsom$metaclustering))), function(i){
  meta_clust <- do.call("rbind",clusters[which(fsom$metaclustering==i)])
  meta_ff<-new("flowFrame",exprs=meta_clust)
  write.FCS(meta_ff, filename = paste0("metacluster",i,".fcs"))
  meta_clust
})



## plot Anova results for metaclusters -------------------------

fsom_test3 <- fsom_test
fsom_test3$FlowSOM$map$codes <- t(sapply(1:30, function(i){means_meta[,i]/max(metacluster_pctgs[,i])}))
fsom_test3$FlowSOM$map$medianValues <- t(sapply(1:30, function(i){means_meta[,i]/max(metacluster_pctgs[,i])}))
fsom_test3$FlowSOM$prettyColnames <- rownames(means_meta)



pdf("Test_anova_meta.pdf", width = 20, height = 20)
PlotNumbers(UpdateNodeSize(fsom_test$FlowSOM, maxNodeSize = 0), view = "grid")
PlotStars(fsom_test$FlowSOM, view = "grid")#,
          #markers=c("CD4","CD8a","CD20","IgM","CD38","CD25","CD3","CD11a","CD19"))
PlotStars(UpdateNodeSize(fsom_test3$FlowSOM, reset =TRUE),
          1:3, view = "grid",
          backgroundValues = factor(p_v_meta < 0.05, levels = c(FALSE, TRUE)),
          backgroundColor = c("#FFFFFF00", "#FF000055"))
plots <- list()
for(i in 1:30){
  data_tmp <- data.frame(Var = pctgs_meta[, i],
                         Visit = group_res$group)
  plots[[i]] <- ggplot(data_tmp) + 
    geom_boxplot(aes(x = Visit, y = Var), #binwidth = 0.01,
                 col = c("black", "red")[1 + (p_v_meta[i] < 0.05)]) + 
    ggtitle(paste0(i," (p = ",round(p_v_meta[i], 3),")"))+
    theme_minimal()
  
  # ggplot(data_tmp) + 
  #   geom_dotplot(aes(x = Visit, y = Var), binaxis = "y",stackdir="center",
  #                col = c("black", "red")[1 + (p_v[i] < 0.05)]) + 
  #   ggtitle(paste0(i," (p = ",round(p_v[i], 3),")"))+
  #   theme_minimal()
}
multiplot(plotlist = plots, 
          layout = matrix(1:36, nrow = 6, byrow = TRUE)[6:1,])
dev.off()



