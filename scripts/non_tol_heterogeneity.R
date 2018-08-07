## investing fsom cluster heterogeneity, especially in the non tolerant group


## on fsom map containing 10000 cells from each file -------------
####################

# per file:
recip<-ff_agg_recip
fsom_recip <- NewData(fsom$FlowSOM, recip)
PlotPies(fsom_recip, as.factor(fsom_recip$data[,"File"]), colorPalette = colorRampPalette(rainbow(49)))
#colorPalette = colorRampPalette(c("ivory","indianred1","green","blue")))
groups_4<-sample_recip$GROUP
groups_4[which((sample_recip$aGVHD==1)&(sample_recip$GROUP=="non_tolerant"))]<-"non_tol_aGVHD"
gr_recip<-rep(groups_4,each=10000)
PlotPies(fsom_recip, gr_recip, colorPalette = colorRampPalette(c("lightgreen","blue","indianred","gold")))


# only non_tol files 
recip<-ff_agg_recip
recip@exprs<- ff_agg_recip@exprs[which(ff_agg_recip@exprs[,"File"]%in%which(sample_recip$GROUP=="non_tolerant")),]
fsom_recip <- NewData(fsom$FlowSOM, recip)
PlotPies(fsom_recip, as.factor(fsom_recip$data[,"File"]), colorPalette = colorRampPalette(rainbow(49)))





## on all cells from original files assigned to fsom clusters  -------------
####################

fsom_test2 <- fsom
fsom_test2$FlowSOM$map$codes <- t(sapply(1:225, function(i){means[,i]/max(pctgs[,i])}))
fsom_test2$FlowSOM$map$medianValues <- t(sapply(1:225, function(i){means[,i]/max(pctgs[,i])}))

fsom_test2$FlowSOM$map$codes <- t(apply(means, 2, function(x){x / max(x)}))
fsom_test2$FlowSOM$map$medianValues <- t(apply(means, 2, function(x){x / max(x)}))
fsom_test2$FlowSOM$prettyColnames <- rownames(means)
PlotStars(fsom_test2$FlowSOM,1:3)

fsom_test2$FlowSOM$map$codes <- t(pctgs)
fsom_test2$FlowSOM$map$medianValues <- t(pctgs)
fsom_test2$FlowSOM$prettyColnames <- names(recip_names)

PlotNumbers(UpdateNodeSize(fsom$FlowSOM, maxNodeSize = 0), view = "grid")
PlotStars(UpdateNodeSize(fsom$FlowSOM, reset =TRUE))
PlotStars(fsom_test2$FlowSOM,1:49)


## only on non tolerant :

fsom_test2$FlowSOM$map$codes <- t(pctgs[which(sample_recip$GROUP=="non_tolerant"),])
fsom_test2$FlowSOM$map$medianValues <- t(pctgs[which(sample_recip$GROUP=="non_tolerant"),])
fsom_test2$FlowSOM$prettyColnames <- names(recip_names[which(sample_recip$GROUP=="non_tolerant")])

PlotStars(fsom_test2$FlowSOM,1:27)




###### only on B cell clusters ---------------------
# clusters 70 and 119

c70.c119<-data.frame(c(ff_agg_recip@exprs[which(fsom$FlowSOM$map$mapping[,1]==70),"File"],
                       ff_agg_recip@exprs[which(fsom$FlowSOM$map$mapping[,1]==119),"File"]))
colnames(c70.c119)<-"files"
c70.c119$cluster <- c(rep(70, length(ff_agg_recip@exprs[which(fsom$FlowSOM$map$mapping[,1]==70),"File"])),
                      rep(119, length(ff_agg_recip@exprs[which(fsom$FlowSOM$map$mapping[,1]==119),"File"])))
c70.c119$cell_id <- c(which(fsom$FlowSOM$map$mapping[,1]==70), which(fsom$FlowSOM$map$mapping[,1]==119))

ggplot(c70.c119, aes(x = as.factor(cluster), y = as.factor(files), fill=as.factor(files))) +
  geom_bar(stat='identity')

ggplot(c70.c119[which(c70.c119$cluster==70),], aes(files, files)) + 
  geom_bar(stat='identity')
names(recip_names)[c(26, 30:31, 48)] # 3 non-tol and 1 tol

barplot(fsom$FlowSOM$map$medianValues[70,colsToUse], cex.names = .5)
barplot(fsom$FlowSOM$map$sdValues[70,colsToUse], cex.names = .5)

ggplot(c70.c119[which(c70.c119$cluster==119),], aes(files, files)) + 
  geom_bar(stat='identity')
names(recip_names)[c(30, 47:49)] # 2 non-tol and 2 tol

barplot(fsom$FlowSOM$map$medianValues[119,colsToUse], cex.names = .5)
barplot(fsom$FlowSOM$map$sdValues[119,colsToUse], cex.names = .5)

