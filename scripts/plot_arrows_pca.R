pca2 <- ade4::dudi.pca(mat4pca, nf = 10, scannf = F)

info_pca <- info_complete_patients %>% column_to_rownames("Id.Cryostem.R")
info_pca <- info_pca[rownames(cyt4pca),]
ade4::s.class(pca_cyto$li[,c(1,6)], info_pca$group, col=c("red", "blue"))

keep1 <- order(abs(pca2$c1$CS1), decreasing = T)#[1:20]
# keep2 <- order(abs(pca_cyto$c1$CS2), decreasing = T)[1:5]
# tokeep <- unique(c(keep1, keep2))
tokeep <- keep1
ade4::s.arrow(pca2$c1[tokeep,c(1,4)], boxes = F, clabel = 0.5)

totest <- pca_2plot[,c(348, 366:374, 376, 390, 400:403)]
totest <- apply(totest, 2, function(x){as.numeric(as.factor(x))})
barplot(abs(cor(pca2$li[,1], totest)), horiz = TRUE, las = 2, cex.names = 0.3, xlim = c(0,0.6))

barplot(abs(cor(pca2$li[,4], totest)), horiz = TRUE, las = 2, cex.names = 0.3, xlim = c(0,0.6))
