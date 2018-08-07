library(dplyr)
library(openxlsx)
library(metabolomics)


#################################
########  PREPROCESSING  ########
#################################

# args = file containing metadata info for patients of interest, with metabo names as rownames
#      = and metabo.xlsx file
# output = matrix where metabo and metadat info is merged
merge_metabo_metadata <- function(metadata_file, metabo_file){
  load(metadata_file)
  recip_meta <- sample_subset[which(!is.na(sample_subset$METABONAME)),]
  rownames(recip_meta) <- sample_subset$METABONAME[which(!is.na(sample_subset$METABONAME))]
  
  ## strange numbers after row 390, I only read 80 first rows (corresponding to patients)
  metabo<- read.xlsx(metabo_file, rows = c(1:82), colNames = F, rowNames = T)
  meta_metabo <- metabo[1:3,]  
  data_metabolites <- metabo[which(rownames(metabo)%in%rownames(recip_meta)),]
  colnames(data_metabolites) <- metabo[which(rownames(metabo)=="metabolite"),] 
  return(list(subset_meta=recip_meta, data_metabolites= data_metabolites, meta_metabo=meta_metabo))
}

merged <- merge_metabo_metadata(metadata_file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/sample_recip.RData",
                                metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis_filtered.xlsx")
recip_meta <- merged$subset_meta
meta_metabo <- merged$meta_metabo
data_metabolites <- merged$data_metabolites

##### filtering #####

data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics

colsums <- by(data_metabolites, recip_meta$GROUP, # identify metabolites <50% in each group
              FUN = function(x) {colSums(is.na(x)) >= nrow(x)*0.5})
high_na <- which(colsums[[1]]&colsums[[2]]&colsums[[3]]) # identify metabolites <50% in all groups
data_metabo <- data_metabolites[ , -high_na]
for(i in 1:ncol(data_metabo)){ #replace remaining NA by 1/2 min column value + noise
  x <- data_metabo[,i]
  to_replace <- data_metabo[is.na(x),i]
  with_noise <- jitter(rep(0.5*(min(as.numeric(x[-which(is.na(x))]))), length(to_replace)))
  data_metabo[is.na(x),i] <- with_noise
}

# ignore metabolites with too small variance:
sds <- apply(data_metabo,2,sd)
plot(sort(sds), type = "l", ylim = c(0,10^7))
abline(h=quantile(sds, 0.23), col="red")
data_metabo <- data_metabo[,which(sds>=quantile(sds, 0.23))]
rnames <- rownames(data_metabo)
data_metabo <- apply(data_metabo,2,as.numeric)
rownames(data_metabo) <- rnames

# logtransform and normalise:
mat2use <- merge.data.frame(as.data.frame(recip_meta[,2], row.names = rownames(recip_meta)),
                            data_metabo, by = "row.names") %>%
  tibble::column_to_rownames(var="Row.names")

logdata <- LogTransform(mat2use)
normdata <- Normalise(logdata$output, method = "median")
colnames(normdata$output) <- colnames(logdata$output)
norm_data <- normdata$output

# merge dataframes:
big_mat <- merge.data.frame(normdata$output, recip_meta, by = "row.names") %>% 
  tibble::column_to_rownames("Row.names")

save(norm_data, file="norm_data.RData")
meta_metabo <- meta_metabo[,which(meta_metabo[3,]%in%colnames(norm_data))] # only selected metabolites
save(meta_metabo, file="meta_metabo.RData")
save(big_mat, file="big_mat.RData")
save(recip_meta, file="recip_meta.RData")



################################
##########  ANALYSIS  ##########
################################

load("outputs/data/metabo/recip/norm_data.RData")
load("outputs/data/metabo/recip/meta_metabo.RData")
load("outputs/data/metabo/recip/big_mat.RData")
load("outputs/data/metabo/recip/recip_meta.RData")

##### pca #####
acp <- prcomp(big_mat[,2:(ncol(big_mat)-ncol(recip_meta))])
plot(acp$x, col=as.factor(big_mat$GROUP))
library(ade4)
library(adegraphics)
pca <- dudi.pca(big_mat[,2:(ncol(big_mat)-ncol(recip_meta))], center = T, scale = T)
3
s.label(pca$li)
s.class(pca$li, as.factor(big_mat$GROUP), col=c("red","green","blue"),clabel=.7)
s.arrow(pca$c1, clabel = .2, boxes = F)
s.class(pca$c1, as.factor(as.character(meta_metabo[2,])))
s.class(pca$c1, as.factor(as.character(meta_metabo[1,])))

# which metabolites/ sub-pathways/ super-pathways drive the pca axes?
metadat <- t(meta_metabo)
rownames(metadat) <- metadat[,3]
princ_axes <- merge.data.frame(metadat, pca$c1, by = "row.names")
sub_pathways <- names(table(princ_axes[,2]))
sub_path_ranks <- lapply(seq_along(sub_pathways), function(i){
  ranks <- which(princ_axes[order(princ_axes$CS1, decreasing = T),2]==sub_pathways[i])
  mean_rank <- sum(ranks)/length(ranks)
  names(mean_rank) <- sub_pathways[i]
  mean_rank
}) # returns mean rank of each sub_pathway

head(unlist(sub_path_ranks[order(unlist(sub_path_ranks))])) #petits nombres = associes ++ a l'axe 1 (gr non-tol ++)
tail(unlist(sub_path_ranks[order(unlist(sub_path_ranks))])) #grands nombres = anti-associes a l'axe 1 (gr tol 1 tol 2)
ordered_subpathways_pc1 <- unlist(sub_path_ranks[order(unlist(sub_path_ranks))])
save(ordered_subpathways_pc1, file="ordered_subpathways_pc1.RData")
ordered_metabolites_pc1 <- princ_axes[order(princ_axes$CS1),c(1,2,5)]
save(ordered_metabolites_pc1, file = "ordered_metabolites_pc1.RData")


## only on tol1 and tol2 ##
mydata <- big_mat[which(big_mat$GROUP!="non_tolerant"),] # select only tol 1 and tol 2
mat2use <- norm_data[which(big_mat$GROUP!="non_tolerant"),] 

pca <- dudi.pca(mat2use[,-1], center = T, scale = T)
4
s.label(pca$li, clabel = .7)
s.class(pca$li, as.factor(mat2use$Group), col=c("blue","green"),clabel=.7)
s.class(pca$li[,2:3], as.factor(mat2use$Group), col=c("blue","green"),clabel=.7)

# which metabolites/ sub-pathways/ super-pathways drive the 2nd pca axis (difference between tol1 and tol2 ++)?
princ_axes <- cbind(metadat, pca$c1)
sub_pathways <- names(table(princ_axes[,1]))
sub_path_ranks <- lapply(seq_along(sub_pathways), function(i){
  ranks <- which(princ_axes[order(princ_axes$CS3),1]==sub_pathways[i])
  mean_rank <- sum(ranks)/length(ranks)
  names(mean_rank) <- sub_pathways[i]
  mean_rank
}) # returns mean rank of each sub_pathway

head(unlist(sub_path_ranks[order(unlist(sub_path_ranks))])) #petits nombres = associes ++ a l'axe 3 (gr 1-tol ++)
tail(unlist(sub_path_ranks[order(unlist(sub_path_ranks))])) #grands nombres = anti-associes a l'axe  (gr tol 2)
ordered_metabolites_pc3 <- princ_axes[order(princ_axes$CS3),c(1,2,6)]
ordered_subpathways_pc3 <- unlist(sub_path_ranks[order(unlist(sub_path_ranks))])
save(ordered_metabolites_pc3, file = "ordered_metabolites_pc3.RData")
save(ordered_subpathways_pc3, file="ordered_subpathways_pc3.RData")



### differential analysis

res <- TwoGroup(mat2use)
VolcanoPlot(res$output[,4], res$output[,2], cexlab = 0.6)
MetBoxPlots(mat2use, "N-methylproline",cols = c("blue","green"),main = "N-methylproline")
MetBoxPlots(mat2use, "chiro-inositol",cols = c("blue","green"),main = "chiro-inositol")
MetBoxPlots(mat2use, "1,7-dimethylurate",cols = c("blue","green"),main = "1,7-dimethylurate")
MetBoxPlots(mat2use, "1-methylurate",cols = c("blue","green"),main = "1-methylurate")

TwoGroupPlots(mat2use[,-1], res$output[,1], foldchanges = res$output[,4], pvalues = res$output[,2],
              padjmethod = "BH", fcutoff = log(2), pcutoff = 0.05)










