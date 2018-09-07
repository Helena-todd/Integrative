library(dplyr)
library(openxlsx)
library(metabolomics)
library(adegraphics)


#################################
########  PREPROCESSING  ########
#################################

sample_subset <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018.xlsx",
                     check.names = FALSE) %>%
  mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
         GROUP = tolower(GROUP)) %>%
  dplyr::filter(!is.na(METABONAME)) %>%
  dplyr::filter(grepl("D",Id.Cryostem.R)) %>%
  dplyr::select_if(~sum(!is.na(.)) > 0)
rownames(sample_subset) <- sample_subset$METABONAME

save(sample_subset, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/donors/sample_donors.RData")

info <- extract_info_metabo(metadata_file = sample_subset,
                            metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis_filtered.xlsx")
donor_meta <- info$subset_meta
meta_metabo <- info$meta_metabo
data_metabolites <- info$data_metabolites

##### filtering #####

data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics

colsums <- by(data_metabolites, sample_subset$GROUP, # identify metabolites <50% in each group
              FUN = function(x) {colSums(is.na(x)) >= nrow(x)*0.5})
high_na <- which(colsums[[1]]&colsums[[2]]&colsums[[3]]) # identify metabolites <50% in all groups
data_metabo <- data_metabolites[ , -high_na]
for(i in 1:ncol(data_metabo)){ #replace remaining NA by 1/2 min column value
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
mat2use <- merge.data.frame(as.data.frame(sample_subset[,2], row.names = rownames(sample_subset)),
                            data_metabo, by = "row.names") %>%
  tibble::column_to_rownames(var="Row.names")

logdata <- LogTransform(mat2use)
save(logdata, file = "logdata_metabo_donors.RData")
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
save(recip_meta, file="donors_meta.RData")



################################
##########  ANALYSIS  ##########
################################

load("outputs/data/metabo/donors/norm_data.RData.RData")
load("outputs/data/metabo/donors/meta_metabo.RData")
load("outputs/data/metabo/donors/big_mat.RData")
load("outputs/data/metabo/donors/donors_meta.RData")

##### pca #####
acp <- prcomp(big_mat[,2:(ncol(big_mat)-ncol(recip_meta))])
plot(acp$x, col=as.factor(big_mat$GROUP))
library(ade4)
pca <- dudi.pca(big_mat[,2:(ncol(big_mat)-ncol(recip_meta))], center = T, scale = T)
4
s.label(pca$li, clabel = .7)
s.class(pca$li, as.factor(big_mat$GROUP), col=c("red","green","blue"),clabel=.7)
s.arrow(pca$c1, clabel = .2, boxes = F)
s.class(pca$c1, as.factor(as.character(meta_metabo[2,])), clabel=.5)
s.class(pca$c1, as.factor(as.character(meta_metabo[1,])),clabel = .5 )

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



## only with tolerant 1 and 2 patients ##
mydata <- big_mat[which(big_mat$GROUP!="non_tolerant"),] # select only tol 1 and tol 2
mat2use <- norm_data[which(big_mat$GROUP!="non_tolerant"),]

pca <- dudi.pca(mat2use[,-1], center = T, scale = T)
2
s.label(pca$li, clabel = .7)
s.class(pca$li, as.factor(mat2use$Group), col=c("blue","green"),clabel=.7)


### differential analysis

res <- TwoGroup(mat2use)
VolcanoPlot(res$output[,4], res$output[,2], cexlab = 0.6)
MetBoxPlots(mat2use, "ursodeoxycholate",cols = c("blue","green"),main = "ursodeoxycholate")

TwoGroupPlots(mat2use[,-1], res$output[,1], foldchanges = res$output[,4], pvalues = res$output[,2],
              padjmethod = "BH", fcutoff = log(2), pcutoff = 0.05)



####################################################
####### Turn metabolites -> subpathway table #######
####################################################

load("outputs/data/metabo/donors/logdata_metabo_donors.RData")
load("outputs/data/metabo/donors/meta_metabo.RData")
load("outputs/data/metabo/donors/sample_donors.RData")

data_metabo <- as.data.frame(logdata$output[,-1])
meta_metabo <- meta_metabo[,which(as.character(meta_metabo[3,])%in%colnames(data_metabo))]
data_metabo <- data_metabo[,which(colnames(data_metabo)%in%as.character(meta_metabo[3,]))]

library(data.table)
t_data_metabo <- t(data_metabo)
subpath_info <- as.character(meta_metabo[1,])
data_meta <- as.data.table(t_data_metabo)
binded <- cbind(data_meta, subpath_info)

subpaths <- binded[,lapply(.SD, median), by=subpath_info]
subpaths <- as.data.frame(subpaths) %>% column_to_rownames("subpath_info")
subpaths <- t(subpaths)
save(subpaths, file = "subpaths_table_donors.RData")

### Normalise :
subpaths_norm <- apply(subpaths, 2, scale)
rownames(subpaths_norm) <- rownames(subpaths)
superpaths <- lapply(seq_along(colnames(subpaths_norm)), function(i){
  superpath <- as.character(meta_metabo[2,which(meta_metabo[1,]==colnames(subpaths_norm)[i])])[1]
})
superpaths <- as.data.frame(unlist(superpaths))
rownames(superpaths) <- colnames(subpaths_norm)

# heatmap of subpaths, annotated by group, gender, and superpaths
pheatmap::pheatmap(subpaths_norm,
                   annotation_row = sample_subset[,c("GROUP", "GENDER")],
                   annotation_colors = list("GROUP"=c("non_tolerant"="#e31a1c90",
                                                      "primary_tolerant"="#00FF0090",
                                                      "secondary_tolerant"="#0000FF90")),
                   annotation_col = superpaths)

## random forest:
subpaths_annot <- cbind.data.frame(sample_subset[,c("GROUP", "GENDER")], subpaths, by = "row_names")
subpaths_annot <- subpaths_annot[,-c(2,80)]
colnames(subpaths_annot) <- make.names(colnames(subpaths_annot))
subpaths_annot$GROUP <- as.factor(subpaths_annot$GROUP)
rf_metabo <- randomForest::randomForest(GROUP~., subpaths_annot)
rf_metabo

subpaths_tol <- subpaths_annot[which(subpaths_annot$GROUP!="non_tolerant"),]
subpaths_tol$GROUP <- as.factor(as.character(subpaths_tol$GROUP))
rf_metabo <- randomForest::randomForest(GROUP~., subpaths_tol)
rf_metabo










