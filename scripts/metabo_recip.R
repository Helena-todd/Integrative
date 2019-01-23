library(dplyr)
library(openxlsx)
library(metabolomics)
library(BioGVHD)
library(tidyverse)
library(data.table)


#################################
########  PREPROCESSING  ########
#################################

samp_rd <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_RNAseq_14:01:19/Data synthesis local cohort Saint-Louis 17012019.xlsx")

info <- extract_info_metabo(metadata_file = samp_rd,
                                metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis_filtered.xlsx")
rd_meta <- info$subset_meta
meta_metabo <- info$meta_metabo
data_metabolites <- info$data_metabolites

##### filtering #####

data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics

colsums <- by(data_metabolites, rd_meta$GROUP, # identify metabolites <50% in each group
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
mat2use <- merge.data.frame(as.data.frame(rd_meta[,2], row.names = rownames(rd_meta)),
                            data_metabo, by = "row.names") %>%
  tibble::column_to_rownames(var="Row.names")

logdata <- LogTransform(mat2use)
normdata <- Normalise(logdata$output, method = "median")
colnames(normdata$output) <- colnames(logdata$output)
norm_data <- normdata$output

# merge dataframes:
big_mat <- merge.data.frame(normdata$output, rd_meta, by = "row.names") %>%
  tibble::column_to_rownames("Row.names")

logdata <- logdata$output
save(logdata, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/logdata.RData")
save(norm_data, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/norm_data.RData")
meta_metabo <- meta_metabo[,which(as.character(meta_metabo[3,])%in%colnames(norm_data))] # only selected metabolites
save(meta_metabo, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/meta_metabo.RData")
save(big_mat, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/big_mat.RData")
save(rd_meta, file="~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/rd_meta.RData")



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

res <- find_subpatways_driving_PC(meta_metabo, pca = pca, PC = 2)
ordered_subpathways_pc2 <- res$ordered_subpathways
princ_axes <- res$princ_axes

save(ordered_subpathways_pc2, file="ordered_subpathways_pc2.RData")
ordered_metabolites_pc1 <- princ_axes[order(princ_axes$CS1),c(1,2,5)]
save(ordered_metabolites_pc1, file = "ordered_metabolites_pc1.RData")


###############################
#### only on tol1 and tol2 ####
###############################

mat2use <- norm_data[which(big_mat$GROUP!="non_tolerant"),]

pca <- dudi.pca(mat2use[,-1], center = T, scale = T)
4
s.label(pca$li, clabel = .7)
s.class(pca$li, as.factor(mat2use$Group), col=c("blue","green"))
s.class(pca$li[,2:3], as.factor(mat2use$Group), col=c("blue","green"))

# which metabolites/ sub-pathways/ super-pathways drive the 3rd pca axis (difference between tol1 and tol2 ++)?
res <- find_subpatways_driving_PC(meta_metabo, pca = pca, PC = 3)
ordered_subpathways_pc3 <- res$ordered_subpathways
head(ordered_subpathways_pc3)
princ_axes <- res$princ_axes
ordered_metabolites_pc3 <- princ_axes[order(princ_axes$CS3),c(1,2,7)]
save(ordered_metabolites_pc3, file = "ordered_metabolites_pc3.RData")
save(ordered_subpathways_pc3, file="ordered_subpathways_pc3.RData")
write.xlsx(ordered_metabolites_pc3, file = "ordered_metabolites_along_PC3.xls")

# which metabolites/ sub-pathways drive the 1st pca axis (variability ++)?
res <- find_subpatways_driving_PC(meta_metabo, pca = pca, PC = 1)
ordered_subpathways_pc1 <- res$ordered_subpathways
ordered_metabolites_pc1 <- res$princ_axes[order(res$princ_axes$CS1),c(1,2,5)]
save(ordered_metabolites_pc1, file = "ordered_metabolites_tol_pc1.RData")
save(ordered_subpathways_pc1, file="ordered_subpathways_tol_pc1.RData")

pca <- dudi.pca(big_mat[which(big_mat$GROUP!="non_tolerant"),2:590], center = T, scale = T)
4
s.class(pca$li, as.factor(big_mat$GROUP[which(big_mat$GROUP!="non_tolerant")]), col=c("blue","green"))
s.class(pca$li, as.factor(big_mat$GENDER[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))
s.class(pca$li, as.factor(big_mat$cGVHD[which(big_mat$GROUP!="non_tolerant")]), col=c("#bfd3e6","#88419d"))
s.class(pca$li, as.factor(big_mat$SOURCEOFGRAFT[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))
s.class(pca$li, as.factor(big_mat$DONORCMV[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))
s.class(pca$li, as.factor(big_mat$DONORSEX[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))
s.class(pca$li, as.factor(big_mat$GROUPE[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))
s.class(pca$li, as.factor(big_mat$CMVStatus[which(big_mat$GROUP!="non_tolerant")]), col=c("red","blue"))


recip_info <- merge.data.frame(pca$li, recip_meta[which(recip_meta$GROUP!="non_tolerant"),],
                               by = "row.names") %>%
  arrange(Axis1)
bla <- lm(Axis1~GROUP*GENDER*CMVStatus*GROUPE*DONORSEX*DONORCMV*DONORGROUPE, recip_info)


### differential analysis

res <- TwoGroup(mat2use)
VolcanoPlot(res$output[,4], res$output[,2], cexlab = 0.6)
MetBoxPlots(mat2use, "N-methylproline",cols = c("blue","green"),main = "N-methylproline")
MetBoxPlots(mat2use, "chiro-inositol",cols = c("blue","green"),main = "chiro-inositol")
MetBoxPlots(mat2use, "1,7-dimethylurate",cols = c("blue","green"),main = "1,7-dimethylurate")
MetBoxPlots(mat2use, "1-methylurate",cols = c("blue","green"),main = "1-methylurate")

TwoGroupPlots(mat2use[,-1], res$output[,1], foldchanges = res$output[,4], pvalues = res$output[,2],
              padjmethod = "BH", fcutoff = log(2), pcutoff = 0.05)






####################################################
####### Turn metabolites -> subpathway table #######
####################################################

load("outputs/data/metabo/recip/logdata.RData")
load("outputs/data/metabo/recip/meta_metabo.RData")
load("outputs/data/metabo/recip/recip_meta.RData")

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
save(subpaths, file = "subpaths_table.RData")

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
                   annotation_row = recip_meta[,c("GROUP", "GENDER")],
                   annotation_colors = list("GROUP"=c("non_tolerant"="#e31a1c90",
                                                      "primary_tolerant"="#00FF0090",
                                                      "secondary_tolerant"="#0000FF90")),
                   annotation_col = superpaths)

## random forest:
subpaths_annot <- cbind.data.frame(recip_meta[,c("GROUP", "GENDER")], subpaths, by = "row_names")
subpaths_annot <- subpaths_annot[,-c(2,80)]
colnames(subpaths_annot) <- make.names(colnames(subpaths_annot))
subpaths_annot$GROUP <- as.factor(subpaths_annot$GROUP)
rf_metabo <- randomForest::randomForest(GROUP~., subpaths_annot)
rf_metabo

subpaths_tol <- subpaths_annot[which(subpaths_annot$GROUP!="non_tolerant"),]
subpaths_tol$GROUP <- as.factor(as.character(subpaths_tol$GROUP))
rf_metabo <- randomForest::randomForest(GROUP~., subpaths_tol)
rf_metabo















###############################################################################
###############################################################################
##################     ANALYSIS ON DONORS AND RECIPIENTS     ##################
###############################################################################
###############################################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/logdata.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/meta_metabo.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/big_mat.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/rd_meta.RData")

### Look at the subpathways instead of metabolites
data_metabo <- as.data.frame(logdata[,-1])
meta_metabo <- meta_metabo[,which(as.character(meta_metabo[3,])%in%colnames(data_metabo))]
data_metabo <- data_metabo[,which(colnames(data_metabo)%in%as.character(meta_metabo[3,]))]

t_data_metabo <- t(data_metabo)
subpath_info <- as.character(meta_metabo[1,])
data_meta <- as.data.table(t_data_metabo)
binded <- cbind(data_meta, subpath_info)

subpaths <- binded[,lapply(.SD, mean), by=subpath_info]
subpaths <- as.data.frame(subpaths) %>% column_to_rownames("subpath_info")
subpaths <- as.data.frame(t(subpaths))

# I add the group and couple nb info in 2 extra columns:
sub_data <- cbind(logdata[,1], subpaths)
colnames(sub_data)[1] <- "group"
sub_data <- cbind(big_mat$COUPLENUMBER, sub_data)
colnames(sub_data)[1] <- "couplenb"

save(sub_data, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/sub_data.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/sub_data.RData")

### paired t-tests:
sub_D <- sub_data[grep(pattern = "D", rownames(sub_data)),]
sub_D <- sub_D[order(sub_D$couplenb),]

sub_R <- sub_data[grep(pattern = "R", rownames(sub_data)),]
sub_R <- sub_R[order(sub_R$couplenb),]

ttests <- function(mat_D, mat_R, group){
  mat_D_gr <- mat_D[which(mat_D$group==group),]
  mat_R_gr <- mat_R[which(mat_R$group==group),]

  ttests_gr <- sapply(3:ncol(mat_D), function(subp){
    paired_ttest <- t.test(as.numeric(mat_D_gr[, subp]),
                           as.numeric(mat_R_gr[, subp]), paired = T)
    pval <- paired_ttest$p.value
    names(pval) <- colnames(sub_D)[subp]
    pval
  })
  return(ttests_gr)
}

ttests_NT <- ttests(sub_D, sub_R, "Non_Tolerant")
# correct for multiple tests
ttests_NT <- p.adjust(ttests_NT, "BH")

DE_NT <- ttests_NT[which(ttests_NT <0.05)]
table_NT <- data.frame(as.numeric(sort(DE_NT)))
rownames(table_NT) <- names(sort(DE_NT))
colnames(table_NT) <- "pvalue"
table_NT
length(DE_NT)
# 33

ttests_1T <- ttests(sub_D, sub_R, "Primary_tolerant")
# correct for multiple tests
ttests_1T <- p.adjust(ttests_1T, "BH")

DE_1T <- ttests_1T[which(ttests_1T <0.05)]
table_1T <- data.frame(as.numeric(sort(DE_1T)))
rownames(table_1T) <- names(sort(DE_1T))
colnames(table_1T) <- "pvalue"
table_1T
length(DE_1T)
# 0

ttests_2T <- ttests(sub_D, sub_R, "Secondary_tolerant")
# correct for multiple tests
ttests_2T <- p.adjust(ttests_2T, "BH")

DE_2T <- ttests_2T[which(ttests_2T <0.05)]
table_2T <- data.frame(as.numeric(sort(DE_2T)))
rownames(table_2T) <- names(sort(DE_2T))
colnames(table_2T) <- "pvalue"
table_2T
length(DE_2T)
# 0





### Normalise? :
subpaths_norm <- apply(subpaths, 2, scale)
rownames(subpaths_norm) <- rownames(subpaths)
