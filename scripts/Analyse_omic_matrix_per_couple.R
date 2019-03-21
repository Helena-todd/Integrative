library(dplyr)
library(openxlsx)
library(metabolomics)
library(BioGVHD)
library(tidyverse)
library(data.table)
library(dendextend)


######## Trying different approaches to analyse the ########
########       Cytof, metabo and RNAseq data        ########

## RF ##
load("/Users/helenatodorov/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/cyto/3_backbones/backbone_2_D&Rall/pctgs_meta_rd_with_metaclust_labels.RData")

# on D-R:
df_tmp <- merge(data.frame(pctgs_meta_rd), samp_rd, by="row.names") %>%
  select(c(make.names(colnames(pctgs_meta_rd)), "Id.Cryostem.R", "GROUP", "COUPLENUMBER")) %>%
  arrange(COUPLENUMBER)

pctgs_meta_couple <- df_tmp %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric, ~.[1]-.[2]) %>%
  mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
  column_to_rownames("COUPLENUMBER")

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple)
rf_subst


# non tol vs tol:
rf <- rf_tol_non_tol(df_orig = pctgs_meta_couple, mtry = 3, ntree = 10000)
rf
# OOB estimate of  error rate: 26.47%
# Confusion matrix:
#   non_tolerant tolerant class.error
# non_tolerant           18        3   0.1428571
# tolerant                6        7   0.4615385
plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 2)
features_idx <- which(rf_subst$importance[,3]>=0.01)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(pctgs_meta_couple[,features_idx], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(pctgs_meta_couple$group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

pctgs_meta_couple <- pctgs_meta_couple %>%
  arrange(group)
pheatmap::pheatmap(pctgs_meta_couple[,features_idx],
                   cluster_rows = F)


get_info_rf_res <- function(df, rf_res, output_directory){
  round(rf_res$importance, 2)
  features_idx <- which(rf_res$importance[,3]>=0.01)
  features_of_interest <-
    rownames(rf_res$importance)[features_idx]

  subst.mds <- cmdscale(1 - rf_res$proximity, eig=TRUE)
  op <- par(pty="s")
  pairs(cbind(df[,features_idx], subst.mds$points), cex=0.6, gap=0,
        col=c("red", "blue")[as.numeric(df$group)],
        main="Predictors and MDS of Proximity Based on RandomForest")
  par(op)
  print(subst.mds$GOF)

  df_per_group <- df %>%
    arrange(group)
  pdf(paste0(output_directory, "pheatmap_main_features.pdf"))
  pheatmap::pheatmap(df_per_group[,features_idx],
                     cluster_rows = F)
  dev.off()

}


# tol 1 vs tol 2:
which(pctgs_meta_couple$group=="non_tolerant")
pctgs_meta_couple_tol <- pctgs_meta_couple[-which(pctgs_meta_couple$group=="non_tolerant"),]
pctgs_meta_couple_tol$group <- factor(as.character(pctgs_meta_couple_tol$group))

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple_tol, mtry = 3,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 53.85%
# Confusion matrix:
#   primary_tolerant secondary_tolerant class.error
# primary_tolerant                  4                  3   0.4285714
# secondary_tolerant                4                  2   0.6666667
plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 2)
features_idx <- which(rf_subst$importance[,3]>=0.01)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(pctgs_meta_couple[,features_idx], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(pctgs_meta_couple$group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

pctgs_meta_couple_tol <- pctgs_meta_couple_tol %>%
  arrange(group)

pheatmap::pheatmap(pctgs_meta_couple_tol[,-42],
                   cluster_rows = F)



df_tmp2 <- df_tmp %>% arrange(GROUP)
rownames(df_tmp2) <- df_tmp2$Id.Cryostem.R

gr <- data.frame(df_tmp2)
rownames(gr) <- rownames(df_tmp2)

pheatmap::pheatmap(df_tmp2[,-c(42:44)], cluster_rows = F, annotation_row = gr, fontsize_row = 6)




# on D/R:
df_tmp <- merge(data.frame(pctgs_meta_rd), samp_rd, by="row.names") %>%
  select(c(make.names(colnames(pctgs_meta_rd)), "Id.Cryostem.R", "GROUP", "COUPLENUMBER")) %>%
  arrange(COUPLENUMBER)

pctgs_meta_couple <- df_tmp %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric, ~.[1]/.[2]) %>%
  mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
  column_to_rownames("COUPLENUMBER")

# replace "inf" values by the maximum of each column + 1:
pctgs_meta_couple_tmp <- do.call(data.frame,lapply(pctgs_meta_couple[,-42],
                                               function(x) replace(x, is.infinite(x),(max(x[!is.infinite(x)])+1))))
pctgs_meta_couple_log <- pctgs_meta_couple_tmp %>%
  mutate_if(is.numeric, funs(log))

pctgs_meta_couple_tmp2 <- do.call(data.frame,lapply(pctgs_meta_couple_log,
                                                   function(x) replace(x, is.infinite(x),(min(x[!is.infinite(x)])-1))))

pctgs_meta_couple_tmp2 <- pctgs_meta_couple_tmp2 %>%
  mutate(group = pctgs_meta_couple$group)

pctgs_meta_couple <- pctgs_meta_couple_tmp2

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple)
rf_subst

# non tol vs tol:
gr_tmp <- as.character(pctgs_meta_couple$group)
gr_tmp[which(gr_tmp %in% c("primary_tolerant", "secondary_tolerant"))] <- "tolerant"
pctgs_meta_couple$group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple, mtry = 3,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 23.53%
# Confusion matrix:
#   non_tolerant tolerant class.error
# non_tolerant           17        4   0.1904762
# tolerant                4        9   0.3076923
plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 2)
features_idx <- which(rf_subst$importance[,3]>=0.01)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(pctgs_meta_couple[,features_idx], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(pctgs_meta_couple$group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

pctgs_meta_couple <- pctgs_meta_couple %>%
  arrange(group)
gr_annot <- as.data.frame(pctgs_meta_couple$group)
rownames(pctgs_meta_couple) <- paste0("couple_", rownames(pctgs_meta_couple))
rownames(gr_annot) <- rownames(pctgs_meta_couple)
pheatmap::pheatmap(pctgs_meta_couple[,features_idx],
                   cluster_rows = F, annotation_row = gr_annot)



# tol 1 vs tol 2:
which(pctgs_meta_couple$group=="non_tolerant")
pctgs_meta_couple_tol <- pctgs_meta_couple[-which(pctgs_meta_couple$group=="non_tolerant"),]
pctgs_meta_couple_tol$group <- factor(as.character(pctgs_meta_couple_tol$group))

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple_tol, mtry = 3,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 53.85%
# Confusion matrix:
#   primary_tolerant secondary_tolerant class.error
# primary_tolerant                  4                  3   0.4285714
# secondary_tolerant                4                  2   0.6666667
plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 2)
features_idx <- which(rf_subst$importance[,3]>=0.01)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(pctgs_meta_couple[,features_idx], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(pctgs_meta_couple$group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

pctgs_meta_couple_tol <- pctgs_meta_couple_tol %>%
  arrange(group)
rownames(pctgs_meta_couple_tol) <- paste0("couple_", rownames(pctgs_meta_couple_tol))

gr_annot <- as.data.frame(pctgs_meta_couple_tol$group)
rownames(gr_annot) <- rownames(pctgs_meta_couple_tol)

pheatmap::pheatmap(pctgs_meta_couple_tol[,-42],
                   cluster_rows = F, annotation_row = gr_annot)



#################################
####  on functional markers  ####
#################################

load("outputs/data/cyto/3_backbones/backbone_2_D&Rall/ratios_funct_meta_rd_old.RData")
load("outputs/data/cyto/3_backbones/backbone_2_D&Rall/funct_mark_res_orig.RData")

funct_big_table <- do.call(cbind, tables_res)
rownames(funct_big_table) <- gsub("^[0-9]*_([^_]*)_.*", "\\1", names(result))

df_tmp <- merge(data.frame(funct_big_table), samp_rd, by="row.names") %>%
  select(c(make.names(colnames(funct_big_table)), "Id.Cryostem.R", "GROUP", "COUPLENUMBER")) %>%
  arrange(COUPLENUMBER)

funct_couple <- df_tmp %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric, ~.[1]-.[2]) %>%
  mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
  column_to_rownames("COUPLENUMBER")

rf_subst <- randomForest::randomForest(group~., funct_couple)
rf_subst

gr_tmp <- as.character(funct_couple$group)
gr_tmp[which(gr_tmp %in% c("primary_tolerant", "secondary_tolerant"))] <- "tolerant"
funct_couple$group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(group~., funct_couple, mtry = 20,
                                       importance = T, proximity = T)
rf_subst






#####################################
##########   METABO data   ##########
#####################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/logdata.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/meta_metabo.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/big_mat.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/rd_meta.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/norm_data.RData")

df_tmp <- big_mat %>%
  select(c(colnames(norm_data)[-1], "Id.Cryostem.R", "GROUP", "COUPLENUMBER")) %>%
  arrange(COUPLENUMBER)

# #### on D-R:
#
# metabo_couple <- df_tmp %>%
#   group_by(COUPLENUMBER) %>%
#   summarise_if(is.numeric, ~.[1]-.[2]) %>%
#   mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
#   column_to_rownames("COUPLENUMBER")

#### on D/R:

metabo_couple <- df_tmp %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric, ~.[1]/.[2]) %>%
  mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
  column_to_rownames("COUPLENUMBER")

# replace "inf" values by the maximum of each column + 1:
metabo_couple_tmp <- do.call(data.frame,lapply(metabo_couple[,-541],
                                                   function(x) replace(x, is.infinite(x),(max(x[!is.infinite(x)])+1))))
metabo_couple_log <- metabo_couple_tmp %>%
  mutate_if(is.numeric, funs(log))

metabo_couple_tmp2 <- do.call(data.frame,lapply(metabo_couple_log,
                                                    function(x) replace(x, is.na(x),(min(x[!is.na(x)])-1))))

metabo_couple_tmp2 <- metabo_couple_tmp2 %>%
  mutate(group = metabo_couple$group)

metabo_couple <- metabo_couple_tmp2


##########################################################
#### This section can be used both on D-R andD on D/R ####
##########################################################

colnames(metabo_couple) <- make.names(colnames(metabo_couple))

rf_subst <- randomForest::randomForest(as.factor(group)~., metabo_couple, mtry = 30,
                                       ntree = 15000)
rf_subst
# OOB estimate of  error rate: 38.24%
# Confusion matrix:
#   Non_Tolerant Primary_tolerant Secondary_tolerant class.error
# Non_Tolerant                 21                0                  0           0
# Primary_tolerant              7                0                  0           1
# Secondary_tolerant            6                0                  0           1


# non tol vs tol:
gr_tmp <- as.character(metabo_couple$group)
gr_tmp[which(gr_tmp %in% c("Primary_tolerant", "Secondary_tolerant"))] <- "tolerant"
metabo_couple$group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(group~., metabo_couple, mtry = 30, ntree = 15000,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 38.24%
# Confusion matrix:
#   Non_Tolerant tolerant class.error
# Non_Tolerant           19        2   0.0952381
# tolerant               11        2   0.8461538

# on D-R:
# OOB estimate of  error rate: 20.59%
# Confusion matrix:
#   Non_Tolerant tolerant class.error
# Non_Tolerant           19        2   0.0952381
# tolerant                5        8   0.3846154

### on reduced_dimensions
pca <- prcomp(metabo_couple[,-541])
plot(pca$x, col = as.factor(metabo_couple$group))

rf_mat <- as.data.frame(pca$x[,1:20]) %>%
  mutate(group = metabo_couple$group)

rf_subst <- randomForest::randomForest(group~., rf_mat, mtry = 3, ntree = 15000,
                                       importance = T, proximity = T)
rf_subst # worth on reduced dimensions than on original dimensions
# OOB estimate of  error rate: 26.47%
# Confusion matrix:
#   Non_Tolerant tolerant class.error
# Non_Tolerant           20        1  0.04761905
# tolerant                8        5  0.61538462


# tol 1 vs tol 2:
metabo_couple_tol <- metabo_couple[-which(metabo_couple$group=="Non_Tolerant"),]
metabo_couple_tol$group <- factor(as.character(metabo_couple_tol$group))
colnames(metabo_couple_tol) <- make.names(colnames(metabo_couple_tol))

rf_subst <- randomForest::randomForest(group~., metabo_couple_tol, mtry = 30, ntree = 10000,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 69.23%
# Confusion matrix:
#   Primary_tolerant Secondary_tolerant class.error
# Primary_tolerant                  3                  4   0.5714286
# Secondary_tolerant                5                  1   0.8333333





######################################
######### on recipients only #########

load("outputs/data/metabo/recip/logdata.RData")
load("outputs/data/metabo/recip/meta_metabo.RData")
load("outputs/data/metabo/recip/big_mat.RData")
load("outputs/data/metabo/recip/norm_data.RData")

colnames(norm_data) <- make.names(colnames(norm_data))

rf_subst <- randomForest::randomForest(as.factor(Group)~., norm_data, mtry = 30,
                                       ntree = 15000)

rf_subst
# OOB estimate of  error rate: 20.59%
# Confusion matrix:
#   Non_Tolerant Primary_tolerant Secondary_tolerant class.error
# Non_Tolerant                 21                0                  0   0.0000000
# Primary_tolerant              2                4                  1   0.4285714
# Secondary_tolerant            3                1                  2   0.6666667



# non tol vs tol:
gr_tmp <- as.character(norm_data$Group)
gr_tmp[which(gr_tmp %in% c("Primary_tolerant", "Secondary_tolerant"))] <- "tolerant"
norm_data$Group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(Group~., norm_data, mtry = 30, ntree = 15000,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 8.82%
# Confusion matrix:
#   Non_Tolerant tolerant class.error
# Non_Tolerant           20        1  0.04761905
# tolerant                2       11  0.15384615

save(rf_subst, file = "outputs/data/metabo/RF/R_non_tol_tol/rf_R_non_tol_tol.Rdata")

plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 3)
features_idx <- which(rf_subst$importance[,3]>=0.005)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(norm_data[,names(features_idx)], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(norm_data$Group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

norm_data_2plot <- norm_data %>%
  rownames_to_column("Id_metabo") %>%
  arrange(Group) %>%
  column_to_rownames("Id_metabo")
colnames(norm_data_2plot) <- make.names(colnames(norm_data_2plot))

gr_annot <- as.data.frame(norm_data_2plot$Group)
rownames(gr_annot) <- rownames(norm_data_2plot)

pheatmap::pheatmap(norm_data_2plot[,names(features_idx)],
                   cluster_rows = F, annotation_row = gr_annot)


# tol1 vs tol2:
norm_data_tol <- norm_data[-which(norm_data$Group=="Non_Tolerant"),]
norm_data_tol$Group <- factor(as.character(norm_data_tol$Group))
colnames(norm_data_tol) <- make.names(colnames(norm_data_tol))

rf_subst <- randomForest::randomForest(Group~., norm_data_tol, mtry = 40, ntree = 15000,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 46.15%
# Confusion matrix:
#   Primary_tolerant Secondary_tolerant class.error
# Primary_tolerant                  5                  2   0.2857143
# Secondary_tolerant                4                  2   0.6666667



# non tol vs tol on reduced_dim:
gr_tmp <- as.character(norm_data$Group)
gr_tmp[which(gr_tmp %in% c("Primary_tolerant", "Secondary_tolerant"))] <- "tolerant"
norm_data$Group <- factor(gr_tmp)
pca <- prcomp(norm_data[,-1])
plot(pca$x, col = norm_data$Group)

rf_mat <- as.data.frame(pca$x[,1:20]) %>%
  mutate(group = norm_data$Group)

rf_subst <- randomForest::randomForest(group~., rf_mat, mtry = 3, ntree = 15000,
                                       importance = T, proximity = T)
rf_subst
# OOB estimate of  error rate: 8.82%
# Confusion matrix:
#   Non_Tolerant tolerant class.error
# Non_Tolerant           20        1  0.04761905
# tolerant                2       11  0.15384615

save(rf_subst, file = "outputs/data/metabo/RF/R_non_tol_tol/rf_R_non_tol_tol.Rdata")

plot(rf_subst)
tree_func(final_model = rf_subst)

round(rf_subst$importance, 3)
features_idx <- which(rf_subst$importance[,3]>=0.005)
features_of_interest <-
  rownames(rf_subst$importance)[features_idx]

subst.mds <- cmdscale(1 - rf_subst$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(norm_data[,names(features_idx)], subst.mds$points), cex=0.6, gap=0,
      col=c("red", "blue")[as.numeric(norm_data$Group)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(subst.mds$GOF)

norm_data_2plot <- norm_data %>%
  rownames_to_column("Id_metabo") %>%
  arrange(Group) %>%
  column_to_rownames("Id_metabo")
colnames(norm_data_2plot) <- make.names(colnames(norm_data_2plot))

gr_annot <- as.data.frame(norm_data_2plot$Group)
rownames(gr_annot) <- rownames(norm_data_2plot)

pheatmap::pheatmap(norm_data_2plot[,names(features_idx)],
                   cluster_rows = F, annotation_row = gr_annot)







#### Distance pca D-R on metabolites:

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/logdata.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/meta_metabo.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/big_mat.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/rd_meta.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/norm_data.RData")

### Plot pctgs' PCA:
pca <- prcomp(norm_data[,-1])
plot(pca$x, col = norm_data$Group)

distances = c()

for (i in names(table(big_mat$COUPLENUMBER))){
  print(i)
  group_status <- big_mat$GROUP[which(big_mat$COUPLENUMBER==i)]
  pt_coord <- pca$x[rownames(big_mat)[which(big_mat$COUPLENUMBER==i)],c(1:2)]
  if(group_status[1]=="Non_Tolerant"){
    lines(pt_coord,col="red")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  } else if(group_status[1]=="Primary_tolerant"){
    lines(pt_coord,col="green")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  } else {
    lines(pt_coord,col="blue")
    points(pt_coord[1,1], pt_coord[1,2], pch = 19)
  }
  distances = c(distances, sqrt((pt_coord[2,1]-pt_coord[1,1])^2 + (pt_coord[2,2]-pt_coord[1,2])^2))
}

list_couples <- names(table(big_mat$COUPLENUMBER))
names(distances) <- paste0("couple_",list_couples)

df_tmp <- big_mat[,c("COUPLENUMBER", "GROUP")] %>%
  unique() %>%
  arrange(COUPLENUMBER)
dis_colors <- c("red","green","blue")[as.factor(df_tmp$GROUP)]
order_dis <- order(distances)

plot(distances[order_dis], col = dis_colors[order_dis], pch=19,
     main = "Distance between donor and recipient from same couple",
     ylab = "Distance D/R", xaxt = "n", xlab = "")
axis(1, at=1:34, labels=FALSE)
text(x=1:34, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels=names(distances[order_dis]), srt=45, adj=1, xpd=TRUE, cex = .7)
legend("topleft", c("non_tolerant","primary_tolerant","secondary_tolerant"),
       col = c("red","green","blue"), pch = 19)

legend("topleft", names(table(samp_rd$AGE)),
       col = c("blue","green","yellow","orange","red","pink"), pch=19)



