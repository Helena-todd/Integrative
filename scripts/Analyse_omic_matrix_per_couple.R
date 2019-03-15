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
gr_tmp <- as.character(pctgs_meta_couple$group)
gr_tmp[which(gr_tmp %in% c("primary_tolerant", "secondary_tolerant"))] <- "tolerant"
pctgs_meta_couple$group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(group~., pctgs_meta_couple, mtry = 3,
                                       importance = T, proximity = T)
rf_subst
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

metabo_couple <- df_tmp %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric, ~.[1]-.[2]) %>%
  mutate(group = unique(df_tmp[,c("GROUP", "COUPLENUMBER")])$GROUP) %>%
  column_to_rownames("COUPLENUMBER")

colnames(metabo_couple) <- make.names(colnames(metabo_couple))

rf_subst <- randomForest::randomForest(as.factor(group)~., metabo_couple)
rf_subst

# non tol vs tol:
gr_tmp <- as.character(metabo_couple$group)
gr_tmp[which(gr_tmp %in% c("primary_tolerant", "secondary_tolerant"))] <- "tolerant"
metabo_couple$group <- factor(gr_tmp)

rf_subst <- randomForest::randomForest(group~., metabo_couple, mtry = 50,
                                       importance = T, proximity = T)
rf_subst

######################################
######### on recipients only #########

load("outputs/data/metabo/recip/logdata.RData")
load("outputs/data/metabo/recip/meta_metabo.RData")
load("outputs/data/metabo/recip/recip_meta.RData")


