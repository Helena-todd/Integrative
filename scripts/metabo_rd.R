library(dplyr)
library(openxlsx)
library(metabolomics)
library(BioGVHD)
library(tidyverse)
library(data.table)
library(dendextend)


#################################
########  PREPROCESSING  ########
#################################

samp_rd <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/Data synthesis local cohort Saint-Louis 17012019.xlsx")

info <- extract_info_metabo(metadata_file = samp_rd,
                            metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/Metabo/Metabolomic local cohort Saint-Louis_filtered.xlsx")
rd_meta <- info$subset_meta
meta_metabo <- info$meta_metabo
data_metabolites <- info$data_metabolites

##### filtering #####

data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics
data_metabo_national <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/National_cohort/Metabo/Metabo NATIONAL cohort CRYOSTEM.xlsx")
national_metabolites <- data_metabo_national$X1[-c(1,2)]
# remove metabolites that are not present in both the national and the St_Louis cohort :
data_metabolites <- data_metabolites[,which(colnames(data_metabolites) %in% national_metabolites)]

##### generate figure similar to David's #####
names_patients <- rownames(data_metabolites)

big_mat <- data_metabolites %>%
  mutate_all(as.numeric) %>%
  mutate(METABONAME = names_patients) %>%
  left_join(rd_meta[,c("GROUP", "COUPLENUMBER", "METABONAME")], by = "METABONAME") %>%
  mutate(status = ifelse(test = grepl("D", METABONAME), yes = "D", no = "R"))

na_pctgs <- big_mat %>%
  group_by(GROUP, status) %>%
  summarize_all(funs(sum(is.na(.)) / length(.)))

melted <- melt(na_pctgs, id = 1:2)

## marche pas:
# to_rm <- melted %>%
#   group_by(variable) %>%
#   summarize_if(is.integer, funs(all(. > 0.5)))

to_rm <- lapply(names(table(melted$variable)), function(metabolite){
  values <- melted$value[which(melted$variable == metabolite)]
  all(values > 0.5)
})

names2rm <- as.character(names(table(melted$variable))[which(to_rm == T)])
plot2rm <- melted[which(melted$variable %in% names2rm),]
table_pctg_na <- plot2rm %>% arrange(status) %>%
  mutate(variable = paste0(variable, GROUP))

library(ggplot2)
ggplot(data = table_pctg_na,
       mapping = aes(x = variable, fill = GROUP,
                     y = ifelse(test = status == "D",
                                yes = -value, no = value))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(table_pctg_na$value) * c(-1,1)) +
  labs(y = "pctg of NA") +
  coord_flip()

##### new filtering : > 50% na in all groups, donors and recipients

high_na <- which(to_rm==T)
data_metabo <- data_metabolites[ , -high_na]

##### previous filtering : >50% in at least one group

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



#####################################################
##########  Generate data recipients only  ##########
#####################################################

samp_rd <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/Data synthesis local cohort Saint-Louis 17012019.xlsx")

info <- extract_info_metabo(metadata_file = samp_rd,
                            metabo_file = "~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/Metabo/Metabolomic local cohort Saint-Louis_filtered.xlsx")
rd_meta <- info$subset_meta
meta_metabo <- info$meta_metabo
data_metabolites <- info$data_metabolites
names_patients <- rownames(data_metabolites)[grep("R", rownames(data_metabolites))]
data_metabo_national <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/National_cohort/Metabo/Metabo NATIONAL cohort CRYOSTEM.xlsx")


##### filtering metabolites #####

res_preprocessing <- metabo_preprocess(patient_type = "recipients", data_metabolites, meta_metabo,
                                       data_metabolites_national, rd_meta, names_patients,
                                       pdf_name = "outputs/plots/metabo/recip/preprocess_removed_metabolites.pdf",
                                       pdf_variance_name = "outputs/plots/metabo/recip/preprocess_variance_cutoff.pdf",
                                       save_results = F, save_res_repository = NULL)


data_metabolites <- data_metabolites[,-which((meta_metabo[2,]=="Xenobiotics")&(meta_metabo[1,]=="Drug"))] # rm drug xenobiotics
data_metabo_national <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/National_cohort/Metabo/Metabo NATIONAL cohort CRYOSTEM.xlsx")
national_metabolites <- data_metabo_national$X1[-c(1,2)]
# remove metabolites that are not present in both the national and the St_Louis cohort :
data_metabolites <- data_metabolites[,which(colnames(data_metabolites) %in% national_metabolites)]

##### generate figure similar to David's #####
names_patients <- rownames(data_metabolites)[grep("R", rownames(data_metabolites))]

big_mat <- data_metabolites[names_patients,] %>%
  mutate_all(as.numeric) %>%
  mutate(METABONAME = names_patients) %>%
  left_join(rd_meta[,c("GROUP", "COUPLENUMBER", "METABONAME")], by = "METABONAME")

na_pctgs <- big_mat %>%
  group_by(GROUP) %>%
  summarize_all(funs(sum(is.na(.)) / length(.)))

melted <- melt(na_pctgs, id = 1:2)

## marche pas:
# to_rm <- melted %>%
#   group_by(variable) %>%
#   summarize_if(is.integer, funs(all(. > 0.5)))

to_rm <- lapply(names(table(melted$variable)), function(metabolite){
  values <- melted$value[which(melted$variable == metabolite)]
  all(values > 0.5)
})

names2rm <- as.character(names(table(melted$variable))[which(to_rm == T)])
plot2rm <- melted[which(melted$variable %in% names2rm),]
table_pctg_na <- plot2rm %>%
  mutate(variable = paste0(variable, GROUP))

library(ggplot2)
ggplot(data = table_pctg_na,
       mapping = aes(x = variable, fill = GROUP,
                     y = value)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = abs, limits = max(table_pctg_na$value) * c(0,1)) +
  labs(y = "pctg of NA") +
  coord_flip()

##### new filtering : > 50% na in all groups, donors and recipients

high_na <- which(to_rm==T)
data_metabo <- data_metabolites[names_patients , -high_na]

##### previous filtering : >50% in at least one group

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
abline(h=quantile(sds, 0.3), col="red")
data_metabo <- data_metabo[,which(sds>=quantile(sds, 0.3))]
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





#################################################
##########  ANALYSIS  recipients only  ##########
#################################################

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
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/norm_data.RData")

### Look at clusters of similar metabolites instead of metabolites directly
data_metabo <- as.data.frame(logdata[,-1])
meta_metabo <- meta_metabo[,which(as.character(meta_metabo[3,])%in%colnames(data_metabo))]

subst_mat <- big_mat %>%
  group_by(COUPLENUMBER) %>%
  summarise_if(is.numeric,~.[1]-.[2])

data_metabo <- subst_mat[,-c(1,542:557)]
rownames(data_metabo) <- subst_mat$COUPLENUMBER
info_gr_couple <- unique(big_mat[,c("Group", "COUPLENUMBER")]) %>%
  rownames_to_column("METABONAME") %>%
  column_to_rownames("COUPLENUMBER")
gr_info <- info_gr_couple[rownames(data_metabo),"Group"]


# kmeans clustering to find groups of metabolites
km100 <- kmeans(t(data_metabo), centers = 100, iter.max = 10000)
km100_info <- rbind(meta_metabo, km100$cluster)
table(as.integer(km100_info[4,]))
clusters_km100 <- lapply(seq_len(100), function(clst){
  km100_info[c(3,1,2),which(km100_info[4,]==clst)]
})
km100_info[c(3,1,2),which(km100_info[4,]==(38))]

write.xlsx(clusters_km100,
           file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/resulting_clusters_km100.xlsx")

## visualise metabolite clusters in patients:
mat_clust <- as.data.frame(t(data_metabo)) %>%
  mutate(clustering = km100$cluster) %>%
  group_by(clustering) %>%
  summarize_all(mean) %>%
  select(-clustering)

matintr <- as.data.frame(t(mat_clust)[rownames(data_metabo),])

mat2heatmap <- matintr %>%
  rownames_to_column("patient") %>%
  mutate(group = gr_info) %>%
  arrange(group) %>%
  column_to_rownames("patient")

annot_row <- as.data.frame(mat2heatmap$group)
rownames(annot_row) <- rownames(mat2heatmap)

pheatmap::pheatmap(mat2heatmap[,-101], cluster_rows = F,
                   annotation_row = annot_row)










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



#######################################################
##########    extract metabolites of inter   ##########
##########           using RF models         ##########
#######################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/norm_data.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/meta_metabo.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/big_mat.RData")
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/rd_meta.RData")

head(big_mat[,1:5])
vars2rm <- colnames(rd_meta)[-c(4, 10)]
grouped_mat <- big_mat %>%
  group_by(COUPLENUMBER) %>%
  select(-vars2rm)

subst <- grouped_mat %>%
  summarise_if(is.numeric, ~.[1]-.[2])# %>%
  #rename_all(funs(paste0("subst_")))
colnames(subst)[-1] <- paste0("subst_", colnames(subst))[-1]

frac <- grouped_mat %>%
  summarise_if(is.numeric, ~log(.[1]/.[2]))
colnames(frac)[-1] <- paste0("frac_", colnames(frac))[-1]

patient_info <- big_mat %>%
  mutate(status = ifelse(test = grepl("D", METABONAME), yes = "D", no = "R"))

donors <- patient_info %>%
  filter(status == "D")
colnames(donors)[-1] <- paste0("donor_", colnames(donors))[-1]

recip <- patient_info %>%
  filter(status == "R") %>%
  select(-c(608:610, 612:663))
colnames(recip)[-c(1,608)] <- paste0("recip_", colnames(recip))[-c(1, 608)]

save(donors, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo_national/donors.RData")
save(recip, file = "~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo_national/recip.RData")

couple_mat <- plyr::join_all(list(subst, frac, donors[,-1], recip), by = "COUPLENUMBER")

head(couple_mat[,1:5])
rownames(couple_mat) <- couple_mat$COUPLENUMBER
couple_mat <- couple_mat[,-1]
new_mat <- couple_mat
colnames(new_mat) <- make.names(colnames(couple_mat))


### RF
rf_r<-randomForest::randomForest(Group~., new_mat, ntree=15000, mtry=100)
rf_r
plot(rf_r)
tree_func(final_model = rf_r)

### RF non tol vs tol
gr <- as.character(new_mat$Group)
gr[which((gr == "Primary_tolerant")|(gr == "Secondary_tolerant"))] <- "Tolerant"
tol_non_tol_mat <- new_mat %>%
  mutate(Group = as.factor(gr))
rf_r<-randomForest::randomForest(Group~., tol_non_tol_mat, ntree=15000,
                                 mtry=50, importance = T)
rf_r
plot(rf_r)
tree_func(final_model = rf_r)
rownames(rf_r$importance)[which(rf_r$importance==max(rf_r$importance))]

topca <- tol_non_tol_mat %>% select_if(is.numeric)
pca <- prcomp(topca)
plot(pca$x,  col = tol_non_tol_mat$Group)
sne <- Rtsne::Rtsne(topca, perplexity = 8)
plot(sne$Y, col = tol_non_tol_mat$Group)



#######################################################
##########    metabo donors and recipients   ##########
##########        compatibility_scores       ##########
#######################################################

load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/donors/subpaths_table_donors.RData")
subpaths_donors <- subpaths
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/recip/subpaths_table.RData")
subpaths_recip <- subpaths

metabo_samples <- sample_subset <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018.xlsx",
                                             check.names = FALSE) %>%
  mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
         GROUP = tolower(GROUP)) %>%
  dplyr::filter(!is.na(METABONAME)) %>%
  dplyr::select_if(~sum(!is.na(.)) > 0)

rd_couples <- metabo_samples$COUPLENUMBER[which(duplicated(metabo_samples$COUPLENUMBER)==TRUE)]
samp_rd <- metabo_samples[which(metabo_samples$COUPLENUMBER %in% rd_couples),]
rownames(samp_rd) <- samp_rd$METABONAME

compatibility_scores <- compute_compatibility_score(samp_rd)

## compute pca :
sub_donors <- subpaths_donors[which(rownames(subpaths_donors)%in%rownames(samp_rd)),]
sub_donors <- sub_donors[,which(colnames(sub_donors)%in%colnames(sub_recip))]
sub_recip <- subpaths_recip[which(rownames(subpaths_recip)%in%rownames(samp_rd)),]
sub_recip <- sub_recip[,which(colnames(sub_recip)%in%colnames(sub_donors))]

mat2use <- rbind.data.frame(sub_donors, sub_recip, by = "col_names")
mat2use <- mat2use[rownames(samp_rd),]
mat2use_num <- apply(mat2use, 2, as.numeric)
rownames(mat2use_num) <- rownames(mat2use)
pca_metabo_rd <- prcomp(mat2use_num)


pca_dist_scores <- compute_pca_distance_scores(pca_metabo_rd, samp_rd)

couple_group <- lapply(unique(samp_rd$COUPLENUMBER), function(x){
  couple <- samp_rd[which(samp_rd$COUPLENUMBER==x),]
  group <- couple$GROUP[1]
})
couple_group <- unlist(couple_group)

table1 <- cbind( compatibility_scores, pca_dist_scores)
table1 <- table1[,-6]
table1 <- as.data.frame(table1) %>%
  mutate(group = couple_group)
save(table1, file = "table1_metabo.RData")

table1 %>% group_by(group) %>% summarise(mean(pc1_diff), sd(pc1_diff))
table1 %>% group_by(group) %>% summarise(mean(pc2_diff), sd(pc2_diff))


# random forest
load("~/Documents/VIB/Projects/Integrative_Paris/Integrative/outputs/data/metabo/r&d/table1_metabo.RData")

table1 <- table1[,-1]
table1$group <- as.factor(table1$group)
rf_r<-randomForest::randomForest(group~., table1, ntree=15000, mtry=3)
rf_r
plot(rf_r)
tree_func(final_model = rf_r)

red_table <- table1[which(table1$group!="non_tolerant"),]
red_table$group <- as.factor(as.character(red_table$group))
rf_r<-randomForest::randomForest(group~., red_table, ntree=15000, mtry=2)
rf_r
plot(rf_r)
tree_func(final_model = rf_r)


