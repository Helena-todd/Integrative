#######################################################
##########    metabo donors and recipients   ##########
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

