### Analysis of the local cohort's RNAseq dataset

library(edgeR)
library("ggplot2")
library("triwise")
library("htmlwidgets")
library("openxlsx")
library('biomaRt')

## Liesbet's functions
###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)

  return(genes_de_sorted)
}

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))

  return(resultMatrix)
}

##### Load raw counts
countData<-read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/RNAseq/Global_table_StLouis_completed.xlsx")
colnames(countData) <- as.character(countData[1,]) # set patient names as colnames
rownames(countData) <- countData$ID # set gene names as rownames
countData_withinfo <- countData[,-1]

## Look at the info on samples:
plot(as.numeric(countData_withinfo[2,]), ylab = "unmapped genes", xlab = "patients")
plot(as.numeric(countData_withinfo[3,]), ylab = rownames(countData)[3], xlab = "patients")
above <- which(as.numeric(countData_withinfo[3,]) > 8*10^6)
colnames(countData_withinfo)[above]
plot(as.numeric(countData_withinfo[4,]), ylab = rownames(countData)[4], xlab = "patients")
plot(as.numeric(countData_withinfo[5,]), ylab = rownames(countData)[5], xlab = "patients")
colnames(countData_withinfo)[which(as.numeric(countData_withinfo[5,]) > 5*10^6)]
# [1] "D708"

countData <- countData[6:nrow(countData), 2:ncol(countData)] # rm the informations on unmapped, ambiguous...
dim(countData)
# [1] 60433    80

### Load meta data
colData<-read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/Local_cohort/Data synthesis local cohort Saint-Louis 17012019.xlsx")
colData <- colData[!is.na(colData$RNAseq.ID),] # rm rows with no RNAseq sample
rownames(colData) <- colData$Id.Cryostem.R
dim(colData)
# 80 55

### Reorder
countData<-countData[,rownames(colData)]
dim(countData)
# 60433    80

# ### Select samples and assign them a patient status
donor_id <- grep("D",colnames(countData))
recip_id <- grep("R",colnames(countData))

counts <- countData[,c(donor_id, recip_id)]

count_nb <- apply(counts,2,as.numeric) # turn characters into numerics
rownames(count_nb) <- rownames(counts) # restore the rownmes which were lost in the transfo
counts <- count_nb

### Preprocessing
y <- DGEList(counts = counts)
keep = rowSums(cpm(y)>1) >= 3
y = y[keep,]
dim(y)
## 18840    80

# Reset lib sizes
y$samples$lib.size = colSums(y$counts)

# Rename genes from ensembl ID -> gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(y$counts)
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=genes,
                mart= mart) # ! not all gene names were retrieved
not_found_at_all <- which(!rownames(y$counts)%in%G_list[[1]])
not_found <- which(G_list[[2]]=="")
not_change <- c(not_found, not_found_at_all)

ensembl_symbols <- G_list[[1]][-which(G_list[[2]]=="")]
gene_symbols <- G_list[[2]][-which(G_list[[2]]=="")]
ensembl_symbols <- ensembl_symbols[-which(duplicated(gene_symbols))]
gene_symbols <- gene_symbols[-which(duplicated(gene_symbols))]

counts_gene_symbol <- y$counts

for ( i in 1:length(gene_symbols)){
  rownames(counts_gene_symbol)[which(rownames(counts_gene_symbol)==ensembl_symbols[i])] <-
    gene_symbols[i]
}

y$counts <- counts_gene_symbol

# compute norm_factors
normfac <- calcNormFactors(y, method = "TMM")

#### Create design matrix
# for the three groups :
patient_groups <- tolower(colData$GROUP)
TS <- patient_groups[c(donor_id, recip_id)]

TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

#voom
#In limma, read counts are converted to log2-counts-per-million (logCPM) and
#the mean-variance relationship is modelled either with precision weights or
#with an empirical Bayes prior trend. The precision weights approach is called “voom”
# lib.size <- counts$samples$lib.size * counts$samples$norm.factors
# E <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
v <- voom(y, design, plot = TRUE)
expTable<-v$E

# check that we have removed size effects
vfit <- lmFit(v)
efit <- eBayes(vfit)
plotSA(efit)

DE_list <- list(v, design)

save(expTable, file = "outputs/data/RNAseq/expression_table.RData")
save(DE_list, file = "outputs/data/RNAseq/data_for_DE.RData")

################################################################################
########## GET DE GENES
################################################################################

#### Fit linear model on data
fit <- lmFit(v, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=non_tolerant-primary_tolerant,
                             group2=non_tolerant-secondary_tolerant,
                             group3=primary_tolerant-secondary_tolerant,
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)


#### Adjust P-values via Benjamini-Hochberg
########################################
##### 1. non-tol vs tol-1
########################################
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.05,2)
dim(DEgenesGroup1)
##35

########################################
##### 2. non-tol vs tol-2
########################################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,2)
dim(DEgenesGroup2)
##1

########################################
##### 3. tol-1 vs tol-2
########################################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.05,2)
dim(DEgenesGroup3)
##0

#######################################################################################################
############################## ------------- TRIWISE PLOT -------------- ##############################
#######################################################################################################

##### All DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2),rownames(DEgenesGroup3)))
length(allDEgenes)
##36


######################################################################
########## PREPARATION
######################################################################
#wantedColors<-c(nodiffall="gray55",diffall="indianred1")
wantedColors<-c(nodiffall="#FFFFFF80",diffall="indianred1")

##### Triwise plots
non_tol_mean<-apply(expTable[,as.numeric(which(design[,1]==1))],1,mean)
tol1_mean<-apply(expTable[,as.numeric(which(design[,2]==1))],1,mean)
tol2_mean<-apply(expTable[,as.numeric(which(design[,3]==1))],1,mean)

expTable_mean<-cbind(non_tol_mean, tol1_mean, tol2_mean)
colnames(expTable_mean)<-c('Non_tol','Primary_tol','Secondary_tol')

######################################################################
########## TRIWISE PLOTS
######################################################################
theBarycoords<-transformBarycentric(expTable_mean)
wantedColors<-c(nodiffall="#FFFFFF00",diffall="indianred1")

barycoords = theBarycoords
Gdiffexp = allDEgenes

###Triwise plot
p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors, showlabels = F)
print(p)
ggsave(p,file="results_allSamples/plots_samplesA/after_combat/triwisePlot.png",dpi=300)

#######################
### Plot only the DEgenes, to "clean up" the plots a bit:
expDE_mean <- expTable_mean[allDEgenes,]
DEBarycoords<-transformBarycentric(expDE_mean)
p<-plotDotplot(DEBarycoords, Gdiffexp=allDEgenes, colorvalues=wantedColors,
               showlabels = F, rmax = 10)
print(p)

## THIS FUNCTION DOESN'T TAKE THE rmax PARAMETER INTO ACCOUNT -> FIX IT FOR MARTIN
p<-interactiveDotplot(expDE_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE,
                      rmax = 5)#, Goi=IL17genes)
print(p)
saveWidget(p,file="~/Documents/VIB/Projects/Immgen/analysis/interactiveTriwisePlot_filtered_DEgenes_all_organs.html") ##needs full path!
#######################

###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)
ggsave(p,file="results_allSamples/plots_samplesA/after_combat/rosePlot.png",dpi=300)


save(expTable, allDEgenes, file = "~/Documents/VIB/Projects/Immgen/analysis/all_organs_expTable_allDEgenes.RData")
