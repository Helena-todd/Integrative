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
countData<-read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_RNAseq_14:01:19/RNAseq Count Saint_Louis (n=80)/Global_Table_Count _Saint_Louis_BASE_01.2019.xlsx")
colnames(countData) <- as.character(countData[1,]) # set patient names as colnames
rownames(countData) <- countData$ID # set gene names as rownames
countData_withinfo <- countData

## Look at the info on samples:
plot(as.numeric(countData_withinfo[2,-1]), ylab = "unmapped genes", xlab = "patients")
plot(as.numeric(countData_withinfo[3,-1]), ylab = rownames(countData)[3], xlab = "patients")
above <- which(as.numeric(countData_withinfo[3,-1]) > 8*10^6)
plot(as.numeric(countData_withinfo[4,-1]), ylab = rownames(countData)[4], xlab = "patients")
plot(as.numeric(countData_withinfo[5,-1]), ylab = rownames(countData)[5], xlab = "patients")
colnames(countData)[which(as.numeric(countData_withinfo[5,-1]) > 5*10^6)]
# [1] "D708"

countData <- countData[6:nrow(countData), 2:ncol(countData)] # rm the informations on unmapped, ambiguous...
dim(countData)
# [1] 60433    80

### Load meta data
colData<-read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_RNAseq_14:01:19/Data synthesis local cohort Saint-Louis 15012019.xlsx")
colData <- colData[!is.na(colData$RNAseq.ID),] # rm rows with no RNAseq sample
rownames(colData) <- colData$Id.Cryostem.R
dim(colData)
# 80 53

### Reorder
countData<-countData[,rownames(colData)]
dim(countData)
# 60433    80

# ### Select samples and assign them a patient status
donor_id <- grep("D",colnames(countData))
recip_id <- grep("R",colnames(countData))

counts <- countData[,c(donor_id, recip_id)]
counts <- counts[,-18] # only NAs after 7143th row
colData <- colData[-which(!rownames(colData)%in%colnames(counts)),] #also remove it in colData

count_nb <- apply(counts,2,as.numeric) # turn characters into numerics
rownames(count_nb) <- rownames(counts) # restore the rownmes which were lost in the transfo
counts <- count_nb

### Preprocessing
y <- DGEList(counts = counts)
keep = rowSums(cpm(y)>1) >= 3
y = y[keep,]
dim(y)
## 18827    79

# Reset lib sizes
y$samples$lib.size = colSums(y$counts)

# Rename genes grom ensembl ID -> gene symbol
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

counts_gene_symbol <- y$counts[which(rownames(y$counts)%in%ensembl_symbols),]
counts_gene_symbol <- as.data.frame(counts_gene_symbol)

which(duplicated(gene_symbols))



which(rownames(counts_gene_symbol)%in%ensembl_symbols)
row.names(counts_gene_symbol[ensembl_symbols,]) <- gene_symbols

rownames(y$counts[ensembl_symbols,]) <- gene_symbols


# compute norm_factors
normfac <- calcNormFactors(y, method = "TMM")

#### Create design matrix
# for the three groups :
patient_groups <- as.character(colData$GROUP)
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
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.05,1)
dim(DEgenesGroup1)
##3

########################################
##### 2. non-tol vs tol-2
########################################
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=2)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,1)
dim(DEgenesGroup2)
##0

########################################
##### 3. tol-1 vs tol-2
########################################
allGenesGroup3<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=3)
DEgenesGroup3<-getDEgenes(allGenesGroup3,0.05,1)
dim(DEgenesGroup3)
##0

#######################################################################################################
############################## ------------- TRIWISE PLOT -------------- ##############################
#######################################################################################################

##### All DE genes #####
allDEgenes<-unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2),rownames(DEgenesGroup3)))
length(allDEgenes)
##5738


######################################################################
########## PREPARATION
######################################################################
#wantedColors<-c(nodiffall="gray55",diffall="indianred1")
wantedColors<-c(nodiffall="#FFFFFF80",diffall="indianred1")

##### Triwise plots
# colsSPA<-grep("A1_|A3_|A4_|A5_|A6_|A7_",colnames(expTable))
# colsRA<-grep("A9_|A11_|A12_",colnames(expTable))
# colsHealthy<-grep("A13_|A14_|A15_|A16_",colnames(expTable))

# for liver:
# cDC1_mean<-apply(expTable[,1:3],1,mean)
# cDC2_mean<-apply(expTable[,4:7],1,mean)
# mac_mean<-apply(expTable[,8:10],1,mean)

# # for spleen :
# cDC1_mean<-apply(expTable[,1:2],1,mean)
# cDC2_mean<-apply(expTable[,3:4],1,mean)
# mac_mean<-apply(expTable[,5:6],1,mean)

# # for thymus :
# cDC1_mean<-apply(expTable[,1:4],1,mean)
# cDC2_mean<-apply(expTable[,5:8],1,mean)
# mac_mean<-apply(expTable[,9:12],1,mean)

# # for lung mac1 :
# cDC1_mean<-apply(expTable[,1:3],1,mean)
# cDC2_mean<-apply(expTable[,4:6],1,mean)
# mac_mean<-apply(expTable[,7:8],1,mean)

# # for lung mac2 :
# cDC1_mean<-apply(expTable[,1:3],1,mean)
# cDC2_mean<-apply(expTable[,4:6],1,mean)
# mac_mean<-apply(expTable[,7:8],1,mean)

# # for lung mac1 and 2 :
# cDC1_mean<-apply(expTable[,1:3],1,mean)
# cDC2_mean<-apply(expTable[,4:6],1,mean)
# mac_mean<-apply(expTable[,7:10],1,mean)

# for all organs :
cDC1_mean<-apply(expTable[,1:12],1,mean)
cDC2_mean<-apply(expTable[,13:25],1,mean)
mac_mean<-apply(expTable[,26:38],1,mean)

expTable_mean<-cbind(cDC1_mean,cDC2_mean,mac_mean)

# for liver:
#colnames(expTable_mean)<-c('Liver_cDC1','Liver_cDC2','Liver_macrophage')

# # for spleen :
# colnames(expTable_mean)<-c('Spleen_cDC1','Spleen_cDC2','Spleen_macrophage')

# # for thymus :
# colnames(expTable_mean)<-c('Thymus_cDC1','Thymus_cDC2','Thymus_macrophage')

# # for lung mac1 :
# colnames(expTable_mean)<-c('Lung_cDC1','Lung_cDC2','Lung_macrophage_1')

# # for lung mac2 :
# colnames(expTable_mean)<-c('Lung_cDC1','Lung_cDC2','Lung_macrophage_2')

# # for lung mac1 and 2 :
# colnames(expTable_mean)<-c('Lung_cDC1','Lung_cDC2','Lung_macrophage_1_and_2')

# for all organs :
colnames(expTable_mean)<-c('All_cDC1','All_cDC2','All_macrophages')

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
                      rmax = 10)#, Goi=IL17genes)
print(p)
saveWidget(p,file="~/Documents/VIB/Projects/Immgen/analysis/interactiveTriwisePlot_filtered_DEgenes_all_organs.html") ##needs full path!
#######################

p<-plotDotplot(theBarycoords, Gdiffexp=allDEgenes, Goi="Clec4f", showlabels = F)
print(p)
ggsave(p,file="results_allSamples/plots_samplesA/after_combat/triwisePlot_IL17.png",dpi=300)

###Rose plot
p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F)
print(p)
ggsave(p,file="results_allSamples/plots_samplesA/after_combat/rosePlot.png",dpi=300)

p<-plotRoseplot(theBarycoords, Gdiffexp=allDEgenes, showlabels = F, Goi=IL17genes)
print(p)
ggsave(p,file="results_allSamples/plots_samplesA/after_combat/rosePlot_IL17.png",dpi=300)

###Interactive plot
p<-interactiveDotplot(expTable_mean, Gdiffexp=allDEgenes, plotLocalEnrichment=FALSE)#, Goi=IL17genes)
print(p)
saveWidget(p,file="~/Documents/VIB/Projects/Immgen/analysis/interactiveTriwisePlot_lung_mac1_and_2.html") ##needs full path!
## Save manually as 'triwisePlot_IL17_highlighted.png'

save(expTable, allDEgenes, file = "~/Documents/VIB/Projects/Immgen/analysis/all_organs_expTable_allDEgenes.RData")
