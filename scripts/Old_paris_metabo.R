##
library(dplyr)
library(openxlsx)
#library(metabomxtr)

library(devtools)
#install_github("xia-lab/MetaboAnalystR")
library(MetaboAnalystR)

## excel file with info about samples
samples <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Data synthesis local cohort Saint-Louis 032018_modified.xlsx",
                     check.names = FALSE) %>%
  mutate(DATEOFCYTOFEXPERIMENT = as.Date(DATEOFCYTOFEXPERIMENT, "%d.%m.%Y"),
         DOB = as.Date(DOB, "%d.%m.%Y"),
         DOG = as.Date(DOG, "%d.%m.%Y"),
         DATEOFSAMPLE = as.Date(DATEOFSAMPLE, "%d.%m.%Y"),
         GROUP = tolower(GROUP)) %>%
  dplyr::filter(!is.na(FCSNAME))
rownames(samples)<- samples$Id.Cryostem.R
samp_recip <- samples[which(rownames(samples) %in% names(recip_names)),]

## strange numbers after row 390, I only read 80 first rows (corresponding to patients)
metabo<- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis_filtered.xlsx",rows = c(4:83))

meta_metabo <- read.xlsx("~/Documents/VIB/Projects/Integrative_Paris/documents_22:02:18/CYTOF_David_Michonneau/Metabolomic local cohort Saint-Louis.xlsx",
                         rows = c(1:2), cols = c(3:842))

sum_nas<-apply(metabo,2,function(x){sum(is.na(x))})
high_na_vars<-which(sum_nas>=40) # 85 columns have more than 50% NA
metabo<-metabo[,-high_na_vars] # remove variables with NAs in more than half patients

met <- apply(metabo,2,function(x){

})

## ------------------- analysis
##recipients:
recip<-metabo[-c(1:34),-which(sum_nas>=25)]
rownames(recip)<-recip$PATIENT.ID
recip$GROUP<-as.factor(recip$GROUP)
recip$GROUP<-relevel(recip$GROUP,ref="RECIPIENT NOTOL") ## set NOTOL as reference group
table(recip$GROUP)
## only two groups -> group tol 1 and tol 2?
recip$GROUP[which(recip$GROUP=="RECIPIENT 1TOL")]<-"RECIPIENT 2TOL"

## mixture model
common<-samples$METABONAME[which(samples$METABONAME%in%recip$PATIENT.ID)]
rectab<-cbind(recip[common,],samples[which(samples$METABONAME%in%recip$PATIENT.ID),colnames(samples)[c(3,10,13:18)]])
rectab[,c("PATIENT.ID","METABONAME")] # control match btw two tables
colnames(samples)[c(3,10,13:18)]
for (i in colnames(samples)[c(3,10,13:18)]){
  print(i)
  rectab[,i]<-as.factor(rectab[,i])
}


fullModel<-~GROUP|GROUP+GENDER+CMVStatus+GROUPE+DONORSEX+DONORCMV+DONORGROUPE+HEMATOLOGICDISORDER
reducedModel<-~1|GENDER+CMVStatus+GROUPE+DONORSEX+DONORCMV+DONORGROUPE+HEMATOLOGICDISORDER
fullModelResults<-mxtrmod(ynames=colnames(rectab)[3:300],mxtrModel=fullModel,data=rectab)
reducedModelResults<-mxtrmod(ynames=colnames(rectab)[3:300],mxtrModel=reducedModel,data=rectab,fullModel=fullModel)
finalResult<-mxtrmodLRT(fullmod=fullModelResults,redmod=reducedModelResults,adj="BH")
finalResult
## ne trouve aucune pvalue ajustee en dessous de 1



## impute missing values ------------------------------
library(impute)
expr<-metabo[,-c(1,2)]

## no column should have more than 80% missing values:
high_na_vars<-which(sum_nas>=40)
expr<-metabo[,-c(1,2,high_na_vars)]
rownames(expr)<-metabo$PATIENT.ID
expr<-t(expr)

iexpr <- cbind(impute.knn(data.matrix(expr))$data,
               stringsAsFactors = FALSE)
head(iexpr[, 1:6])
iexpr<-iexpr[,-80]

pca<-prcomp(t(iexpr),scale. = T)
library(car)
colrs<-metabo$GROUP %>% recode("'DONOR NOTOL'='yellow';'DONOR TOL1'='orange';'DONOR TOL2'='red';'RECIPIENT 1TOL'='green';'RECIPIENT 2TOL'='blue';'RECIPIENT NOTOL'='black'")
plot(pca$x[,1:2], col=colrs)

library(ade4)
acp<-dudi.pca(t(iexpr),scale=T)
scatter(acp)
s.class(acp$li,fac = as.factor(metabo$GROUP),col=rainbow(6)) ## on the right: recip no tol ++
s.class(acp$li,fac = as.factor(metabo$GROUP),col=rainbow(6))
s.arrow(acp$c1, clabel = .3, boxes = F)
s.label(acp$li,label = metabo$PATIENT.ID, boxes=F, clabel = .5)
rownames(acp$c1)[which(acp$c1[,1]==max(acp$c1[,1]))] # myristate in recip no tol ++
rownames(acp$c1)[which(acp$c1[,1]==min(acp$c1[,1]))] # androstenediol in recip no tol --
rownames(acp$c1)[which(acp$c1[,1]>=0.07)] # lipides +++ palmitoyl, oleoyl, laurate, myristate

library(tsne)
library(ggplot2)
tsne1<-tsne(iexpr[,-757],k = 10)
plot(tsne1, col=colrs)

metabo_tsne <- metabo %>% dplyr::mutate(tsne_1 = tsne1[,1],
                                          tsne_2 = tsne1[,2])
ggplot(metabo_tsne) +
  geom_point(aes(x = tsne_1, y = tsne_2, col = as.factor(GROUP)), size = 4) +
  geom_text(aes(x = tsne_1, y = tsne_2, label= PATIENT.ID)) +
  theme_minimal() +
  theme(text = element_text(size = 12)) # tSNE separe pas vraiment patients par groupe


## apprentissage ----------------------------

library(randomForest)
iiexpr<- as.data.frame(t(iexpr)) %>% dplyr::mutate(group = metabo$GROUP)
iiexpr$group<-as.factor(iiexpr$group)

## changer noms pour que rf soit content: trop chiant
# colnames(iiexpr)<-gsub("^[[:digit:]]{+}-","",colnames(iiexpr),ignore.case=T)
# colnames(iiexpr)<-gsub("[(,),.,+,-,:,/]","",colnames(iiexpr),ignore.case=T)
# head(colnames(iiexpr))
#
# bli<-colnames(t(iexpr))[which(duplicated(colnames(iiexpr)))]
# bli2<-bli[grep("^[[:digit:]][-^,]",bli)]
# blii2<-unlist(lapply(bli2,function(x) paste(strsplit(x,"-")[[1]][2],strsplit(x,"-")[[1]][1],sep="")))

colnames(iiexpr)<-paste0("a",c(1:757))
colnames(iiexpr)[757]<-"group"
rf<-randomForest(group ~ ., data=iiexpr, ntree=3000, mtry=300)
var_imp<-round(importance(rf), 2)
colnames(t(iexpr))[which(var_imp==max(var_imp))] ## best var to separate classes
partialPlot(rf, iiexpr, "a331", "RECIPIENT NOTOL")

## pathways ---------------------------------





