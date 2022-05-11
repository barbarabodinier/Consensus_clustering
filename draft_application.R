rm(list=ls())
setwd("~/Dropbox/Consensus_clustering")

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(abind)
library(rCOSA)
library(plotrix)

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for (name in r) eval(parse(text = paste0(name, "<-sharp:::", name)))

# Loading all additional functions
myfunctions <- list.files("Scripts/Functions/")
myfunctions <- myfunctions[myfunctions != "Former"]
for (k in 1:length(myfunctions)) {
  source(paste0("Scripts/Functions/", myfunctions[k]))
}

source("Scripts/additional_functions_specific_to_comparisons.R")

# Model parameters
K <- 100
nc_max <- 20
noit <- 5
niter <- 10

mydata=readRDS("Data/ukb_extracted.rds")
mydata=mydata[,grep("\\.0\\.0", colnames(mydata))]
colnames(mydata)=gsub("\\.0\\.0", "", colnames(mydata))
# codes_centre=data.table::fread("Data/codes_10.txt", data.table = FALSE)
# city_not_england=c("Cardiff", "Glasgow", "Edinburgh", "")

features=c("no2","nox","pm10","pm25","pm25_abs", "pm", "noise_day","greenspace", "garden", "water", "natural")
xdata=mydata[,features]
for (j in 1:ncol(xdata)){
  xdata[,j]=log(xdata[,j])
  xdata[which(is.infinite(xdata[,j])),j]=NA
}
xdata=na.exclude(xdata)
xdata=xdata[sample(nrow(xdata), size = 500),]

mydata=mydata[rownames(xdata),]

stab=Clustering(xdata=xdata, 
                nc=1:nc_max,
                implementation=HierarchicalClustering)
CalibrationCurve(stab)
Argmax(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))
table(Clusters(stab))
table(Clusters(stab), mydata$urban)

stab=Clustering(xdata=xdata, 
                K=K,
                nc=1:nc_max,
                implementation=COSAClustering,
                Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
                noit = noit,
                niter = niter)
CalibrationCurve(stab)
Argmax(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))
par(mar=c(21,5,1,1))
WeightBoxplot(stab, boxwex = 0.05)
table(Clusters(stab))
table(Clusters(stab), mydata$urban)
table(Clusters(stab), mydata$centre)





mydata=readRDS("Data/Covid_dataset_baseline_final_selection.rds")
xdata=mydata[,c(23:26)]
xdata=na.exclude(xdata)
xdata=xdata[sample(nrow(xdata), size=300),]
stab=Clustering(xdata=xdata, 
                nc=1:nc_max,
                implementation=HierarchicalClustering)
CalibrationCurve(stab)
Argmax(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))

stab=Clustering(xdata=xdata, 
                K=K,
                nc=1:nc_max,
                implementation=COSAClustering,
                Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
                noit = noit,
                niter = niter)
CalibrationCurve(stab)
Argmax(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))
par(mar=c(21,5,1,1))
WeightBoxplot(stab, boxwex = 0.05)


mydata=readRDS("Data/NOWAC_MTT_smoking_159.rds")
nc_max=50
stab=Clustering(xdata=mydata, 
                nc=1:nc_max,
                implementation=HierarchicalClustering)
CalibrationCurve(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))



# Loading the data
proteins=readRDS("Data/Proteins_selected_denoised_re.rds")
covars=readRDS("Data/Covariates_selected_proteins.rds")
covars=covars[rownames(proteins),]
print(all(rownames(proteins)==rownames(covars)))

metab=readRDS("Data/Metabolites_positive_imputed_summarised.rds")
ids=intersect(rownames(metab), rownames(covars))
metab=metab[ids,]
covars=covars[ids,]

source("~/Dropbox/Metabolomics_lung_cancer/Scripts/functions.R")
covars$smoking_binary=ifelse(covars$smoking_status=="Current", yes=1, no=0)
mymodels=RunLogistic(covars=covars, 
                     predictors=metab, 
                     outcome="smoking_binary", 
                     conf=c("age.sample", "gender"))
ids=rownames(mymodels)[which(as.numeric(mymodels$pval)<0.001)]
nc_max=50
stab=Clustering(xdata=metab[,ids], 
                nc=1:nc_max,
                implementation=HierarchicalClustering)
CalibrationCurve(stab)
pheatmap::pheatmap(ConsensusMatrix(stab, argmax_id = 3))
table(Clusters(stab), covars$smoking_status)

# covars$packyears=scale(covars$packyears)
# covars$smok_intensity=scale(covars$smok_intensity)
# covars$smok_duration=scale(covars$smok_duration)
# covars$CSI=scale(covars$CSI)
# covars$age.sample=scale(covars$age.sample)
# covars$bmi=scale(covars$bmi)

proteins=proteins[covars$gender=="Female",]
covars=covars[rownames(proteins),]

proteins=proteins[covars$case=="0",]
covars=covars[rownames(proteins),]

# Consensus clustering
stab=Clustering(xdata=proteins, 
                nc=1:nc_max,
                implementation=HierarchicalClustering)
CalibrationCurve(stab)
pheatmap::pheatmap(ConsensusMatrix(stab))
table(Clusters(stab), covars$smoking_status)
table(Clusters(stab), covars$case)

# Consensus weighted clustering
stab=Clustering(xdata=proteins, 
                K=K,
                nc=1:nc_max,
                implementation=COSAClustering,
                Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
                noit = noit,
                niter = niter)

CalibrationCurve(stab)
pheatmap::pheatmap(ConsensusMatrix(stab, argmax_id = 4))







mysamples=data.table::fread("Data/GSE3790-GPL96_series_matrix.txt", skip = 51)
mysamples=mysamples[c(1,7,9),]

mydata=data.table::fread("Data/GSE3790-GPL96_series_matrix.txt", skip = 83)


# mydata=readRDS("Data/dataset_LA5.rds")
# ydata=mydata[,1]
# xdata=mydata[,-1]
# # xdata=xdata[,which(apply(xdata,2,sd)>0.5)]
# stab=Clustering(xdata=xdata, 
#                 nc=1:20,
#                 implementation=HierarchicalClustering)
# CalibrationCurve(stab)
# pheatmap::pheatmap(ConsensusMatrix(stab))
# table(Clusters(stab), ydata)
# mydata=readRDS("Data/dataset_LA5.rds")
# ydata=mydata[,1]
# xdata=mydata[,-1]
# xdata=xdata[,which(apply(xdata,2,sd)>1.5)]
# stab=Clustering(xdata=xdata, 
#                 nc=1:20,
#                 implementation=COSAClustering,
#                 Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 5),
#                 noit = noit,
#                 niter = niter)
# CalibrationCurve(stab)
# pheatmap::pheatmap(ConsensusMatrix(stab))
# table(Clusters(stab), ydata)
# median_weights=apply(stab$Beta[ArgmaxId(stab)[1],,],1,median)
# ids=sort.list(median_weights)
# WeightBoxplot(stab, at=ids)

