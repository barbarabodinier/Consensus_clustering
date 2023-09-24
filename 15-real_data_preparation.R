rm(list = ls())

library(datasetsICR)
library(data.table)
library(sharp)


# https://archive.ics.uci.edu/dataset/53/iris

## Loading the data
x <- scale(iris[, 1:4])
y <- iris[, 5]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_1_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_1_y.rds")


# Seeds data also in https://archive.ics.uci.edu/dataset/236/seeds

## Loading the data
data("seeds")
mydata <- seeds

## Extracting attributes and target
y <- mydata$variety
x <- mydata[, 1:7]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_2_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_2_y.rds")


# Wine data also in https://archive.ics.uci.edu/dataset/109/wine

## Loading the data
data(wine)
mydata <- wine

## Extracting attributes and target
y <- mydata[, 1]
x <- mydata[, -1]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_3_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_3_y.rds")


# R package palmerpenguins

## Loading the data
library(palmerpenguins)

## Extracting attributes and target
attribute_names <- c(
  "bill_length_mm",
  "bill_depth_mm",
  "flipper_length_mm",
  "body_mass_g"
)
x <- as.matrix(penguins[, attribute_names])
y <- penguins$species

## Excluding missing values (N=2)
ids <- which(apply(x, 1, FUN = function(z) {
  !any(is.na(z))
}))
x <- x[ids, ]
y <- y[ids]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_4_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_4_y.rds")


# Hawks data

## Loading the data
library(Stat2Data)
data("Hawks")
mydata <- Hawks

## Extracting attributes and target
y <- mydata$Species
x <- mydata[, c(10:13, 15)]

## Excluding missing values
ids <- which(apply(x, 1, FUN = function(z) {
  !any(is.na(z))
}))
x <- x[ids, ]
y <- y[ids]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_5_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_5_y.rds")


# Microarray data

## Loading the data and preparing the data
mydata <- fread("Data/Lung_cancer/fig1tree.cdt.tsv", data.table = FALSE)
mydata <- mydata[-c(1, 2), ]
rownames(mydata) <- mydata[, 2]
mydata <- mydata[, -c(1:4)]
for (k in 1:ncol(mydata)) {
  mydata[, k] <- as.numeric(mydata[, k])
}
x <- t(mydata)

## Checking median transformation of items
apply(mydata, 2, median)
apply(mydata, 2, FUN = function(z) {
  sum(z^2)
})

## Removing N=2 outliers
x <- x[-c(9, 166), ]

## Defining the classes
covars <- data.frame(
  normal = ifelse(grepl("^NL", rownames(x)), yes = 1, no = 0),
  adenocarcinoma = ifelse(grepl("^AD", rownames(x)), yes = 1, no = 0),
  carcinoid = ifelse(grepl("^COID", rownames(x)), yes = 1, no = 0),
  small_cell = ifelse(grepl("^SMCL", rownames(x)), yes = 1, no = 0),
  squamous = ifelse(grepl("^SQ", rownames(x)), yes = 1, no = 0)
)
subtypes <- sharp:::DummyToCategories(covars)
names(subtypes) <- rownames(x)
print(all(apply(covars, 1, sum) == 1))

## Excluding adenocarcinomas
x <- x[which(subtypes != 2), ]
y <- subtypes[which(subtypes != 2)]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_6_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_6_y.rds")


# https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq

## Loading gene expression data
mydata <- read.table("Data/ICU_datasets/raw/TCGA-PANCAN-HiSeq-801x20531/data.csv",
  sep = ",", header = TRUE
)
rownames(mydata) <- mydata[, 1]
mydata <- mydata[, -1]

## RNAseq data preparation
nvar <- 1000
x <- log2(mydata + 1)
myvariances <- apply(x, 2, var)
x <- x[, sort.list(myvariances, decreasing = TRUE)[1:nvar]]
x <- scale(x)

## Loading labels
mylabels <- read.table("Data/ICU_datasets/raw/TCGA-PANCAN-HiSeq-801x20531/labels.csv",
  sep = ",", header = TRUE
)
y <- mylabels[, 2]
names(y) <- mylabels[, 1]

## Excluding breast cancer samples
ids <- which(y == "BRCA")
y <- y[-ids]
x <- x[-ids, ]

## Random sampling of observations
set.seed(1)
nobs <- 200
row_ids <- Resample(data = y, family = "multinomial", tau = nobs / nrow(x))
x <- x[row_ids, ]
y <- y[row_ids]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_7_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_7_y.rds")


# Single cell Tabula Muris

## Loading the data
library(scCCESS)
data(sce, package = "scCCESS")

## Extracting attributes and target
x <- SingleCellExperiment::counts(sce)
x <- t(prefilter(x))
xvar <- apply(x, 2, var)
n_genes <- 1000
x <- x[, sort.list(xvar, decreasing = TRUE)[1:n_genes]]
y <- sce$cellTypes[rownames(x)]

## Excluding immune cells
ids <- which(y %in% c(
  "basal cell",
  "enterocyte of epithelium of large intestine",
  "fibroblast",
  "keratinocyte stem cell",
  "myeloid cell",
  "skeletal muscle satellite cell"
))
y <- y[ids]
x <- x[ids, ]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_8_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_8_y.rds")


# Single cell data

## Loading the data
mydata <- readRDS("Data/ICU_datasets/other/single_cell_10x_5cl.rds")

## Extracting attributes and target
y <- mydata[, 1]
x <- mydata[, 3:ncol(mydata)]

## Formatting the data
x <- as.matrix(x)
rownames(x) <- names(y) <- paste0("obs", 1:nrow(x))

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_9_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_9_y.rds")


# Other single cell data

## Loading the data
library(SC3)
x <- t(yan)
y <- as.character(ann[, 1])
x <- log2(x + 1)

xvar <- apply(x, 2, var)
n_genes <- 5000
x <- x[, sort.list(xvar, decreasing = TRUE)[1:n_genes]]

## Saving prepared datasets
saveRDS(x, "Data/ICU_datasets/prepared/data_10_x.rds")
saveRDS(y, "Data/ICU_datasets/prepared/data_10_y.rds")
