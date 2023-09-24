rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(doSNOW)
library(foreach)
library(ggplot2)
library(M3C)
library(abind)
library(sparcl)
library(cluster)
library(mclust)

setwd("../../")

source("Scripts/additional_functions_specific_to_comparisons.R")

# Dataset id
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- as.numeric(args[1])

# Model parameters
K <- 100
algo <- "hclust"
linkage <- "complete"
iters <- 25 # default 25, recommended 5-100 (for M3C)
nc_max <- 20

# Printing messages
print(paste0("Clustering method: ", algo))

# Creating folder of public data results
dir.create("Results", showWarnings = FALSE)
dir.create(paste0("Results/Real_data_consensus_unweighted_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Real_data_consensus_unweighted_", algo, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)
print(paste0("Dataset: ", dataset_id))

# Loading the data
xdata <- readRDS(paste0("Data/ICU_datasets/prepared/data_", dataset_id, "_x.rds"))
ydata <- readRDS(paste0("Data/ICU_datasets/prepared/data_", dataset_id, "_y.rds"))
simul <- list()
simul$data <- scale(xdata)
simul$theta <- ydata
n <- table(simul$theta)

# Clustering with G* (unweighted)
if (algo == "hclust") {
  tmptime <- system.time({
    mydist <- dist(simul$data)
    myhclust <- hclust(d = mydist, method = "complete")
    myclusters <- cutree(myhclust, k = length(n))
  })
}
nperf <- c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Silhouette score
silhouette <- InternalCalibration(
  xdata = simul$data,
  method = algo, linkage = linkage,
  index = "silhouette"
)
id <- ManualArgmaxId(silhouette)
if (algo == "hclust") {
  myclusters <- cutree(myhclust, k = id)
}
if (algo == "pam") {
  myclusters <- cluster::pam(x = mydistance, k = id, diss = TRUE, cluster.only = TRUE)
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# CH score
ch <- InternalCalibration(
  xdata = simul$data,
  method = algo, linkage = linkage,
  index = "ch"
)
id <- ManualArgmaxId(ch)
if (algo == "hclust") {
  myclusters <- cutree(myhclust, k = id)
}
if (algo == "pam") {
  myclusters <- cluster::pam(x = mydistance, k = id, diss = TRUE, cluster.only = TRUE)
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# DB score
db <- InternalCalibration(
  xdata = simul$data,
  method = algo, linkage = linkage,
  index = "db"
)
id <- ManualArgmaxId(-db)
if (algo == "hclust") {
  myclusters <- cutree(myhclust, k = id)
}
if (algo == "pam") {
  myclusters <- cluster::pam(x = mydistance, k = id, diss = TRUE, cluster.only = TRUE)
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# GAP statistic
set.seed(1)
tmptime <- system.time({
  out <- GapStatistic(
    xdata = simul$data,
    method = algo, linkage = linkage
  )
})
gap <- out$gap
id <- ManualArgmaxId(gap)
if (algo == "hclust") {
  myclusters <- cutree(myhclust, k = id)
}
if (algo == "pam") {
  myclusters <- cluster::pam(x = mydistance, k = id, diss = TRUE, cluster.only = TRUE)
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Consensus clustering (unweighted)
tmptime <- system.time({
  stab <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    verbose = FALSE,
    nc = 1:nc_max,
  )
})

# Saving the stability object
saveRDS(stab, paste0(filepath, "Stability_unweighted_", dataset_id, ".rds"))

# Consensus clustering with G*
nperf <- rbind(
  nperf,
  c(
    G = length(n),
    ClusteringPerformance(theta = Clusters(stab, argmax_id = length(n)), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Silhouette score
silhouette <- InternalCalibration(
  xdata = simul$data, stability = stab,
  index = "silhouette"
)
id <- ManualArgmaxId(silhouette)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# CH score
ch <- InternalCalibration(
  xdata = simul$data, stability = stab,
  index = "ch"
)
id <- ManualArgmaxId(ch)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# DB score
db <- InternalCalibration(
  xdata = simul$data, stability = stab,
  index = "db"
)
id <- ManualArgmaxId(-db)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# GAP statistic
gap <- InternalCalibration(
  xdata = simul$data, stability = stab,
  index = "gap"
)
id <- ManualArgmaxId(gap)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Delta score
delta <- DeltaAreaCDF(stab)
id <- ManualArgmaxId(delta)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# PAC score
pac <- PAC(stab)
id <- ManualArgmaxId(-pac)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# PINS discrepancy
discrepancy <- PINSDiscrepancy(
  x = simul$data, stab,
  method = algo, linkage = linkage
)
id <- ManualArgmaxId(discrepancy)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# M3C score (PAC)
tmptime2 <- system.time({
  scores <- MonteCarloScore(
    x = simul$data, stab,
    method = algo, linkage = linkage,
    objective = "PAC", iters = iters
  )
})
rcsi_pac <- scores$RCSI # criterion to define assignments in their code
id <- ManualArgmaxId(rcsi_pac)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = ifelse(scores$BETA_P[id] < 0.05, yes = 1, no = 0),
    time = as.numeric(tmptime2[1] + (K - 2) / K * tmptime[1])
  )
)

# M3C score (entropy)
tmptime3 <- system.time({
  scores <- MonteCarloScore(
    x = simul$data, stab,
    method = algo, linkage = linkage,
    objective = "entropy", iters = iters
  )
})
rcsi_entropy <- scores$RCSI # criterion to define assignments in their code
id <- ManualArgmaxId(rcsi_entropy)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul$theta),
    signif = ifelse(scores$NORM_P[id] < 0.05, yes = 1, no = 0),
    time = as.numeric(tmptime3[1] + (K - 2) / K * tmptime[1])
  )
)

# Consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    ClusteringPerformance(theta = stab, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap",
  "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db", "consensus_gap",
  "delta", "pac", "pins_discrepancy", "rcsi_pac", "rcsi_entropy", "consensus_score"
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", dataset_id, ".rds"))
