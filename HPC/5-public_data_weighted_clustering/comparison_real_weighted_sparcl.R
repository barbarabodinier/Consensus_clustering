rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(abind)
library(sparcl)
library(cluster)

setwd("../../")

source("Scripts/additional_functions_specific_to_comparisons.R")

# Dataset id
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- as.numeric(args[1])
algo=as.character(args[2])

# Model parameters
K <- 100
linkage <- "complete"
nc_max <- 20
noit <- 20
niter <- 10
n_lambda <- 10

# Printing messages
print(paste0("Clustering method: ", algo))

# Creating folder of public data results
dir.create("Results", showWarnings = FALSE)
dir.create(paste0("Results/Real_data_consensus_sparcl_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Real_data_consensus_sparcl_", algo, "/")
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
if (algo == "kmeans") {
  tmptime <- system.time({
    set.seed(1)
    mykmeans <- stats::kmeans(x = simul$data, centers = length(n))
    myclusters <- mykmeans$cluster
  })
}
nperf <- c(
  G = length(n),
  lambda = "",
  q = "",
  ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# GAP statistic
set.seed(1)
tmptime <- system.time({
  out <- GapStatistic(
    xdata = simul$data,
    method = algo,
    linkage = linkage,
    nc_max = nc_max
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
if (algo == "kmeans") {
  set.seed(1)
  mykmeans <- stats::kmeans(x = simul$data, centers = id)
  myclusters <- mykmeans$cluster
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    lambda = "",
    q = "",
    ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Sparse clustering with G* and different lambdas
Lambda <- LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
if (algo == "hclust") {
  tmptime <- system.time({
    single_run <- SparseHierarchicalClustering(
      xdata = simul$data,
      nc = length(n),
      Lambda = Lambda
    )
  })
}
if (algo == "kmeans") {
  tmptime <- system.time({
    single_run <- SparseKMeansClustering(
      xdata = simul$data,
      nc = length(n),
      Lambda = Lambda
    )
  })
}

# Clustering performance with G* and different lambdas
for (i in 1:n_lambda) {
  pseudo_dist <- as.dist(1 - single_run$comembership[, , i])
  myclusters <- cutree(hclust(d = pseudo_dist), k = length(n))
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = Lambda[i],
      q = sum(single_run$weight[i, ] != 0),
      ClusteringPerformance(theta = myclusters, theta_star = simul$theta),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Consensus clustering (unweighted)
if (algo=="hclust"){
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = HierarchicalClustering,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
    )
  })
}
if (algo=="kmeans") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = KMeansClustering,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
    )
  })
}

# Consensus clustering with G* (unweighted)
nperf <- rbind(
  nperf,
  c(
    G = length(n),
    lambda = "",
    q = "",
    ClusteringPerformance(theta = Clusters(stab, argmax_id = length(n)), theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    lambda = "",
    q = "",
    ClusteringPerformance(theta = stab, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Consensus sparse clustering
if (algo == "hclust") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = SparseHierarchicalClustering,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
      Lambda = LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
    )
  })
}
if (algo == "kmeans") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = SparseKMeansClustering,
      K = K,
      verbose = FALSE,
      nc = 2:nc_max,
      Lambda = LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
    )
  })
}
stab_saved <- stab

# Consensus sparse clustering with G*
for (argmax_id in which(stab$nc == length(n))) {
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = stab$Lambda[argmax_id],
      q = stab$Q[argmax_id],
      ClusteringPerformance(theta = Clusters(stab, argmax_id = argmax_id), theta_star = simul$theta),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Consensus score
stab <- stab_saved
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    lambda = Argmax(stab)[2],
    q = stab$Q[ArgmaxId(stab)[1]],
    ClusteringPerformance(theta = stab, theta_star = simul$theta),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star",
  "single_run_gap",
  paste0("single_run_sparcl_star_", 1:n_lambda),
  "consensus_star", "consensus",
  paste0("consensus_sparcl_star_", 1:n_lambda),
  "consensus_sparcl"
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", dataset_id, ".rds"))
