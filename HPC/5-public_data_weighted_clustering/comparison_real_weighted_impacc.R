rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(abind)
library(cluster)
library(IMPACC)

setwd("../../")

source("Scripts/additional_functions_specific_to_comparisons.R")

# Dataset id
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- as.numeric(args[1])

# Model parameters
K <- 100
algo <- "hclust"
linkage <- "complete"
nc_max <- 20
noit <- 20
niter <- 10
n_lambda <- 10

# Printing messages
print(paste0("Clustering method: ", algo))

# Creating folder of public data results
dir.create("Results", showWarnings = FALSE)
dir.create(paste0("Results/Real_data_consensus_impacc_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Real_data_consensus_impacc_", algo, "/")
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
rownames(simul$data) <- paste0("obs", 1:nrow(simul$data))
colnames(simul$data) <- paste0("var", 1:ncol(simul$data))

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

# IMPACC with G*
for (pFeature in seq(0.1, 0.9, by = 0.1)) {
  print(pFeature)

  set.seed(1)
  tmptime <- system.time({
    impacc <- tryCatch(
      IMPACC(
        d = t(simul$data),
        K = length(n),
        reps = K,
        pItem = 0.5,
        pFeature = pFeature,
        innerLinkage = linkage,
        distance = "euclidean",
        finalAlgorithm = "hclust",
        finalLinkage = "complete",
        verbose = FALSE
      ),
      error = function(e) {
        message("Not run.")
      }
    )
  })

  # Clustering performance with G*
  if (!is.null(impacc)) {
    myclusters <- impacc$labels
    nperf <- rbind(
      nperf,
      c(
        G = length(n),
        lambda = NA,
        q = "",
        ClusteringPerformance(theta = as.numeric(myclusters), theta_star = simul$theta),
        signif = NA,
        time = as.numeric(tmptime[1])
      )
    )
  } else {
    nperf <- rbind(nperf, NA)
  }
}

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star",
  "single_run_gap",
  "consensus_star", "consensus",
  paste0("impacc_", seq(1, 9))
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", dataset_id, ".rds"))
