library(fake)
library(sharp)
library(aricode)
library(M3C)
library(abind)
library(cluster)
library(sparcl)
library(rCOSA)

setwd("../../")

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation study parameters
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
params_id <- as.numeric(args[2])
seed <- as.numeric(args[3])
simulation_id <- paste0(params_id, "_", seed)

# Extracting simulation parameters
params_list <- read.table(paste0("Scripts/Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Model parameters
K <- 100
nc_max <- 20
noit <- 20
niter <- 10
n_lambda <- 10

# Extracting simulation parameters
nc <- params_list[params_id, "nc"]
equal_size <- params_list[params_id, "equal_size"]
n_tot <- params_list[params_id, "n_tot"]
p <- params_list[params_id, "p"]
ev_xc <- params_list[params_id, "ev_xc"]
nu_xc <- params_list[params_id, "nu_xc"]
v_min <- params_list[params_id, "v_min"]
v_max <- params_list[params_id, "v_max"]

# Printing messages
print(paste0("Simulation study ", simul_study_id))
print(paste0("Number of items: ", n_tot))
print(paste0("Number of features: ", p))
print(paste0("Explained variance per feature: ", ev_xc))
print(paste0("Proportion of contributing features: ", nu_xc))
print(paste0("v_min: ", v_min))
print(paste0("v_max: ", v_max))
print(paste("Simulation ID:", simulation_id))

# Creating folder of simulation study
dir.create("Results", showWarnings = FALSE)
dir.create("Results/Simulations_consensus_weighted_hierarchical", showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_weighted_hierarchical/Simulations_", simul_study_id, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Data simulation
set.seed(seed)
if (equal_size) {
  n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
} else {
  n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
}
pk <- round(rep(0.2, 5) * p)
q <- round(nu_xc * p)
theta_xc <- c(rep(1, q), rep(0, p - q))
simul <- SimulateClustering(
  n = n,
  pk = pk,
  ev_xc = ev_xc,
  nu_within = 1,
  nu_between = 0,
  v_within = c(v_min, v_max),
  v_between = 0,
  v_sign = -1,
  pd_strategy = "min_eigenvalue",
  theta_xc = theta_xc,
  output_matrices = TRUE
)
simul$data <- scale(simul$data)

# Hierarchical clustering with G*
tmptime <- system.time({
  mydist <- dist(simul$data)
  myhclust <- hclust(d = mydist, method = "complete")
  myclusters <- cutree(myhclust, k = length(n))
})
nperf <- c(
  G = length(n),
  lambda = "",
  q = "",
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Hierarchical clustering with max silhouette score
silhouette <- SilhouetteScore(x = simul$data, method = "hclust")
id <- ManualArgmaxId(silhouette)
myclusters <- cutree(myhclust, k = id)
nperf <- rbind(
  nperf,
  c(
    G = id,
    lambda = "",
    q = "",
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Hierarchical clustering with max GAP statistic
tmptime <- system.time({
  out <- GapStatistic(xdata = simul$data, method = "hclust")
})
gap <- out$gap
id <- ManualArgmaxId(gap)
myclusters <- cutree(myhclust, k = id)
nperf <- rbind(
  nperf,
  c(
    G = id,
    lambda = "",
    q = "",
    ClusteringPerformance(theta = myclusters, theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = length(n)), theta_star = simul),
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
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Consensus clustering (sparcl)
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

# Consensus clustering with G* (sparcl)
for (argmax_id in which(stab$nc == length(n))) {
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = stab$Lambda[argmax_id],
      q = stab$Q[argmax_id],
      ClusteringPerformance(theta = Clusters(stab, argmax_id = argmax_id), theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Selection with G* (sparcl)
selperf <- NULL
for (argmax_id in which(stab$nc == length(n))) {
  selprop <- stab$selprop[argmax_id, ]
  median_weights <- apply(stab$Beta[argmax_id, , ], 1, median)
  tmp <- cbind(selprop, median_weights)
  tmp <- tmp[order(selprop, median_weights, decreasing = TRUE), ]
  selected <- rep(0, ncol(stab$Beta))
  names(selected) <- colnames(stab$Beta)
  selected[rownames(tmp)[1:q]] <- 1
  selperf <- rbind(
    selperf, SelectionPerformance(
      theta = selected,
      theta_star = theta_xc
    )
  )
}

# Consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    lambda = Argmax(stab)[2],
    q = stab$Q[ArgmaxId(stab)[1]],
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Selection with consensus score (sparcl)
argmax_id <- ArgmaxId(stab)[1]
selprop <- stab$selprop[argmax_id, ]
median_weights <- apply(stab$Beta[argmax_id, , ], 1, median)
tmp <- cbind(selprop, median_weights)
tmp <- tmp[order(selprop, median_weights, decreasing = TRUE), ]
selected <- rep(0, ncol(stab$Beta))
names(selected) <- colnames(stab$Beta)
selected[rownames(tmp)[1:q]] <- 1
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = selected,
    theta_star = theta_xc
  )
)

# Consensus clustering (cosa)
tmptime <- system.time({
  stab <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    verbose = FALSE,
    nc = 1:nc_max,
    Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = n_lambda),
    noit = noit,
    niter = niter
  )
})

# Consensus clustering with G* (cosa)
for (argmax_id in which(stab$nc == length(n))) {
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = stab$Lambda[argmax_id],
      q = "",
      ClusteringPerformance(theta = Clusters(stab, argmax_id = argmax_id), theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Selection with G* (cosa)
for (argmax_id in which(stab$nc == length(n))) {
  median_weights <- apply(stab$Beta[argmax_id, , ], 1, median)
  selected <- rep(0, ncol(stab$Beta))
  names(selected) <- colnames(stab$Beta)
  selected[names(sort(median_weights, decreasing = TRUE))[1:q]] <- 1
  selperf <- rbind(
    selperf, SelectionPerformance(
      theta = selected,
      theta_star = theta_xc
    )
  )
}

# Consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    lambda = Argmax(stab)[2],
    q = "",
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Selection with consensus score
argmax_id <- ArgmaxId(stab)[1]
median_weights <- apply(stab$Beta[argmax_id, , ], 1, median)
selected <- rep(0, ncol(stab$Beta))
names(selected) <- colnames(stab$Beta)
selected[names(sort(median_weights, decreasing = TRUE))[1:q]] <- 1
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = selected,
    theta_star = theta_xc
  )
)

# Re-formatting output objects
rownames(nperf) <- c(
  "hclust_star", "hclust_silhouette", "hclust_gap",
  "unweighted_star", "unweighted",
  paste0("sparcl_star_", 1:n_lambda), "sparcl",
  paste0("cosa_star_", 1:n_lambda), "cosa"
)
rownames(selperf) <- c(
  paste0("sparcl_star_", 1:n_lambda), "sparcl",
  paste0("cosa_star_", 1:n_lambda), "cosa"
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
saveRDS(selperf, paste0(filepath, "Selection_performances_", simulation_id, ".rds"))
