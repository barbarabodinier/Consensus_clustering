library(fake)
library(rCOSA)
library(sharp)
library(abind)
library(sparcl)
library(cluster)

setwd("../../")

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation study parameters
args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
params_id <- as.numeric(args[2])
seed <- as.numeric(args[3])
simulation_id <- paste0(params_id, "_", seed)
algo <- as.character(args[4])

# Extracting simulation parameters
params_list <- read.table(paste0("Scripts/Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Model parameters
K <- 100
linkage <- "complete"
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
print(paste0("Clustering method: ", algo))

# Creating folder of simulation study
dir.create("Results", showWarnings = FALSE)
dir.create(paste0("Results/Simulations_consensus_sparcl_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_sparcl_", algo, "/Simulations_", simul_study_id, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Data simulation
set.seed(seed)
if (equal_size) {
  n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
} else {
  if (nc == 5) {
    n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
  } else {
    n <- round(c(500, 300, 150, 50) / sum(c(500, 300, 150, 50)) * n_tot)
  }
}
pk <- round(rep(0.2, 5) * p)
q <- round(nu_xc * p)
theta_xc <- c(rep(1, q), rep(0, p - q))
sigma <- SimulateCorrelation(
  pk = pk,
  nu_within = 1,
  nu_between = 0,
  v_within = c(v_min, v_max),
  v_between = 0,
  v_sign = -1,
  pd_strategy = "min_eigenvalue"
)$sigma
simul <- SimulateClustering(
  n = n,
  pk = pk,
  sigma = sigma,
  ev_xc = ev_xc,
  theta_xc = theta_xc,
  output_matrices = TRUE
)
simul$data <- scale(simul$data)

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
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
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
      ClusteringPerformance(theta = myclusters, theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Selection performance with G* and different lambdas
selperf <- NULL
for (i in 1:n_lambda) {
  median_weights <- single_run$weight[i, ]
  selected <- rep(0, ncol(single_run$weight))
  names(selected) <- colnames(single_run$weight)
  selected[names(sort(median_weights, decreasing = TRUE))[1:q]] <- 1
  selperf <- rbind(
    selperf, SelectionPerformance(
      theta = selected,
      theta_star = theta_xc
    )
  )
}

# Sparse clustering with G* and lambda calibrated by gap
Lambda <- LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
if (algo == "hclust") {
  tmptime <- system.time({
    single_run <- SparseHierarchicalClustering(
      xdata = simul$data,
      nc = length(n),
      Lambda = Lambda,
      use_permutations = TRUE
    )
  })
}
if (algo == "kmeans") {
  tmptime <- system.time({
    single_run <- SparseKMeansClustering(
      xdata = simul$data,
      nc = length(n),
      Lambda = Lambda,
      use_permutations = TRUE
    )
  })
}

# Clustering performance with G* and lambda calibrated by gap
i <- 1
pseudo_dist <- as.dist(1 - single_run$comembership[, , i])
myclusters <- cutree(hclust(d = pseudo_dist), k = length(n))
nperf <- rbind(
  nperf,
  c(
    G = length(n),
    lambda = single_run$Lambda,
    q = sum(single_run$weight != 0),
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Selection performance with G* and lambda calibrated by gap
i <- 1
median_weights <- single_run$weight[i, ]
selected <- rep(0, ncol(single_run$weight))
names(selected) <- colnames(single_run$weight)
selected[names(sort(median_weights, decreasing = TRUE))[1:q]] <- 1
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = selected,
    theta_star = theta_xc
  )
)

# Consensus clustering (unweighted)
if (algo == "hclust") {
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
if (algo == "kmeans") {
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = length(n)), theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)
selperf <- rbind(selperf, NA)

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
selperf <- rbind(selperf, NA)

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
      ClusteringPerformance(theta = Clusters(stab, argmax_id = argmax_id), theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Selection with G*
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

# Consensus sparse clustering with lambda calibrated by gap
if (algo == "hclust") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = SparseHierarchicalClustering,
      use_permutations = TRUE,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
      Lambda = LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
    )
  })

  # Fixing the output as sharp is not designed for this
  stab$coprop <- stab$coprop[, , 1:nc_max]
  stab$nc <- cbind(1:nc_max)
  stab$Lambda <- cbind(rep(NA, nc_max))
}
if (algo == "kmeans") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = SparseKMeansClustering,
      use_permutations = TRUE,
      K = K,
      verbose = FALSE,
      nc = 2:nc_max,
      Lambda = LambdaSequence(lmax = 10, lmin = 1.1, cardinal = n_lambda)
    )
  })

  # Fixing the output as sharp is not designed for this
  stab$coprop <- stab$coprop[, , 1:nc_max]
  stab$nc <- cbind(1:nc_max)
  stab$Lambda <- cbind(rep(NA, nc_max))
}

# Consensus sparse clustering with G*
for (argmax_id in which(stab$nc == length(n))) {
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = "",
      q = stab$Q[argmax_id],
      ClusteringPerformance(theta = Clusters(stab, argmax_id = argmax_id), theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Selection with G*
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
stab <- stab_saved
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

# Selection with consensus score
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

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star",
  paste0("single_run_sparcl_star_", 1:n_lambda),
  "single_run_sparcl_star_gap",
  "consensus_star", "consensus",
  paste0("consensus_sparcl_star_", 1:n_lambda),
  "consensus_sparcl_star_gap",
  "consensus_sparcl"
)
rownames(selperf) <- c(
  paste0("single_run_sparcl_star_", 1:n_lambda),
  "single_run_sparcl_star_gap",
  "consensus_star", "consensus",
  paste0("consensus_sparcl_star_", 1:n_lambda),
  "consensus_sparcl_star_gap",
  "consensus_sparcl"
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
saveRDS(selperf, paste0(filepath, "Selection_performances_", simulation_id, ".rds"))
