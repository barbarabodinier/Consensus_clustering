library(fake)
library(rCOSA)
library(sharp)
library(abind)
library(IMPACC)

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
dir.create(paste0("Results/Simulations_consensus_impacc_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_impacc_", algo, "/Simulations_", simul_study_id, "/")
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
nperf <- c(
  G = length(n),
  lambda = "",
  q = "",
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
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

# IMPACC with G*
selperf <- NULL
for (pFeature in seq(0.1, 0.9, by = 0.1)) {
  set.seed(seed)
  tmptime <- system.time({
    impacc <- IMPACC(
      d = t(simul$data),
      K = length(n),
      reps = K,
      pItem = 0.5,
      pFeature = pFeature,
      innerLinkage = linkage,
      distance = "euclidean",
      finalAlgorithm = algo,
      finalLinkage = "complete",
      verbose = FALSE
    )
  })

  # Clustering performance with G*
  myclusters <- impacc$labels
  nperf <- rbind(
    nperf,
    c(
      G = length(n),
      lambda = NA,
      q = "",
      ClusteringPerformance(theta = as.numeric(myclusters), theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )

  # Selection with G*
  median_weights <- impacc$feature_importance
  selected <- rep(0, ncol(simul$data))
  names(selected) <- colnames(simul$data)
  selected[names(sort(median_weights, decreasing = TRUE))[1:q]] <- 1
  selperf <- rbind(
    selperf,
    SelectionPerformance(
      theta = selected,
      theta_star = theta_xc
    )
  )
}

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star",
  "consensus_star", "consensus",
  paste0("impacc_", seq(1, 9))
)
rownames(selperf) <- c(
  paste0("impacc_", seq(1, 9))
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
saveRDS(selperf, paste0(filepath, "Selection_performances_", simulation_id, ".rds"))
