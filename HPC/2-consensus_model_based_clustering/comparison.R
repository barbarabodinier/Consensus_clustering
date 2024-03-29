setwd("../../")

library(fake)
library(sharp)
library(doSNOW)
library(foreach)
library(ggplot2)
library(M3C)
library(abind)
library(cluster)
library(mclust)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Checking package version
print(packageVersion("M3C"))
print(packageVersion("fake"))
print(packageVersion("sharp"))

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
iters <- 25 # default 25, recommended 5-100 (for M3C)
nc_max <- 20 # maximum number of clusters

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
dir.create(paste0("Results/Simulations_consensus_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_", algo, "/Simulations_", simul_study_id, "/")
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
  nu_xc = nu_xc,
  output_matrices = TRUE
)
simul$data <- scale(simul$data)

# Clustering with G*
if (algo == "kmeans") {
  tmptime <- system.time({
    set.seed(1)
    mykmeans <- stats::kmeans(x = simul$data, centers = length(n))
    myclusters <- mykmeans$cluster
  })
}
if (algo == "gmm") {
  tmptime <- system.time({
    myclust <- mclust::Mclust(data = simul$data, G = length(n), verbose = FALSE)
    myclusters <- myclust$classification
  })
}
nperf <- c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Silhouette score
silhouette <- InternalCalibration(
  xdata = simul$data,
  method = algo,
  index = "silhouette"
)
id <- ManualArgmaxId(silhouette)
if (algo == "kmeans") {
  set.seed(1)
  mykmeans <- stats::kmeans(x = simul$data, centers = id)
  myclusters <- mykmeans$cluster
}
if (algo == "gmm") {
  myclust <- mclust::Mclust(data = simul$data, G = id, verbose = FALSE)
  myclusters <- myclust$classification
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# CH score
ch <- InternalCalibration(
  xdata = simul$data,
  method = algo,
  index = "ch"
)
id <- ManualArgmaxId(ch)
if (algo == "kmeans") {
  set.seed(1)
  mykmeans <- stats::kmeans(x = simul$data, centers = id)
  myclusters <- mykmeans$cluster
}
if (algo == "gmm") {
  myclust <- mclust::Mclust(data = simul$data, G = id, verbose = FALSE)
  myclusters <- myclust$classification
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# DB score
db <- InternalCalibration(
  xdata = simul$data,
  method = algo,
  index = "db"
)
id <- ManualArgmaxId(-db)
if (algo == "kmeans") {
  set.seed(1)
  mykmeans <- stats::kmeans(x = simul$data, centers = id)
  myclusters <- mykmeans$cluster
}
if (algo == "gmm") {
  myclust <- mclust::Mclust(data = simul$data, G = id, verbose = FALSE)
  myclusters <- myclust$classification
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# GAP statistic
tmptime <- system.time({
  out <- GapStatistic(
    xdata = simul$data,
    method = algo
  )
})
gap <- out$gap
id <- ManualArgmaxId(gap)
if (algo == "kmeans") {
  set.seed(1)
  mykmeans <- stats::kmeans(x = simul$data, centers = id)
  myclusters <- mykmeans$cluster
}
if (algo == "gmm") {
  myclust <- mclust::Mclust(data = simul$data, G = id, verbose = FALSE)
  myclusters <- myclust$classification
}
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# BIC
if (algo == "gmm") {
  bic <- InformationTheoryCalibration(xdata = simul$data, index = "bic")
  id <- ManualArgmaxId(bic)
  myclust <- mclust::Mclust(data = simul$data, G = id, verbose = FALSE)
  myclusters <- myclust$classification
  nperf <- rbind(
    nperf,
    c(
      G = id,
      ClusteringPerformance(theta = myclusters, theta_star = simul),
      signif = NA,
      time = as.numeric(tmptime[1])
    )
  )
}

# Consensus clustering
if (algo == "kmeans") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = KMeansClustering,
      linkage = linkage,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
    )
  })
}
if (algo == "gmm") {
  tmptime <- system.time({
    stab <- Clustering(
      xdata = simul$data,
      implementation = GMMClustering,
      linkage = linkage,
      K = K,
      verbose = FALSE,
      nc = 1:nc_max,
    )
  })
}

# Consensus clustering with G*
nperf <- rbind(
  nperf,
  c(
    G = length(n),
    ClusteringPerformance(theta = Clusters(stab, argmax_id = length(n)), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
    ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

if (algo == "kmeans") {
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
      ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
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
      ClusteringPerformance(theta = Clusters(stab, argmax_id = id), theta_star = simul),
      signif = ifelse(scores$NORM_P[id] < 0.05, yes = 1, no = 0),
      time = as.numeric(tmptime3[1] + (K - 2) / K * tmptime[1])
    )
  )
}

# Consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Re-formatting output object
if (algo == "kmeans") {
  rownames(nperf) <- c(
    "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap",
    "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db", "consensus_gap",
    "delta", "pac", "pins_discrepancy", "rcsi_pac", "rcsi_entropy", "consensus_score"
  )
}
if (algo == "gmm") {
  rownames(nperf) <- c(
    "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap", "single_run_bic",
    "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db", "consensus_gap",
    "delta", "pac", "pins_discrepancy", "consensus_score"
  )
}

# Saving output object
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
