library(fake)
library(rCOSA)
library(sharp)
library(abind)
library(sparcl)
library(cluster)
library(clustvarsel)
library(VarSelLCM)

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
dir.create(paste0("Results/Simulations_consensus_sparse_", algo), showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_sparse_", algo, "/Simulations_", simul_study_id, "/")
dir.create(filepath, showWarnings = FALSE)
print("Path to results:")
print(filepath)

# Data simulation
set.seed(seed)
if (equal_size) {
  n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
} else {
  if (nc==5) {
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

# Clustering with G*
tmptime <- system.time({
  myclust <- mclust::Mclust(data = simul$data, G = length(n), verbose = FALSE)
  myclusters <- myclust$classification
})
nperf <- c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Clustering calibrated with BIC
tmptime <- system.time({
  out <- GMM0Calibration(xdata = simul$data, nc_max = nc_max, index = "bic")
})
id <- ManualArgmaxId(out$scores)
nperf <- rbind(nperf, c(
  G = id,
  ClusteringPerformance(theta = out$clusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
))

# Clustvarsel with G*
tmptime <- system.time({
  myclustvarsel <- clustvarsel(simul$data,
    search = "headlong",
    G = length(n),
    verbose = FALSE
  )
  myclusters <- myclustvarsel$model$classification
})

# Clustering performance with G*
nperf <- rbind(nperf, c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
))

# Selection performance with G*
selected <- rep(0, ncol(simul$data))
names(selected) <- colnames(simul$data)
selected[names(myclustvarsel$subset)] <- 1
selperf <- NULL
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = selected,
    theta_star = theta_xc
  )
)

# Clustvarsel calibrated by BIC
tmptime <- system.time({
  out <- GMM1Calibration(
    xdata = simul$data,
    nc_max = nc_max,
    search = "headlong",
    index = "bic"
  )
})

# Clustering performance when calibrated by BIC
id <- ManualArgmaxId(out$scores)
nperf <- rbind(nperf, c(
  G = id,
  ClusteringPerformance(theta = out$clusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
))

# Selection performance when calibrated by BIC
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = out$selected,
    theta_star = theta_xc
  )
)

# VarSelLCM with G*
tmptime <- system.time({
  myvarsellcm <- VarSelCluster(
    x = simul$data,
    gvals = length(n),
    vbleSelec = TRUE,
    crit.varsel = "MICL"
  )
  myclusters <- fitted(myvarsellcm)
})

# Clustering performance with G*
nperf <- rbind(nperf, c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
))

# Selection performance with G*
selected <- rep(0, ncol(simul$data))
names(selected) <- colnames(simul$data)
selected[myvarsellcm@model@names.relevant] <- 1
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = selected,
    theta_star = theta_xc
  )
)

# VarSelLCM calibrated by MICL
tmptime <- system.time({
  out <- GMM2Calibration(
    xdata = simul$data,
    nc_max = nc_max,
    index = "micl"
  )
})

# Clustering performance when calibrated by MICL
id <- ManualArgmaxId(out$scores)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = out$clusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Selection performance when calibrated by MICL
selperf <- rbind(
  selperf, SelectionPerformance(
    theta = out$selected,
    theta_star = theta_xc
  )
)

# Re-formatting output objects
rownames(nperf) <- c(
  "single_run_star",
  "single_run_bic",
  "clustvarsel_star",
  "clustvarsel_bic",
  "varsellcm_star",
  "varsellcm_micl"
)
rownames(selperf) <- c(
  "clustvarsel_star",
  "clustvarsel_bic",
  "varsellcm_star",
  "varsellcm_micl"
)

# Saving output objects
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))
saveRDS(selperf, paste0(filepath, "Selection_performances_", simulation_id, ".rds"))
