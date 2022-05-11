library(sharp)
library(aricode)
library(M3C)
library(abind)
library(cluster)

setwd("../../")

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for (name in r) eval(parse(text = paste0(name, "<-sharp:::", name)))

# Loading all additional functions
myfunctions <- list.files("Scripts/Functions/")
for (k in 1:length(myfunctions)) {
  source(paste0("Scripts/Functions/", myfunctions[k]))
}

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

# Creating folder of simulation study
dir.create("Results", showWarnings = FALSE)
dir.create("Results/Simulations_consensus_hierarchical", showWarnings = FALSE)
filepath <- paste0("Results/Simulations_consensus_hierarchical/Simulations_", simul_study_id, "/")
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
  nu_xc = nu_xc,
  output_matrices = TRUE
)

# Hierarchical clustering with G*
tmptime <- system.time({
  mydist <- dist(simul$data)
  myhclust <- hclust(d = mydist, method = "complete")
  myclusters <- cutree(myhclust, k = length(n))
})
nperf <- c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Hierarchical clustering with max silhouette score
silhouette <- SilhouetteScore(mydist, myhclust)
id <- ManualArgmaxId(silhouette)
myclusters <- cutree(myhclust, k = id)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Hierarchical clustering with max GAP statistic
tmptime <- system.time({
  out <- GapStatistic(xdata = simul$data)
})
gap <- out$gap
id <- ManualArgmaxId(gap)
myclusters <- cutree(myhclust, k = id)
nperf <- rbind(
  nperf,
  c(
    G = id,
    ClusteringPerformance(theta = myclusters, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Consensus clustering
tmptime <- system.time({
  stab <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    verbose = FALSE,
    nc = 1:nc_max,
  )
})

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

# M3C score (PAC)
tmptime2 <- system.time({
  scores <- MonteCarloScore(x = simul$data, stab, objective = "PAC", iters = iters)
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
  scores <- MonteCarloScore(x = simul$data, stab, objective = "entropy", iters = iters)
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
rownames(nperf) <- c(
  "hclust_star", "hclust_silhouette", "hclust_gap",
  "consensus_star", "delta", "pac", "rcsi_pac", "rcsi_entropy", "consensus_score"
)

# Saving output object
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))

# Saving Spearman's correlation
spearman <- c(
  rcsi_pac = cor(stab$Sc, rcsi_pac, method = "spearman", use = "complete.obs"),
  rcsi_entropy = cor(stab$Sc, rcsi_entropy, method = "spearman", use = "complete.obs")
)
saveRDS(spearman, paste0(filepath, "Correlations_", simulation_id, ".rds"))
