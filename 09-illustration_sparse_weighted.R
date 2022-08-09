rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

library(fake)
library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(FactoMineR)
library(diceR)
library(ConsensusClusterPlus)
library(M3C)
library(abind)
library(rCOSA)
library(plotrix)
library(openxlsx)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Model parameters
K <- 100
nc_max <- 20
noit <- 20
niter <- 10

# Simulation of data with clusters
simul_study_id <- 1
print(paste0("Simulation study ", simul_study_id))
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
params_id <- 1
n_tot <- params_list[params_id, "n_tot"]
v_min <- params_list[params_id, "v_min"]
v_max <- params_list[params_id, "v_max"]

# Manually set for sparse/weighted
p <- 100
nu_xc <- 0.2
ev_xc <- 0.8

set.seed(0)
n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
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
simul$data <- scale(simul$data)

# Heatmap of pairwise distances
Heatmap(as.matrix(dist(simul$data)))

# Consensus clustering (unweighted)
tmptime <- system.time({
  stab_unweighted <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    nc = 1:nc_max
  )
})
print(tmptime)

# Consensus clustering (sparcl)
tmptime <- system.time({
  stab_sparcl <- Clustering(
    xdata = simul$data,
    implementation = SparseHierarchicalClustering,
    K = K,
    nc = 1:nc_max,
    Lambda = LambdaSequence(lmax = 10, lmin = 1.1, cardinal = 10)
  )
})
print(tmptime)

# Consensus clustering (cosa)
tmptime <- system.time({
  stab_cosa <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    nc = 1:nc_max,
    Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
    noit = noit,
    niter = niter,
    verbose = TRUE
  )
})
print(tmptime)

# Making figure
{
  pdf("Figures/Example_sparse_weighted.pdf",
    height = 11, width = 16
  )
  layout(
    mat = matrix(c(1:6), byrow = TRUE, ncol = 2),
    widths = c(1, 3)
  )
  par(mar = c(5, 5, 1, 1))
  mycolours <- ifelse(apply(simul$theta_xc, 2, sum) > 0, yes = "red", no = "grey")
  CalibrationPlot(stab_unweighted, ylab = "Unweighted", legend = FALSE)
  plot.new()
  CalibrationPlot(stab_sparcl,
    ylab = "Sparse (sparcl)",
    ncol = 2, cex.legend = 1
  )
  WeightBoxplot(stab_sparcl, col = mycolours, frame = "T")
  for (i in 1:ncol(stab_cosa$Beta)) {
    axis(
      side = 1, at = i, las = 2,
      labels = colnames(stab_cosa$Beta)[i],
      col.axis = mycolours[i],
      font = ifelse(mycolours[i] == "red", yes = 2, no = 1)
    )
  }
  for (i in 1:ncol(stab_sparcl$Beta)) {
    axis(
      side = 3, at = i, las = 2,
      labels = formatC(stab_sparcl$selprop[ArgmaxId(stab_sparcl)[1], i],
        format = "f", digits = 2
      ),
      col.axis = ifelse(simul$theta_xc[i] == 1, yes = "red", no = "black"),
      font = ifelse(simul$theta_xc[i] == 1, yes = 2, no = 1)
    )
  }
  CalibrationPlot(stab_cosa,
    xlab = "Number of clusters",
    ylab = "Weighted (COSA)",
    ncol = 2, cex.legend = 1
  )
  WeightBoxplot(stab_cosa, col = mycolours, frame = "T")
  for (i in 1:ncol(stab_cosa$Beta)) {
    axis(
      side = 1, at = i, las = 2,
      labels = colnames(stab_cosa$Beta)[i],
      col.axis = ifelse(simul$theta_xc[i] == 1, yes = "red", no = "black"),
      font = ifelse(simul$theta_xc[i] == 1, yes = 2, no = 1)
    )
  }
  dev.off()
}

# Clustering performances
continuous_metrics <- c("rand", "ari", "jaccard")
mytable <- rbind(
  c(
    Argmax(stab_unweighted)[1],
    "",
    formatC(as.numeric(ClusteringPerformance(stab_unweighted, simul)[, continuous_metrics]),
      format = "f", digits = 2
    )
  ),
  c(
    Argmax(stab_sparcl)[1],
    formatC(as.numeric(Argmax(stab_sparcl)[2]), format = "f", digits = 2),
    formatC(as.numeric(ClusteringPerformance(stab_sparcl, simul)[, continuous_metrics]),
      format = "f", digits = 2
    )
  ),
  c(
    Argmax(stab_cosa)[1],
    formatC(as.numeric(Argmax(stab_cosa)[2]), format = "f", digits = 2),
    formatC(as.numeric(ClusteringPerformance(stab_cosa, simul)[, continuous_metrics]),
      format = "f", digits = 2
    )
  )
)
colnames(mytable) <- c("G", "lambda", continuous_metrics)
rownames(mytable) <- c("unweighted", "sparcl", "cosa")

write.xlsx(as.data.frame(mytable),
  "Tables/Example_sparse_weighted.xlsx",
  overwrite = TRUE,
  colNames = TRUE, rowNames = TRUE
)
