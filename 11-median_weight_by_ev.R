rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(igraph)
library(colorspace)
library(diceR)
library(M3C)
library(abind)
library(plotrix)

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
nu_xc <- 0.3

set.seed(1)
n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
pk <- round(rep(0.2, 5) * p)
q <- round(nu_xc * p)
theta_xc <- c(rep(1, q), rep(0, p - q))
ev_xc <- c(seq(0.2, 0.99, length.out = q), rep(0, p - q))
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

# Heatmap of pairwise distances
Heatmap(as.matrix(dist(simul$data)))

stab <- Clustering(
  xdata = simul$data,
  implementation = HierarchicalClustering,
  K = K,
  nc = 1:nc_max,
  Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
  noit = noit,
  niter = niter,
  verbose = TRUE
)

{
  pdf("Figures/Distribution_median_weight_by_ev.pdf",
    width = 18, height = 6
  )
  layout(
    mat = matrix(c(1, 2), byrow = TRUE, nrow = 1),
    widths = c(1, 2)
  )
  par(mar = c(5, 5, 1, 1))
  CalibrationPlot(stab, xlab = "Number of clusters", ylab = "sharp score")
  WeightBoxplot(stab,
    at = simul$ev, col = ifelse(theta_xc == 1, yes = "red", no = "grey"),
    frame = TRUE, boxwex = 0.007,
    xlab = "Proportion of explained variance",
    ylab = "Median weight"
  )
  mycolours <- c("grey", "grey30", "grey30", "grey")
  quantile_list <- c(0, 0.25, 0.75, 1)
  for (i in 1:length(quantile_list)) {
    myquantile <- quantile_list[i]
    z <- apply(stab$Beta[ArgmaxId(stab)[1], , ], 1, quantile, probs = myquantile)
    if (myquantile > 0.5) {
      abline(h = max(z[which(simul$ev == 0)]), lty = 3, col = mycolours[i])
    } else {
      abline(h = min(z[which(simul$ev == 0)]), lty = 3, col = mycolours[i])
    }
  }
  dev.off()
}
