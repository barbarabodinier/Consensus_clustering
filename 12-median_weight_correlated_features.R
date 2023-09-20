rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(igraph)
library(colorspace)
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
n_tot <- 300
p <- 100
nu_xc <- 0.3
ev_xc <- 0.7
v_min <- 1
v_max <- 1

set.seed(0)
n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
pk <- round(rep(0.2, 5) * p)
q <- round(nu_xc * p)
theta_xc <- sample(c(rep(1, q), rep(0, p - q)))
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
Heatmap(cor(simul$data[which(simul$theta == 2), ]),
  col = c("navy", "white", "darkred"),
  legend_range = c(-1, 1)
)

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
  pdf("Figures/Distribution_median_weight_correlated.pdf",
    width = 15, height = 8
  )
  layout(
    mat = matrix(c(1, 2, 3, 4, 4, 4), byrow = TRUE, nrow = 2),
  )
  par(mar = c(5, 5, 1, 7))
  Heatmap(as.matrix(dist(simul$data)))
  Heatmap(cor(simul$data[which(simul$theta == 2), ]),
    col = c("navy", "white", "darkred"),
    legend_range = c(-1, 1)
  )
  CalibrationPlot(stab, xlab = "Number of clusters", ylab = "Consensus score")
  WeightBoxplot(stab,
    col = ifelse(theta_xc == 1, yes = "red", no = "grey"),
    frame = TRUE, boxwex = 0.2,
    ylab = "Median weight",
  )
  for (i in 1:ncol(stab$Beta)) {
    axis(
      side = 1, at = i, las = 2,
      labels = colnames(stab$Beta)[i],
      col.axis = ifelse(theta_xc[i] == 1, yes = "red", no = "grey")
    )
  }
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
