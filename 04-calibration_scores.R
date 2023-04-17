rm(list = ls())

library(fake)
library(sharp)
library(colorspace)
library(diceR)
library(PINSPlus)
library(doSNOW)
library(foreach)
library(ggplot2)
library(M3C)
library(abind)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
n_tot <- 150
p <- 10
n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
pk <- round(rep(0.2, 5) * p)
ev_xc <- 0.6
v_min <- v_max <- 0
nu_xc <- 1
set.seed(1)
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
x <- scale(simul$data)
par(mfrow = c(1, 1), mar = c(5, 5, 1, 6))
Heatmap(as.matrix(dist(x)))

# Hierarchical clustering
myhclust <- hclust(d = dist(x), method = "complete")
myclusters <- cutree(myhclust, k = 3)
ClusteringPerformance(theta = myclusters, theta_star = simul)

# Stability
stab <- Clustering(xdata = x, implementation = HierarchicalClustering, nc = 1:20)

# Delta score
delta <- DeltaAreaCDF(stab)

# PAC score
pac <- PAC(stab)

# PINS discrepancy score
discrepancy <- PINSDiscrepancy(x = simul$data, stab)

# M3C score (PAC)
scores <- MonteCarloScore(x = simul$data, stab, objective = "PAC")
rcsi_pac <- scores$RCSI # criterion to define assignments in their code

# M3C score (entropy)
scores <- MonteCarloScore(x = simul$data, stab, objective = "entropy")
rcsi_entropy <- scores$RCSI # criterion to define assignments in their code

# Clustering performances for different numbers of clusters
perf <- AllPerf(stab)

# Figures
{
  pdf("Figures/Score_vs_performance_simul.pdf",
    width = 6, height = 6
  )
  par(mar = rep(5, 4))
  Heatmap(as.matrix(dist(x)))
  dev.off()
}

{
  pdf("Figures/Score_vs_performance.pdf",
    width = 11, height = 7
  )
  par(mfrow = c(2, 3))
  par(mar = c(5, 5, 0.5, 0.5))
  ScatterPerf(x = delta, perf = perf, xlab = expression(Delta))
  ScatterPerf(x = -pac, perf = perf, xlab = "- PAC")
  ScatterPerf(x = discrepancy, perf = perf, xlab = "PINS discrepancy")
  ScatterPerf(x = rcsi_pac, perf = perf, xlab = "RCSI (PAC)")
  ScatterPerf(x = rcsi_entropy, perf = perf, xlab = "RCSI (entropy)")
  ScatterPerf(x = stab$Sc, perf = perf, xlab = "Consensus score")
  dev.off()
}
