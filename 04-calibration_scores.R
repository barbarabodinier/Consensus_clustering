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

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
set.seed(0)
n <- c(20, 50, 30, 10, 40)
simul <- SimulateClustering(
  n = n,
  pk = 20,
  nu_xc = 1,
  ev_xc = 0.5
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

# Consensus score
plot(stab$Sc)
which.max(stab$Sc)

# Delta score
delta <- DeltaAreaCDF(stab)
plot(delta)
which.max(delta)

# PAC score
pac <- PAC(stab)
plot(pac)
which.min(pac)

# M3C score (PAC)
scores <- MonteCarloScore(x = simul$data, stab, objective = "PAC")
rcsi_pac <- scores$RCSI # criterion to define assignments in their code
plot(rcsi_pac)
which.max(rcsi_pac)

# M3C score (entropy)
scores <- MonteCarloScore(x = simul$data, stab, objective = "entropy")
rcsi_entropy <- scores$RCSI # criterion to define assignments in their code
plot(rcsi_entropy)
which.max(rcsi_entropy)

# Clustering performances for different numbers of clusters
perf <- AllPerf(stab)

# Figure
{
  pdf("Figures/Score_vs_performance.pdf",
    width = 12, height = 7
  )
  par(mfrow = c(2, 3))
  par(mar = c(5, 5, 0.5, 4))
  Heatmap(as.matrix(dist(x)))
  par(mar = c(5, 5, 0.5, 0.5))
  ScatterPerf(x = delta, perf = perf, xlab = expression(Delta))
  ScatterPerf(x = -pac, perf = perf, xlab = "- PAC")
  ScatterPerf(x = rcsi_pac, perf = perf, xlab = "RCSI (PAC)")
  ScatterPerf(x = rcsi_entropy, perf = perf, xlab = "RCSI (entropy)")
  ScatterPerf(x = stab$Sc, perf = perf, xlab = "Consensus score")
  dev.off()
}
