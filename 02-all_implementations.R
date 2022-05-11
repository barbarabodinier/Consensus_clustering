rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

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

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for (name in r) eval(parse(text = paste0(name, "<-sharp:::", name)))

# Loading all additional functions
myfunctions <- list.files("Scripts/Functions/")
myfunctions <- myfunctions[myfunctions != "Former"]
for (k in 1:length(myfunctions)) {
  source(paste0("Scripts/Functions/", myfunctions[k]))
}

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
set.seed(1)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 10, nu_xc = 1,
  ev_xc = 0.7
)
x <- simul$data


# Hierarchical clustering
myhclust <- hclust(d = dist(x), method = "complete")
myclusters <- cutree(myhclust, k = 3)
ClusteringPerformance(theta = myclusters, theta_star = simul)

# Stability
stab <- Clustering(xdata = x, implementation = HierarchicalClustering)

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
scores <- MonteCarloScore(stab, objective = "PAC")
rcsi_pac <- scores$RCSI # criterion to define assignments in their code
plot(rcsi_pac)
which.max(rcsi_pac)

# M3C score (entropy)
scores <- MonteCarloScore(stab, objective = "entropy")
rcsi_entropy <- scores$RCSI # criterion to define assignments in their code
plot(rcsi_entropy)
which.max(rcsi_entropy)

# Clustering performances for different numbers of clusters
perf <- AllPerf(stab)

par(mfrow = c(2, 3))
plot(delta, perf$ari)
plot(-pac, perf$ari)
plot.new()
plot(rcsi_pac, perf$ari)
plot(rcsi_entropy, perf$ari)
plot(stab$Sc, perf$ari)







# M3C with objective entropy from their code
objective <- "PAC"
objective <- "entropy"

out <- M3C(mydata = t(x), clusteralg = "hc", maxK = 10, pItem = 0.5, repsreal = 100, seed = 1, objective = objective)

real <- out$scores[, 1:2]
colnames(real)[2] <- "PAC_REAL"
ls <- out$refpacscores
real$PAC_REF <- colMeans(ls)

## if PAC/entropy is zero set it to really really small
ptemp <- real$PAC_REAL
ptemp[ptemp == 0] <- 0.0001 ## changed
pacreal <- ptemp

diffM <- sweep(log(ls), 2, log(pacreal))
real$RCSI <- colMeans(diffM)
real$RCSI_SE <- (apply(diffM, 2, sd)) / sqrt(nrow(ls))

## usual p value derivation
pvals <- vapply(seq_len(ncol(ls)), function(i) {
  distribution <- as.numeric(ls[, i])
  ((length(distribution[distribution < real$PAC_REAL[i]])) + 1) / (iters + 1) # (b+1)/(m+1)=pval
}, numeric(1))
real$MONTECARLO_P <- pvals

if (objective == "PAC") {
  ## estimate p values using a beta distribution
  variance <- apply(ls, 2, var)
  pvals2 <- vapply(seq_len(nrow(real)), function(i) {
    mean <- real$PAC_REF[i]
    var <- variance[[i]]
    realpac <- real$PAC_REAL[i]
    params2 <- M3C:::estBetaParams(mu = mean, var = var)
    pbeta(realpac, params2[[1]], params2[[2]])
  }, numeric(1))
  real$BETA_P <- pvals2
  real$P_SCORE <- -log10(real$BETA_P)
} else if (objective == "entropy") {
  ## estimate p values using a beta distribution
  variance <- apply(ls, 2, sd)
  pvals2 <- vapply(seq_len(nrow(real)), function(i) {
    mean <- real$PAC_REF[i]
    var <- variance[[i]]
    realpac <- real$PAC_REAL[i]
    pnorm(realpac, mean = mean, sd = var)
  }, numeric(1))
  real$NORM_P <- pvals2
  real$P_SCORE <- -log10(real$NORM_P)
  # fix names
  colnames(real)[2:3] <- c("ENTROPY_REAL", "ENTROPY_REF")
}

real

out$scores


# M3C with our consensus matrix
# objective="PAC"
objective <- "entropy"

MonteCarloScore(stab, objective = "PAC")

MonteCarloScore <- function(stab, K_ref = 25, objective = "entropy") {
  # Running M3C for reference distribution
  out <- M3C(
    mydata = t(x),
    iters = K_ref,
    clusteralg = "hc",
    maxK = max(stab$nc),
    pItem = stab$params$tau,
    repsreal = 2,
    seed = 1,
    objective = objective
  )
  # print(out$scores)
  #
  # # Storing consensus matrices
  # stab=NULL
  # stab$coprop=CoMembership(groups=rep(1, nrow(out$realdataresults[2][[1]]$consensus_matrix)))
  # for (i in 2:length(out$realdataresults)){
  #   stab$coprop=abind(stab$coprop, out$realdataresults[i][[1]]$consensus_matrix, along=3)
  # }
  # stab$nc=cbind(1:length(out$realdataresults))

  # Calculation of PAC scores
  if (objective == "PAC") {
    real <- data.frame(K = stab$nc, PAC_REAL = PAC(stab))
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Calculation of entropy scores
  if (objective == "entropy") {
    entropies <- rep(NA, dim(stab$coprop)[3])
    for (i in 2:dim(stab$coprop)[3]) {
      entropies[i] <- M3C:::entropy(stab$coprop[, , i])
    }
    real <- data.frame(K = stab$nc, PAC_REAL = entropies)
    real <- real[-1, ]
    rownames(real) <- 1:nrow(real)
  }

  # Storing the reference
  ls <- out$refpacscores

  # Calculating reference mean
  real$PAC_REF <- colMeans(ls)

  # Checking PAC values
  ptemp <- real$PAC_REAL
  ptemp[ptemp == 0] <- 0.0001
  pacreal <- ptemp

  # Calculating RCSI score
  diffM <- sweep(log(ls), 2, log(pacreal))
  real$RCSI <- colMeans(diffM)
  real$RCSI_SE <- (apply(diffM, 2, sd)) / sqrt(nrow(ls))

  # Calculating p-value
  pvals <- vapply(seq_len(ncol(ls)), function(i) {
    distribution <- as.numeric(ls[, i])
    ((length(distribution[distribution < real$PAC_REAL[i]])) + 1) / (iters + 1) # (b+1)/(m+1)=pval
  }, numeric(1))
  real$MONTECARLO_P <- pvals

  if (objective == "PAC") {
    variance <- apply(ls, 2, var)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      params2 <- M3C:::estBetaParams(mu = mean, var = var)
      pbeta(realpac, params2[[1]], params2[[2]])
    }, numeric(1))
    real$BETA_P <- pvals2
    real$P_SCORE <- -log10(real$BETA_P)
  } else if (objective == "entropy") {
    variance <- apply(ls, 2, sd)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      pnorm(realpac, mean = mean, sd = var)
    }, numeric(1))
    real$NORM_P <- pvals2
    real$P_SCORE <- -log10(real$NORM_P)
    # fix names
    colnames(real)[2:3] <- c("ENTROPY_REAL", "ENTROPY_REF")
  }

  return(real)
}

out$scores










real <- data.frame(PAC_REAL = out$scores$PAC_REAL)
ls <- out$refpacscores
iters <- nrow(ls)
pvals <- vapply(seq_len(ncol(ls)), function(i) {
  distribution <- as.numeric(ls[, i])
  ((length(distribution[distribution < real$PAC_REAL[i]])) + 1) / (iters + 1) # (b+1)/(m+1)=pval
}, numeric(1))
real$MONTECARLO_P <- pvals

real$PAC_REF <- colMeans(ls)
variance <- apply(ls, 2, sd)
pvals2 <- vapply(seq_len(nrow(real)), function(i) {
  mean <- real$PAC_REF[i]
  var <- variance[[i]]
  realpac <- real$PAC_REAL[i]
  pnorm(realpac, mean = mean, sd = var)
}, numeric(1))
real$NORM_P <- pvals2
real$P_SCORE <- -log10(real$NORM_P)
# fix names
colnames(real)[2:3] <- c("ENTROPY_REAL", "ENTROPY_REF")

M3C(mydata = t(x), clusteralg = "hc", maxK = 10, pItem = 0.5, repsreal = 100, seed = 1, objective = "entropy")$scores



plot(PAC, stab$Sc[2:10],
  xlab = "PAC", ylab = "Stability score"
)
text(PAC, stab$Sc[2:10], labels = Kvec, pos = 4)










areas <- rep(NA, dim(stab$coprop)[3])
for (k in 2:dim(stab$coprop)[3]) {
  areas[k] <- AreaUnderCDF(M = stab$coprop[, , k])$area
}
names(areas) <- stab$nc[, 1]
plot(areas)
plot(DeltaArea(areas[-1]))




# diceR
# https://cran.r-project.org/web/packages/diceR/vignettes/overview.html
out <- consensus_cluster(data = x, algorithms = "hc", reps = 100, p.item = 0.5, nk = 2:10)

# Consensus Cluster Plus
out <- ConsensusClusterPlus(
  d = t(x), maxK = 10, reps = 100, pItem = 0.5,
  distance = "euclidean", seed = 1,
  writeTable = TRUE,
  innerLinkage = "complete", finalLinkage = "complete"
)
out[[3]] # visual in plot
out[[3]]$consensusClass
ClusteringPerformance(theta = out[[3]]$consensusClass, theta_star = simul$theta)
calc <- calcICL(out)

results <- out
maxK <- length(results)
Kvec <- 2:maxK
x1 <- 0.1
x2 <- 0.9 # threshold defining the intermediate sub-interval
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
for (i in Kvec) {
  M <- results[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
} # end for i
# The optimal K
optK <- Kvec[which.min(PAC)]


results <- stab
maxK <- 10
Kvec <- 2:maxK
x1 <- 0.1
x2 <- 0.9 # threshold defining the intermediate sub-interval
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
for (i in Kvec) {
  M <- results$coprop[, , i]
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
} # end for i
# The optimal K
optK <- Kvec[which.min(PAC)]

plot(PAC, stab$Sc[2:10],
  xlab = "PAC", ylab = "Stability score"
)
text(PAC, stab$Sc[2:10], labels = Kvec, pos = 4)

par(mfrow = c(1, 2))
Heatmap(ConsensusMatrix(stab, argmax_id = matrix(c(2, 1), nrow = 1)))
Heatmap(ConsensusMatrix(stab, argmax_id = matrix(c(3, 1), nrow = 1)))

hat_k <- which.max(PAC) + 1
hat_k <- 3
myclusters <- cutree(hclust(as.dist(1 - ConsensusMatrix(stab, argmax_id = matrix(c(hat_k, 1), nrow = 1)))), k = hat_k)
ClusteringPerformance(theta = myclusters, theta_star = simul$theta, implementation = igraph::cluster_walktrap)

# M3C
out <- M3C(mydata = t(x), clusteralg = "hc", maxK = 10, pItem = 0.5, repsreal = 100, seed = 1)
out$scores
out$assignments
ClusteringPerformance(theta = out$assignments, theta_star = simul$theta)
ClusteringPerformance(theta = out$realdataresults[[3]]$assignments, theta_star = simul$theta)

m3c_cm <- out$realdataresults[[3]]$consensus_matrix
rownames(m3c_cm) <- colnames(m3c_cm) <- rownames(out$realdataresults[[3]]$ordered_annotation)
m3c_cm <- m3c_cm[rownames(ConsensusMatrix(stab)), rownames(ConsensusMatrix(stab))]
Heatmap(m3c_cm)
Heatmap(ConsensusMatrix(stab))


# Stability
stab <- Clustering(xdata = x, implementation = HierarchicalClustering)
ClusteringPerformance(theta = stab, theta_star = simul$theta, implementation = igraph::cluster_walktrap)
ClusteringPerformance(theta = stab, theta_star = simul$theta, implementation = igraph::cluster_louvain)
ClusteringPerformance(theta = stab, theta_star = simul$theta, implementation = igraph::components)

M <- ConsensusMatrix(stability = stab)
myhclust <- hclust(as.dist(1 - M))
cutree(myhclust, h = Argmax(stab, clustering = TRUE)[2])
ClusteringPerformance(theta = cutree(myhclust, h = Argmax(stab, clustering = TRUE)[2]), theta_star = simul$theta)

ClusteringPerformance(theta = cutree(myhclust, k = Argmax(stab, clustering = TRUE)[1]), theta_star = simul$theta)

Graph(Adjacency(stab))
