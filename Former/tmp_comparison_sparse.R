rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

library(FactoMineR)
library(sharp)
library(aricode)
library(M3C)
library(abind)
library(colorspace)

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

# Simulation study parameters
simul_study_id <- 1
params_id <- 1
seed <- 1
simulation_id <- paste0(params_id, "_", seed)

# Extracting simulation parameters
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Model parameters
K <- 100
iters <- 25 # default 25, recommended 5-100

# Extracting simulation parameters
nc <- params_list[params_id, "nc"]
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


n=300
ev_xc=0.9
p=50
nu_xc=0.25
equal_size=TRUE



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

# Principal Component Analysis
mycolours=darken(c("navy", "red", "forestgreen", "orange", "purple"), amount=0.2)
par(mfrow=c(1,3), mar = c(5, 4.5, 2, 4.5))
mypca <- PCA(simul$data, graph = FALSE)
tmp <- matrix(c(1, 2, 1, 3, 2, 3), byrow = TRUE, ncol = 2)
for (k in 1:nrow(tmp)) {
  xcomp <- tmp[k, 1]
  ycomp <- tmp[k, 2]
  plot(NA,
       xlim = range(mypca$ind$coord[, 1:max(tmp)]),
       ylim = range(mypca$ind$coord[, 1:max(tmp)]),
       las = 1, cex.lab = 1.5, cex.lab = 1.5,
       xlab = paste0("Comp ", xcomp, " (", round(mypca$eig[xcomp, 2], digits = 2), "% e.v.)"),
       ylab = paste0("Comp ", ycomp, " (", round(mypca$eig[ycomp, 2], digits = 2), "% e.v.)")
  )
  abline(v = axTicks(1), lty = 3, col = "grey")
  abline(h = axTicks(2), lty = 3, col = "grey")
  text(mypca$ind$coord[, xcomp], mypca$ind$coord[, ycomp],
       pch = 19, cex = 0.7, las = 1,
       labels = gsub("obs", "", rownames(mypca$ind$coord)),
       col = mycolours[simul$theta]
  )
}

# Hierarchical clustering with G*
tmptime <- system.time({
  myhclust <- hclust(d = dist(simul$data), method = "complete")
  myclusters <- cutree(myhclust, k = length(n))
})
nperf <- c(
  G = length(n),
  ClusteringPerformance(theta = myclusters, theta_star = simul),
  signif = NA,
  time = as.numeric(tmptime[1])
)

# Consensus clustering
tmptime <- system.time({
  stab <- Clustering(
    xdata = simul$data,
    implementation = HierarchicalClustering,
    K = K,
    verbose = FALSE,
    nc = 1:20,
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

# Consensus clustering with consensus score
nperf <- rbind(
  nperf,
  c(
    G = length(n),
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

# Sparse consensus clustering
tmptime <- system.time({
  stab <- Clustering(
    xdata = simul$data,
    implementation = SparseHierarchicalClustering,
    K = K,
    verbose = FALSE,
    nc = 1:20,
    Lambda=LambdaSequence(lmax=10, lmin=1+1e-5, cardinal=50)
  )
})

# Sparse consensus clustering with consensus score
nperf <- rbind(
  nperf,
  c(
    G = Argmax(stab)[1],
    ClusteringPerformance(theta = stab, theta_star = simul),
    signif = NA,
    time = as.numeric(tmptime[1])
  )
)

Argmax(stab)

par(mfrow=c(1,1))
plot(stab$Sc, col=ifelse(stab$nc==5, yes="red", no="black"))

i=1
plot(stab$Sc[(i-1)*20+seq(1, 20)], type="l",
     ylim=c(0,1))
for (i in 2:20){
  lines(stab$Sc[(i-1)*20+seq(1, 20)])
}

plot(apply(stab$Beta[405,,],1,median), stab$selprop[405,], col=ifelse(simul$theta_xc==1, yes="red", no="black"))
abline(h=Argmax(stab)[1, "pi"], lty=2, col="darkred")







mybetas=stab_sparse$Beta[ArgmaxId(stab_sparse)[1],,]
j=1
plot(mybetas[j,], type="l", ylim=c(0,0.4), 
     col=ifelse(simul$theta_xc[j]==1, yes="red", no="black"))
for (j in 1:nrow(mybetas)){
  lines(mybetas[j,], type="l",
        col=ifelse(simul$theta_xc[j]==1, yes="red", no="black"))
}


i=1
plot(stab$nc[(i-1)*20+seq(1, 20)], stab$Sc[(i-1)*20+seq(1, 20)], type="l",
     ylim=c(0,1))
for (i in 2:20){
  lines(stab$nc[(i-1)*20+seq(1, 20)], stab$Sc[(i-1)*20+seq(1, 20)],)
}
lines(stab$nc[which(stab$S_2d==max(stab$S_2d, na.rm=TRUE), arr.ind = TRUE)[,1]],
      stab$Sc[which(stab$S_2d==max(stab$S_2d, na.rm=TRUE), arr.ind = TRUE)[,1]], 
      col="red")

i=1
plot(stab$Q[seq(i,length(stab$Sc),by=30)],
     stab$Sc[seq(i,length(stab$Sc),by=30)], type="l",
     ylim=c(0,1))
for (i in 2:30){
  lines(stab$Q[seq(i,length(stab$Sc),by=30)], 
        stab$Sc[seq(i,length(stab$Sc),by=30)])
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
rownames(nperf) <- c("hclust_star", "consensus_star", "delta", "pac", "rcsi_pac", "rcsi_entropy", "consensus_score")

# Saving output object
saveRDS(nperf, paste0(filepath, "Performances_", simulation_id, ".rds"))

# Saving Spearman's correlation
spearman <- c(
  rcsi_pac = cor(stab$Sc, rcsi_pac, method = "spearman", use = "complete.obs"),
  rcsi_entropy = cor(stab$Sc, rcsi_entropy, method = "spearman", use = "complete.obs")
)
saveRDS(spearman, paste0(filepath, "Correlations_", simulation_id, ".rds"))
