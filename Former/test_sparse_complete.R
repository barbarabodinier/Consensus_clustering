rm(list=ls())
setwd("~/Dropbox/Consensus_clustering")

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for(name in r) eval(parse(text=paste0(name, '<-sharp:::', name)))

# Loading all additional functions
myfunctions=list.files('Scripts/Functions/')
myfunctions=myfunctions[myfunctions!="Former"]
for (k in 1:length(myfunctions)){
  source(paste0('Scripts/Functions/', myfunctions[k]))
}

# # Simulation of data with clusters # XXX corresponding to a good example
# set.seed(1)
# n=rep(20,5)
# simul=SimulateClustering(n=n,
#                          pk=100, 
#                          nu_xc=0.2,
#                          ev_xc=0.9)
# x=simul$data

# Simulation of data with clusters
set.seed(1)
n=rep(30,5)
simul=SimulateClustering(n=n,
                         theta_xc = c(rep(1, 10), rep(0, 20)),
                         ev_xc=0.8)
x=simul$data

# Heatmap of pairwise distances
par(mfrow=c(1,2), mar=c(5,5,1,5))
thr=max(as.matrix(dist(simul$data)))
Heatmap(as.matrix(dist(simul$data)), 
        legend_range = c(0,thr))
thr=max(as.matrix(dist(simul$data[,which(simul$theta_xc==1)])))
Heatmap(as.matrix(dist(simul$data[,which(simul$theta_xc==1)])),
        legend_range = c(0,thr))

# Consensus clustering
stab=Clustering(xdata=simul$data, implementation = HierarchicalClustering)
plot(stab$Sc, ylim=c(0,1))
stab$nc[which.max(stab$Sc)]
Heatmap(ConsensusMatrix(stab))
ClusteringPerformance(stab, simul)

stab=Clustering(xdata=simul$data[,which(simul$theta_xc==1)], implementation = HierarchicalClustering)
plot(stab$Sc,ylim=c(0,1))

Heatmap(ConsensusMatrix(stab, argmax_id=which.max(stab$Sc)))
Heatmap(ConsensusMatrix(stab, argmax_id=5))

stab=Clustering(xdata=simul$data, 
                implementation = SparseHierarchicalClustering,
                Lambda=LambdaSequence(lmax=10, lmin=1.01, cardinal = 10))
plot(stab$Q)
plot(stab$Sc,ylim=c(0,1))
stab$nc[which.max(stab$Sc)]
stab$Lambda[which.max(stab$Sc)]
stab$Q[which.max(stab$Sc)]
Heatmap(ConsensusMatrix(stab, argmax_id=which.max(stab$Sc)))
ClusteringPerformance(stab, simul)
plot(stab$selprop[which.max(stab$Sc),], col=ifelse(simul$theta_xc==1, yes="red", no="grey30"))
abline(h=seq(0.6, 0.9, by = 0.01)[which.max(StabilityScore(stab$selprop[which.max(stab$Sc),], K=100, 
                         pi_list = seq(0.6, 0.9, by = 0.01)))], col="darkred", lty=2)
plot(apply(stab$Beta[which.max(stab$Sc),,], 1, mean),
  stab$selprop[which.max(stab$Sc),], col=ifelse(simul$theta_xc==1, yes="red", no="grey30"))

