rm(list=ls())
setwd("~/Dropbox/Consensus_clustering")

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(FactoMineR)
library(corpcor)
library(M3C)

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for(name in r) eval(parse(text=paste0(name, '<-sharp:::', name)))

# Loading all additional functions
myfunctions=list.files('Scripts/Functions/')
for (k in 1:length(myfunctions)){
  source(paste0('Scripts/Functions/', myfunctions[k]))
}

# Simulation of data with clusters
set.seed(1)
n=rep(20,5)
simul=SimulateClustering(n=n, pk=10, 
                         nu_xc = 1, ev_xc=0.5)

# Heatmap of pairwise distances
layout(mat=matrix(c(1,1,1,2,3,4), 
                  nrow=2, byrow=TRUE), 
       heights = c(3,1))
par(mar=c(5,5,5,1))
plot(simul)
x=simul$data

# Principal Component Analysis
par(mar=c(5,5,1,1))
mypca=PCA(simul$data, graph = FALSE)
tmp=matrix(c(1,2,1,3,2,3), byrow=TRUE, ncol=2)
for (k in 1:nrow(tmp)){
  xcomp=tmp[k,1]
  ycomp=tmp[k,2]
  plot(NA, 
       xlim=range(mypca$ind$coord[,1:max(tmp)]), 
       ylim=range(mypca$ind$coord[,1:max(tmp)]),
       xlab=paste0("Comp ",xcomp," (", round(mypca$eig[xcomp,2], digits=2), "% e.v.)"),
       ylab=paste0("Comp ",ycomp," (", round(mypca$eig[ycomp,2], digits=2), "% e.v.)"))
  abline(v=axTicks(1), lty=3, col="grey")
  abline(h=axTicks(2), lty=3, col="grey")
  text(mypca$ind$coord[,xcomp], mypca$ind$coord[,ycomp], pch=19, cex=0.5, las=1, 
       labels=gsub("obs", "", rownames(mypca$ind$coord)),
       col=simul$theta)
}

stab=Clustering(xdata=simul$data, implementation=HierarchicalClustering, tau=0.5, K=1000)

# pi_list=seq(0.51, 0.99, by = 0.01)
# pi_list=0.99
scores=matrix(NA, nrow=dim(stab$coprop)[3], ncol=1)
for (nc in 2:dim(stab$coprop)[3]){
  print(nc)
  scores[nc, ]=StabilityScoreClustering(selprop=ConsensusMatrix(stab, argmax_id = nc), nc=nc, pi_list=pi_list, K=1000)
  cat("\n")
}
rownames(scores)=1:dim(stab$coprop)[3]
colnames(scores)=pi_list
Heatmap(scores)
plot(apply(scores,1,max))
hat_nc=which.max(apply(scores,1,max))

Heatmap(ConsensusMatrix(stab, argmax_id=hat_nc))
Heatmap(ConsensusMatrix(stab, argmax_id=5))
ClusteringPerformance(theta=cutree(hclust(as.dist(1-ConsensusMatrix(stab, argmax_id=hat_nc)), method="complete"), k=hat_nc), 
                      theta_star = simul)
ClusteringPerformance(theta=cutree(hclust(as.dist(1-ConsensusMatrix(stab, argmax_id=hat_nc)), method="complete"), k=5), 
                      theta_star = simul)

ari=rep(NA, dim(stab$coprop)[3])
for (nc in 2:dim(stab$coprop)[3]){
  theta=cutree(hclust(as.dist(1-ConsensusMatrix(stab, argmax_id = nc))), k=nc)
  ari[nc]=ClusteringPerformance(theta=theta, theta_star=simul)$ari
}

par(mar=c(5,5,1,1))
plot(NA, las=1,
     xlab="Stability Ratio", 
     ylab="Adjusted Rand Index",
     cex.lab=1.5,
     xlim=range(scores, na.rm=TRUE), 
     ylim=range(ari, na.rm=TRUE))
text(scores, ari, labels = 1:length(scores))


# M3C
out=M3C(mydata=t(simul$data), clusteralg = "hc", maxK=10, pItem = 0.5, repsreal = 100, seed=1)
out$scores
Heatmap(out$realdataresults[[5]]$consensus_matrix)
ClusteringPerformance(theta=out$realdataresults[[5]]$assignments, theta_star=simul$theta)





