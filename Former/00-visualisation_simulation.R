rm(list=ls())
setwd("~/Dropbox/Consensus_clustering")

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(FactoMineR)
library(corpcor)

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
n=c(20,50,30)
simul=SimulateClustering(n=n, pk=10, 
                         nu_xc = 1, ev_xc=0.3)

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

# Simulation of data with clusters
dev.off()
set.seed(1)
n=c(20,50,30)
simul=SimulateClustering(n=n, pk=10, 
                         theta_xc = c(1,1, rep(0,8)),
                         ev_xc=0.8)
plot(simul$data[,1], simul$data[,2], col=simul$theta)
plot(simul$data[,1], simul$data[,3], col=simul$theta)
plot(simul$data[,3], simul$data[,4], col=simul$theta)
par(mar=c(5,5,1,5))
Heatmap(cor(simul$data), col = c("blue","white","red"))


# Simulation of data with clusters
dev.off()
set.seed(1)
n=c(20,50,30)
pk=10
A=matrix(0,ncol=sum(pk), nrow=sum(pk))
A[1,2]=1
A=A+t(A)
simul=SimulateClustering(n=n, pk=pk, 
                         theta_xc = c(1,1, rep(0,8)),
                         adjacency = A,
                         ev_xc=0.8, output_matrices = TRUE)
par(mar=c(5,5,1,5))
Heatmap(cor(simul$data), col = c("blue","white","red"))
plot(simul$data[,1], simul$data[,2], col=simul$theta)
plot(simul$data[,1], simul$data[,3], col=simul$theta)
plot(simul$data[,3], simul$data[,4], col=simul$theta)

simul$phi
cor2pcor(cov2cor(simul$sigma))

apply(simul$data, 2, mean)
cor(simul$data)
cor2pcor(cor(simul$data))
simul$data[,1]=simul$data[,1]+5
simul$data[,2]=simul$data[,2]+10
simul$data[,3]=simul$data[,3]+2
apply(simul$data, 2, mean)
cor2pcor(cor(simul$data))
cor(simul$data)




