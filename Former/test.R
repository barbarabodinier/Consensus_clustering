rm(list=ls())
setwd("~/Dropbox/Consensus_clustering")

library(corpcor)
library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(FactoMineR)
library(diceR)
library(ConsensusClusterPlus)
library(M3C)

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for(name in r) eval(parse(text=paste0(name, '<-sharp:::', name)))

for (filename in list.files("Scripts/Functions/")){
  source(paste0("Scripts/Functions/", filename))
}

# https://www.stat.cmu.edu/~arinaldo/Teaching/36710/F18/Scribed_Lectures/Nov12.pdf
set.seed(3)
n=rep(30,5)
simul=SimulateClustering(n=n,
                         pk=10, nu_xc=1,
                         ev_xc=0.5)
plot(simul)
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
stab=Clustering(simul$data, implementation = HierarchicalClustering, nc = 1:10, K=100)

nc=2
M=ConsensusMatrix(stab, argmax_id = nc)
K=100
q=sum(M)/2
N=nrow(M)*(nrow(M)-1)/2
Heatmap(M)
plot(sh_clust)
L=2

Score(M=M, L=6)

s=matrix(NA,10,10)
for (nc in 2:10){
  M=ConsensusMatrix(stab, argmax_id = nc)
  for (L in 2:10){
    print(paste(nc," - ",L))
    s[nc, L]=Score(M=M, L=L)
  }
}
s[is.infinite(s)]=NA
Heatmap(s)
plot(apply(s,1,max,na.rm=TRUE))

tmp=apply(s,1,max,na.rm=TRUE)
tmp[is.infinite(tmp)]=NA
nc_hat=which.max(tmp)
M=ConsensusMatrix(stab, argmax_id = nc_hat)
sh_clust=hclust(d=as.dist(1-M), method = "complete")
L_hat=which.max(s[nc_hat,])
membership=cutree(sh_clust, k=L_hat)
table(membership)
ClusteringPerformance(theta=membership, theta_star = simul)

M=ConsensusMatrix(stab, argmax_id = 4)
sh_clust=hclust(d=as.dist(1-M), method = "complete")
membership=cutree(sh_clust, k=4)
ClusteringPerformance(theta=membership, theta_star = simul)

M=ConsensusMatrix(stab, argmax_id = 5)
sh_clust=hclust(d=as.dist(1-M), method = "complete")
membership=cutree(sh_clust, k=5)
ClusteringPerformance(theta=membership, theta_star = simul)

M=ConsensusMatrix(stab, argmax_id = 4)
sh_clust=hclust(d=as.dist(1-M), method = "complete")
membership=cutree(sh_clust, k=4)
ClusteringPerformance(theta=membership, theta_star = simul)


Score=function(M, L){
  K=100
  # q=sum(M)/2
  # N=nrow(M)*(nrow(M)-1)/2
  # print(q/N)
  sh_clust=hclust(d=as.dist(1-M), method = "complete")
  sm=cutree(sh_clust, k=L)
  H=matrix(0,L,L)
  for (l in 1:L){
    ids1=which(sm==l)
    for (k in 1:L){
      # print(paste(k, "-", l))
      if (k==l){
        # Case k==l
        # H[l,l]=sum(M[ids1,ids1]/2*K)
        tmpM=M[ids1,ids1]
        H[l,l]=mean(tmpM[upper.tri(tmpM)])
      } else {
        # Case k!=l
        ids2=which(sm==k)
        # H[l,k]=H[k,l]=sum(M[ids1,ids2]*K)
        tmpM=M[ids1,ids2]
        H[k,l]=H[l,k]=mean(tmpM)
      }
    }
  }
  
  # q=sum(H[upper.tri(H,diag = TRUE)])
  # N=sum(upper.tri(H,diag = TRUE))
  # print(q)
  # print(N)
  
  q=sum(M[upper.tri(M)])
  N=sum(upper.tri(M))
  
  loglik=0
  for (l in 1:L){
    ids1=which(sm==l)
    for (k in 1:L){
      print(paste(k, "-", l))
      if (k==l){
        # Probability of observing something at least as extreme under assumption of equal probabilities (>=)
        # thr=floor(H[l,k]*K)
        # logp <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
        
        # thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2*K)
        # logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2*K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
        
        thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2)
        logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2, prob = 1/L, lower.tail = FALSE, log.p = TRUE)
      } else {
        ids2=which(sm==k)
        # Probability of observing something at least as extreme under assumption of equal probabilities (<=)
        # thr=ceiling(H[l,k]*K)
        # logp <- stats::pbinom(thr, size = K, prob = q / N, log.p = TRUE)
        
        # thr=ceiling(H[l,k]*length(ids1)*length(ids2)*K)
        # logp <- stats::pbinom(thr, size = length(ids1)*length(ids2)*K, prob = q / N, log.p = TRUE)
        thr=ceiling(H[l,k]*length(ids1)*length(ids2))
        logp <- stats::pbinom(thr, size = length(ids1)*length(ids2), prob = 1/L, log.p = TRUE)
      }
      print(logp)
      
      # Adding to the log-likelihood
      loglik=loglik+logp
    }
  }
  Sc=-loglik
  return(Sc)
}

# Score=function(M, L){
#   K=100
#   # q=sum(M)/2
#   # N=nrow(M)*(nrow(M)-1)/2
#   # print(q/N)
#   sh_clust=hclust(d=as.dist(1-M), method = "complete")
#   sm=cutree(sh_clust, k=L)
#   H=matrix(0,L,L)
#   for (l in 1:L){
#     ids1=which(sm==l)
#     for (k in 1:L){
#       # print(paste(k, "-", l))
#       if (k==l){
#         # Case k==l
#         # H[l,l]=sum(M[ids1,ids1]/2*K)
#         tmpM=M[ids1,ids1]
#         H[l,l]=mean(tmpM[upper.tri(tmpM)])
#       } else {
#         # Case k!=l
#         ids2=which(sm==k)
#         # H[l,k]=H[k,l]=sum(M[ids1,ids2]*K)
#         tmpM=M[ids1,ids2]
#         H[k,l]=H[l,k]=mean(tmpM)
#       }
#     }
#   }
# 
#   # q=sum(H[upper.tri(H,diag = TRUE)])
#   # N=sum(upper.tri(H,diag = TRUE))
#   # print(q)
#   # print(N)
# 
#   q=sum(M[upper.tri(M)])
#   N=sum(upper.tri(M))
# 
#   loglik=0
#   for (l in 1:L){
#     ids1=which(sm==l)
#     for (k in 1:L){
#       print(paste(k, "-", l))
#       if (k==l){
#         # Probability of observing something at least as extreme under assumption of equal probabilities (>=)
#         # thr=floor(H[l,k]*K)
#         # logp <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
# 
#         # thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2*K)
#         # logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2*K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
# 
#         thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2)
#         logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2, prob = q / N, lower.tail = FALSE, log.p = TRUE)
#       } else {
#         ids2=which(sm==k)
#         # Probability of observing something at least as extreme under assumption of equal probabilities (<=)
#         # thr=ceiling(H[l,k]*K)
#         # logp <- stats::pbinom(thr, size = K, prob = q / N, log.p = TRUE)
# 
#         # thr=ceiling(H[l,k]*length(ids1)*length(ids2)*K)
#         # logp <- stats::pbinom(thr, size = length(ids1)*length(ids2)*K, prob = q / N, log.p = TRUE)
#         thr=ceiling(H[l,k]*length(ids1)*length(ids2))
#         logp <- stats::pbinom(thr, size = length(ids1)*length(ids2), prob = q / N, log.p = TRUE)
#       }
#       print(logp)
# 
#       # Adding to the log-likelihood
#       loglik=loglik+logp
#     }
#   }
#   Sc=-loglik
#   return(Sc)
# }


# Score=function(M, L){
#   K=100
#   # q=sum(M)/2
#   # N=nrow(M)*(nrow(M)-1)/2
#   # print(q/N)
#   sh_clust=hclust(d=as.dist(1-M), method = "complete")
#   sm=cutree(sh_clust, k=L)
#   H=matrix(0,L,L)
#   # for (l in 1:L){
#   ids1=which(sm==l)
#   for (k in 1:L){
#     # print(paste(k, "-", l))
#     if (k==l){
#       # Case k==l
#       # H[l,l]=sum(M[ids1,ids1]/2*K)
#       tmpM=M[ids1,ids1]
#       H[l,l]=mean(tmpM[upper.tri(tmpM)])
#     } else {
#       # Case k!=l
#       ids2=which(sm==k)
#       # H[l,k]=H[k,l]=sum(M[ids1,ids2]*K)
#       tmpM=M[ids1,ids2]
#       H[k,l]=H[l,k]=mean(tmpM)
#     }
#   }
#   # }
#   
#   # q=sum(H[upper.tri(H,diag = TRUE)])
#   # N=sum(upper.tri(H,diag = TRUE))
#   # print(q)
#   # print(N)
#   
#   q=sum(M[upper.tri(M)])
#   N=sum(upper.tri(M))
#   
#   loglik=0
#   for (l in 1:L){
#     ids1=which(sm==l)
#     for (k in 1:L){
#       print(paste(k, "-", l))
#       if (k==l){
#         # Probability of observing something at least as extreme under assumption of equal probabilities (>=)
#         # thr=floor(H[l,k]*K)
#         # logp <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
#         
#         # thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2*K)
#         # logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2*K, prob = q / N, lower.tail = FALSE, log.p = TRUE)
#         
#         # thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2)
#         # logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2, prob = q / N, lower.tail = FALSE, log.p = TRUE)
#         
#         thr=floor(H[l,k]*length(ids1)*(length(ids1)-1)/2)
#         logp <- stats::pbinom(thr - 1, size = length(ids1)*(length(ids1)-1)/2, prob = q / N, lower.tail = FALSE, log.p = TRUE)
#       } else {
#         ids2=which(sm==k)
#         # Probability of observing something at least as extreme under assumption of equal probabilities (<=)
#         # thr=ceiling(H[l,k]*K)
#         # logp <- stats::pbinom(thr, size = K, prob = q / N, log.p = TRUE)
#         
#         # thr=ceiling(H[l,k]*length(ids1)*length(ids2)*K)
#         # logp <- stats::pbinom(thr, size = length(ids1)*length(ids2)*K, prob = q / N, log.p = TRUE)
#         thr=ceiling(H[l,k]*length(ids1)*length(ids2))
#         logp <- stats::pbinom(thr, size = length(ids1)*length(ids2), prob = q / N, log.p = TRUE)
#       }
#       print(logp)
#       
#       # Adding to the log-likelihood
#       loglik=loglik+logp
#     }
#   }
#   Sc=-loglik
#   return(Sc)
# }



q=NULL
for (nc in 2:10){
  M=ConsensusMatrix(stab, argmax_id = nc)
  q=c(q, sum(M))
}

# # Probability of observing a selection count below thr_down under the null (uniform selection)
# p_0 <- stats::pbinom(thr - 1, size = K, prob = q / N, log.p = TRUE) # proportion < pi
# 
# # Probability of observing a selection count above thr_up under the null
# p_1 <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE) # proportion >= pi



mygraph=Graph(Adjacency(stab), node_colour=simul$theta)
par(mar=rep(0,4))
plot(mygraph, layout=layout_with_fr(mygraph))
Heatmap(Adjacency(stab))


A=ConsensusMatrix(stab, argmax_id=5)
Heatmap(A)
D=diag(rowSums(A)/2)
L=D-A
eig=eigen(L)
eig$vectors
plot(eig$values)

# M=ConsensusMatrix(stab)
# 
# K=stab$params$K
# H=K*M
# 
# sh_clust=hclust(as.dist(1-M))
# plot(sh_clust)
# pi=0.9
# 
# M=ConsensusMatrix(stab, argmax_id = 3)
# 
# pi_list = seq(0.51, 0.99, by = 0.01)
# scores=matrix(NA, nrow=dim(stab$coprop)[3], ncol=length(pi_list))
# for (k in 1:dim(stab$coprop)[3]){
#   M=ConsensusMatrix(stab, argmax_id = k)
#   scores[k,]=ClusteringStability(M, pi_list = pi_list)
# }
# which(scores==max(scores, na.rm=TRUE), arr.ind = TRUE)
# Heatmap(t(scores))

pi_list=seq(0.51, 0.99, by = 0.01)
metrics=StabilityMetrics(selprop=stab$coprop, pi_list=pi_list, clustering=TRUE)
mat=metrics$S_2d
rownames(mat)=stab$nc
colnames(mat)=pi_list
Heatmap(t(mat))
argmax_id=which(metrics$S_2d==max(metrics$S_2d, na.rm=TRUE), arr.ind = TRUE)
Heatmap(ConsensusMatrix(stability = stab, argmax_id = argmax_id))
Heatmap(ConsensusMatrix(stability = stab, argmax_id = 5))

# Heatmap(ConsensusMatrix(stability = stab, argmax_id = 7))
sh_clust=hclust(as.dist(1-ConsensusMatrix(stability = stab, argmax_id = argmax_id)), method="complete")
plot(sh_clust)
A=ifelse(ConsensusMatrix(stability = stab, argmax_id = argmax_id)>=pi_list[argmax_id[2]], yes=1, no=0)
Heatmap(A)
Z=cutree(sh_clust, h=1-pi_list[argmax_id[2]])
sh_clust$height
Zbis=cutree(sh_clust, k=argmax_id[1])

sh_clust=hclust(as.dist(1-ConsensusMatrix(stability = stab, argmax_id = 5)), method="complete")
Zthree=cutree(sh_clust, k=5)


ClusteringPerformance(theta=Z, theta_star=simul)
ClusteringPerformance(theta=Zbis, theta_star=simul)
ClusteringPerformance(theta=Zthree, theta_star=simul)


scores=NULL
nc_list=2:10
for (nc in 2:10){
  scores=c(scores, ClusteringStability(coprop=ConsensusMatrix(stability = stab, argmax_id = nc), nc=nc))
}
names(scores)=nc_list
plot(as.numeric(names(scores)), scores)

argmax_id=as.numeric(names(scores)[which.max(scores)])
sh_clust=hclust(as.dist(1-ConsensusMatrix(stability = stab, argmax_id = argmax_id)), method="complete")
Z=cutree(sh_clust, k=argmax_id)
Heatmap(ConsensusMatrix(stability = stab, argmax_id = argmax_id))
Heatmap(ConsensusMatrix(stability = stab, argmax_id = 5))

Heatmap(t(stab$Sc_2d))
ArgmaxId(stab, clustering=TRUE)

ClusteringPerformance(theta=Z, theta_star=simul)


# M3C
out=M3C(mydata=t(simul$data), clusteralg = "hc", maxK=10, pItem = 0.5, repsreal = 100, seed=1)
out$scores
Heatmap(out$realdataresults[[5]]$consensus_matrix)
ClusteringPerformance(theta=out$realdataresults[[5]]$assignments, theta_star=simul$theta)




