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
                         pk=100, 
                         nu_xc=0.2,
                         ev_xc=0.8)
x=simul$data

# Heatmap of pairwise distances
par(mfrow=c(1,2), mar=c(5,5,1,5))
thr=max(as.matrix(dist(simul$data)))
Heatmap(as.matrix(dist(simul$data)), 
        legend_range = c(0,thr))
Heatmap(as.matrix(dist(simul$data[,which(simul$theta_xc==1)])),
        legend_range = c(0,thr))

# Consensus clustering
stab=Clustering(xdata=simul$data, implementation = HierarchicalClustering)

# mycolours=colorRampPalette(c("gold","darkred"))(10)
# plot(stab$coprobnull[1,], ylim=c(0,1))
# for (i in 1:10){
#   lines(stab$coprobnull[i,], col=mycolours[i])
# }
# 
# id=5
# coprop=stab$coprop[,,id]
# coprobnull=stab$coprobnull[id,]
# 
# 
# dpoisbinom(x=1:100, pp = coprobnull, log_d = FALSE)
# dpbinom(x=1:100, probs = coprobnull, log = FALSE)
# 
# ConsensusScore(coprop=coprop, nc=5, coprobnull = coprobnull)
# 
# K=length(coprobnull)
# loglik=0
# tmp=0
# for (i in 1:(nrow(coprop)-1)){
#   for (j in (i+1):nrow(coprop)){
#     tmp=tmp+1
#     # loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
#     loglik=loglik+dpoisbinom(x=round(coprop[i,j]*K), pp = coprobnull, log_d = TRUE)
#   }
# } # XXX maybe a problem with different cluster sizes (ratio above 1)







plot(stab$Sc)
Heatmap(ConsensusMatrix(stab))
which.max(stab$Sc)
ClusteringPerformance(theta=stab, theta_star=simul)

# Sparse consensus clustering
stab=Clustering(xdata=simul$data, 
                implementation = SparseHierarchicalClustering, 
                Lambda=LambdaSequence(lmax=10, lmin=1+1e-6, cardinal = 20))
Heatmap(ConsensusMatrix(stab, argmax_id = 5))




coprop=ConsensusMatrix(stab, argmax_id=3)
Heatmap(coprop)
coprobnull=stab$coprobnull[3,]
plot(density(coprobnull))
K=length(coprobnull)

plot(0:K, PoissonBinomial::dpbinom(x=0:K, probs = coprobnull, log=TRUE), type="b")
points(0:K, dbinom(x=0:K, size=K, prob=mean(coprobnull), log=TRUE), col="red", type="b")





# Calculating log-likelihood for observed consensus matrix (under expected p under the null)
loglik=0
tmp=0
for (i in 1:(nrow(coprop)-1)){
  for (j in (i+1):nrow(coprop)){
    tmp=tmp+1
    # loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
    loglik=loglik+PoissonBinomial::dpbinom(x=round(coprop[i,j]*K), probs = coprobnull, log = TRUE)
  }
}
print(loglik)

coprobnull=sort(coprobnull)
n_exp=rep(0, K+1)
names(n_exp)=as.character(K:0)
N=sum(upper.tri(coprop))

n_exp[1]=coprobnull[1]*N
for (i in 1:(K-1)){
  n_exp[i+1]=(coprobnull[i+1]-coprobnull[i])*N
}
n_exp[K+1]=N-sum(n_exp)

best_score=0
for (i in 1:length(n_exp)){
  best_score=best_score+n_exp[i]*PoissonBinomial::dpbinom(x=as.numeric(names(n_exp[i])), probs = coprobnull, log = TRUE)
}



# Calculating log-likelihood for most stable (binary) consensus matrix (under the expected p under the null)
N=length(theta)*(length(theta)-1)/2
p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
n_1=p_unif*N
n_3=(1-p_unif)*N
print(n_1)
print(n_3)
p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
best_score=(n_1*p_1+n_3*p_3)
print(best_score)






stabtest=Clustering(xdata=simul$data, 
                implementation = SparseHierarchicalClustering, 
                Lambda=1e10)
Heatmap(ConsensusMatrix(stabtest, argmax_id = 5))
stab_nonsparse=Clustering(xdata=simul$data, 
                          implementation = HierarchicalClustering)
Heatmap(ConsensusMatrix(stab_nonsparse, argmax_id = 5))

mycolours=colorRampPalette(c("gold","darkred"))(3)
plot(stab$coprobnull[1,], ylim=c(0,1))
for (i in 2:3){
  lines(stab$coprobnull[i,], col=mycolours[i])
}

stab$Sc[2:3]
stab$coprobnull[2:3,]

stab$Q[2]

plot(stab$coprobnull[2,], ylim=c(0,1), type="l", col="darkred")
lines(stab$coprobnull[3,], col="blue")
lines(stab$coprobnull[4,], col="forestgreen")
lines(stab$coprobnull[5,], col="red")
lines(stab$coprobnull[6,], col="navy")
lines(stab$coprobnull[7,], col="orange")

plot(density(stab$coprobnull[2,], from=0, to=1), col="darkred", ylim=c(0,25))
lines(density(stab$coprobnull[3,], from=0, to=1), col="blue")
lines(density(stab$coprobnull[4,], from=0, to=1), col="forestgreen")
lines(density(stab$coprobnull[5,], from=0, to=1), col="red")
lines(density(stab$coprobnull[6,], from=0, to=1), col="navy")
lines(density(stab$coprobnull[7,], from=0, to=1), col="orange")

ConsensusScore(coprop=ConsensusMatrix(stab, argmax_id=2), nc = 2, coprobnull = stab$coprobnull[2,])


plot(stab$Lambda, stab$Q)
# plot(stab$Sc)
# mat=matrix(NA, nrow=length(unique(stab$nc)), ncol=length(unique(stab$Lambda)))
# rownames(mat)=unique(stab$nc)
# colnames(mat)=unique(stab$Lambda)
# mat[cbind(as.character(stab$nc), as.character(stab$Lambda))]=stab$Sc
# Heatmap(mat)
# 
# stab$nc[which.max(stab$Sc)]
# stab$Lambda[which.max(stab$Sc)]
# plot(stab$selprop[which.max(stab$Sc),],  
#      pch=19,
#      col=ifelse(simul$theta_xc==1, yes="red", no="grey"))
# Heatmap(ConsensusMatrix(stab))
# Heatmap(ConsensusMatrix(stab, argmax_id=105))
# table(Clusters(stab))
# ClusteringPerformance(theta=stab, theta_star=simul)

# Two-step approach: 1-selection, 2-clustering
par(mfrow=c(1,1))
plot(stab$S, ylim=c(0,1))
lines(stab$Sc)
text(which(!duplicated(stab$Q)), unique(stab$S), labels = unique(stab$Q), col="red")

par(mfrow=c(3,2))
plot(stab$Q, stab$S)
plot(stab$Q, apply(stab$S_2d,1,max))
# plot(stab$Sc, type="l")
ids_best_sel=which(stab$S==max(stab$S, na.rm=TRUE))
# ids_best_sel=which(apply(stab$S_2d,1,max,na.rm=TRUE)==max(apply(stab$S_2d,1,max,na.rm=TRUE)))
Sc_best_sel=stab$Sc[ids_best_sel]
argmax_id=ids_best_sel[which.max(Sc_best_sel)]
stab$nc[argmax_id]
stab$Lambda[argmax_id]

Heatmap(ConsensusMatrix(stab, argmax_id=argmax_id))
plot(stab$selprop[argmax_id,], col=ifelse(simul$theta_xc==1, yes="red", no="grey"))
abline(h=stab$P[argmax_id], col="darkred", lty=2)

ClusteringPerformance(Clusters(stab, argmax_id=argmax_id), theta_star=simul)



par(mfrow=c(1,2))

ConsensusScore(ConsensusMatrix(stab, argmax_id=2), nc=2, K=100)
plot(dbinom(seq(0,100), size=100, prob=0.9866667, log = TRUE), type="h", lwd=3)
points(dbinom(seq(0,100), size=100, prob=0.9390899, log=TRUE), col="red", pch=19)

ConsensusScore(ConsensusMatrix(stab, argmax_id=10), nc=10, K=100)
plot(dbinom(seq(0,100), size=100, prob=0.1746756, log = TRUE), type="h", lwd=3)
points(dbinom(seq(0,100), size=100, prob=0.2014615, log=TRUE), col="red", pch=19)



par(mfrow=c(1,2))
plot(dbinom(seq(0,100), size=100, prob=0.9, log = FALSE), type="h", lwd=3)
points(dbinom(seq(0,100), size=100, prob=0.8, log=FALSE), col="red", pch=19)

plot(dbinom(seq(0,100), size=100, prob=0.9, log = TRUE), type="h", lwd=3)
points(dbinom(seq(0,100), size=100, prob=0.8, log=TRUE), col="red", pch=19)


sum(dbinom(seq(0,100), size=100, prob=0.9, log = FALSE))
sum(dbinom(seq(0,100), size=100, prob=0.8, log = FALSE))

sum(dbinom(seq(0,100), size=100, prob=0.9, log = TRUE))
sum(dbinom(seq(0,100), size=100, prob=0.8, log = TRUE))









coprop=ConsensusMatrix(stab, argmax_id = 2)
nc=2
K=100
linkage="complete"





# Clustering on the consensus matrix
sh_clust=hclust(as.dist(1-coprop), method = linkage)

# Identifying stable clusters
theta=cutree(sh_clust, k=nc)
print(table(theta))

# Probability that items i and j belong to the same cluster
N=length(theta)*(length(theta)-1)/2
p_unif=sum(table(theta)*(table(theta)-1))/(2*N) # P(i) * P(j | i)
print(p_unif)
print(sum(coprop[upper.tri(coprop)])/N)
print(mean(coprop[upper.tri(coprop)]))

mean(coprop[upper.tri(coprop)]*p_unif/(sum(coprop[upper.tri(coprop)])/N))

p_unif



m=round(sum(coprop[upper.tri(coprop)]))
plot(dbinom(seq(K*N*(p_unif-0.01),K*N*(p_unif+0.01)), size = K*N, prob = p_unif), type="h")
plot(dbinom(seq(0,K), size = K, prob = p_unif), type="h")







m=round(sum(coprop[upper.tri(coprop)]))
n=N-m
# p_unif=sum(coprop[upper.tri(coprop)])/N
# print(p_unif)
print(N)
print(m)
print(n)

plot(dbinom(seq(0,100), size = K, prob = p_unif), type="h", lwd=3, ylim=c(0,1))
plot(dbinom(seq(0,100), size = K, prob = sum(coprop[upper.tri(coprop)])/N), type="h", lwd=3, ylim=c(0,1))
points(dhyper(seq(0,100), m = m, n = n, k=K), col="red", pch=19)


# Calculating log-likelihood for observed consensus matrix
loglik=0
tmp=0
for (i in 1:(nrow(coprop)-1)){
  for (j in (i+1):nrow(coprop)){
    tmp=tmp+1
    loglik=loglik+dbinom(round(coprop[i,j]*K), size=K, prob=p_unif, log = TRUE)
  }
} # XXX maybe a problem with different cluster sizes (ratio above 1)
print(loglik)

# Calculating log-likelihood for most stable (binary) consensus matrix
n_1=p_unif*N
n_3=(1-p_unif)*N
print(n_1)
print(n_3)
p_1=dbinom(K, size=K, prob=p_unif, log = TRUE)
p_3=dbinom(0, size=K, prob=p_unif, log = TRUE)
best_score=(n_1*p_1+n_3*p_3)
print(loglik)

# Calculating consensus score as likelihood ratio
score=loglik/best_score












