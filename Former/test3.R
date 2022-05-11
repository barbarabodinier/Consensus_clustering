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
set.seed(2)
n=rep(20,5)
simul=SimulateClustering(n=n, pk=20, 
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

stab=Clustering(xdata=simul$data, implementation=HierarchicalClustering, tau=0.5)

scores=rep(NA, 10)
for (nc in 2:dim(stab$coprop)[3]){
  print(nc)
  scores[nc]=Score(ConsensusMatrix(stab, argmax_id = nc), nc=nc, pi=0.9)
  cat("\n")
}
plot(scores)

Heatmap(ConsensusMatrix(stab, argmax_id = which.max(scores)))
Heatmap(ConsensusMatrix(stab, argmax_id = 5))

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  loglik=0
  for (cluster_id in 1:nc){
    ids=which(theta==cluster_id)
    if (length(ids)>1){
      tmpM=M[ids,ids]
      tmpM=tmpM[upper.tri(tmpM)]
      # N=length(tmpM)
      N=1
      # pi=min(tmpM)
      # print(pi)
      log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=2)
      loglik=loglik+1/N*(sum(tmpM<pi)*log_proba$p_0+sum(tmpM>=pi)*log_proba$p_1)
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  loglik=0
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if (theta[i]==theta[j]){
        log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=2)
        loglik=loglik+ifelse(M[i,j]>=pi, 
                             yes=log_proba$p_1, 
                             no=log_proba$p_0)
      } else {
        log_proba=BinomialProbabilities(q=nc-1, N=nc, pi=pi, K=100, n_cat=2)
        loglik=loglik+ifelse(M[i,j]>=pi, 
                             yes=log_proba$p_1, 
                             no=log_proba$p_0)
      }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  loglik=0
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if (theta[i]==theta[j]){
        log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        if (M[i,j]>=pi){
          loglik=loglik+log_proba$p_3
        }
        if (M[i,j]<=(1-pi)){
          loglik=loglik+log_proba$p_1
        }
        if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
          loglik=loglik+log_proba$p_2
        }
      } else {
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        # if (M[i,j]>=pi){
        #   loglik=loglik+log_proba$p_3
        # }
        # if (M[i,j]<=(1-pi)){
        #   loglik=loglik+log_proba$p_1
        # }
        # if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
        #   loglik=loglik+log_proba$p_2
        # }
      }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  loglik=0
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if (theta[i]==theta[j]){
        # print(table(theta)[as.character(theta[i])]/length(theta))
        log_proba=BinomialProbabilities(q=(table(theta)[as.character(theta[i])]/length(theta))^2, N=nc, pi=pi, K=100, n_cat=3)
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        if (M[i,j]>=pi){
          loglik=loglik+log_proba$p_3
        }
        if (M[i,j]<=(1-pi)){
          loglik=loglik+log_proba$p_1
        }
        if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
          loglik=loglik+log_proba$p_2
        }
      } else {
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        # if (M[i,j]>=pi){
        #   loglik=loglik+log_proba$p_3
        # }
        # if (M[i,j]<=(1-pi)){
        #   loglik=loglik+log_proba$p_1
        # }
        # if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
        #   loglik=loglik+log_proba$p_2
        # }
      }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  print(table(theta))
  print(sum((table(theta)/length(theta))^2),)
  loglik=0
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if (theta[i]==theta[j]){
        # print(table(theta)[as.character(theta[i])]/length(theta))
        log_proba=BinomialProbabilities(q=sum((table(theta)/length(theta))^2), N=nc, pi=pi, K=100, n_cat=3)
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        if (M[i,j]>=pi){
          loglik=loglik+log_proba$p_3
        }
        if (M[i,j]<=(1-pi)){
          loglik=loglik+log_proba$p_1
        }
        if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
          loglik=loglik+log_proba$p_2
        }
      } else {
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        # if (M[i,j]>=pi){
        #   loglik=loglik+log_proba$p_3
        # }
        # if (M[i,j]<=(1-pi)){
        #   loglik=loglik+log_proba$p_1
        # }
        # if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
        #   loglik=loglik+log_proba$p_2
        # }
      }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  print(table(theta))
  print(sum((table(theta)/length(theta))^2),)
  loglik=0
  for (i in 1:(nrow(M)-1)){
    for (j in (i+1):ncol(M)){
      if (theta[i]==theta[j]){
        # print(table(theta)[as.character(theta[i])]/length(theta))
        log_proba=BinomialProbabilities(q=sum((table(theta)/length(theta))^2), N=nc, pi=pi, K=100, n_cat=3)
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        if (M[i,j]>=pi){
          loglik=loglik+log_proba$p_3
        }
        if (M[i,j]<=(1-pi)){
          loglik=loglik+log_proba$p_1
        }
        if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
          loglik=loglik+log_proba$p_2
        }
      } else {
        # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
        # if (M[i,j]>=pi){
        #   loglik=loglik+log_proba$p_3
        # }
        # if (M[i,j]<=(1-pi)){
        #   loglik=loglik+log_proba$p_1
        # }
        # if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
        #   loglik=loglik+log_proba$p_2
        # }
      }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.99){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  print(table(theta))
  print(sum((table(theta)/length(theta))^2))
  loglik=0
  for (i in 1:(nrow(M)-1)){
    for (j in (i+1):ncol(M)){
      # if (theta[i]==theta[j]){
      # print(table(theta)[as.character(theta[i])]/length(theta))
      log_proba=BinomialProbabilities(q=sum((table(theta)/length(theta))^2), N=1, pi=pi, K=100, n_cat=3)
      # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
      if (M[i,j]>=pi){
        loglik=loglik+log_proba$p_3
      }
      if (M[i,j]<=(1-pi)){
        loglik=loglik+log_proba$p_1
      }
      if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
        loglik=loglik+log_proba$p_2
      }
      # } else {
      # log_proba=BinomialProbabilities(q=1, N=nc, pi=pi, K=100, n_cat=3)
      # if (M[i,j]>=pi){
      #   loglik=loglik+log_proba$p_3
      # }
      # if (M[i,j]<=(1-pi)){
      #   loglik=loglik+log_proba$p_1
      # }
      # if ((M[i,j]>(1-pi))&((M[i,j]<pi))){
      #   loglik=loglik+log_proba$p_2
      # }
      # }
    }
  }
  return(-loglik)
}

Score=function(M, nc, pi=0.9){
  sh_clust=hclust(as.dist(1-M), method = "complete")
  theta=cutree(sh_clust, k=nc)
  print(table(theta))
  print(sum((table(theta)/length(theta))^2))
  
  n_1=sum(M[upper.tri(M)]<=(1-pi))
  n_2=sum((M[upper.tri(M)]>(1-pi))&(M[upper.tri(M)]<pi))
  n_3=sum(M[upper.tri(M)]>=pi)
  N=nrow(M)*(nrow(M)-1)/2
  if ((n_1+n_2+n_3-N)!=0){
    stop("Inconsistent.")
  }
  print(n_1)
  print(n_2)
  print(n_3)
  
  log_proba=BinomialProbabilities(q=sum((table(theta)/length(theta))^2), N=1, pi=pi, K=100, n_cat=3)
  
  print(log_proba$p_1)
  print(log_proba$p_2)
  print(log_proba$p_3)
  
  loglik=n_1*log_proba$p_1+n_2*log_proba$p_2+n_3*log_proba$p_3

  return(-loglik)
}



