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

# Simulation of data with clusters
set.seed(1)
n=c(20,50,30)
simul=SimulateClustering(n=n,
                         pk=10, nu_xc=1,
                         ev_xc=0.5)
x=simul$data

# Consensus clustering
stab=Clustering(xdata=x, implementation=HierarchicalClustering)

# Initialising figure
{pdf("Figures/Illustration.pdf", width=23, height=20)
  layout(mat=matrix(c(1,5,9,13,
                      2,6,10,14,
                      3,7,11,15,
                      4,8,12,16), ncol=4, byrow = TRUE))
  
  # Distance metric on different subsamples
  par(mar=c(5,5,1,7))
  for (k in 1:3){
    set.seed(k)
    s=Resample(data=x)
    xsub=x[s,]
    mydist=as.matrix(dist(xsub, upper = TRUE))
    mydist_full=matrix(NA, ncol=nrow(x), nrow=nrow(x))
    rownames(mydist_full)=colnames(mydist_full)=rownames(x)
    for (i in 1:nrow(mydist_full)){
      row_id=rownames(mydist_full)[i]
      if (row_id %in% rownames(mydist)){
        for (j in 1:ncol(mydist_full)){
          col_id=colnames(mydist_full)[j]
          if (col_id %in% colnames(mydist)){
            mydist_full[row_id, col_id]=mydist[row_id, col_id]
          }
        } 
      }
    }
    Heatmap(mydist_full, legend_range = c(0,9),
            legend=ifelse(k==3, yes=TRUE, no=FALSE))
    
    for (l in c(2,3,4)){
      myhclust=hclust(as.dist(mydist), method="complete")
      members=CoMembership(cutree(myhclust, k=l))
      members_full=matrix(NA, ncol=nrow(x), nrow=nrow(x))
      rownames(members_full)=colnames(members_full)=rownames(x)
      for (i in 1:nrow(members_full)){
        row_id=rownames(members_full)[i]
        if (row_id %in% rownames(members)){
          for (j in 1:ncol(members_full)){
            col_id=colnames(members_full)[j]
            if (col_id %in% colnames(members)){
              members_full[row_id, col_id]=members[row_id, col_id]
            }
          } 
        }
      }
      Heatmap(members_full, col=c("royalblue", "darkred"), legend = FALSE)
    }
  }
  
  # Consensus matrices
  plot.new()
  for (l in c(2,3,4)){
    Heatmap(ConsensusMatrix(stab, argmax_id=l))
  }
  dev.off()}


# Calibration plot
max_N=10
{pdf("Figures/Illustration_calibration.pdf", width=15, height=5)
par(mar=c(1,5,5,1))
plot(stab$nc[1:max_N], stab$Sc[1:max_N], pch=19, col="navy",
     panel.first=c(abline(h=stab$Sc, lty=3, col="grey"),
                   abline(v=stab$nc, lty=3, col="grey")),
     xlab="", ylab="Consensus Score", xaxt="n",
     cex.lab=1.5, las=3)
axis(side = 3, at = stab$nc[1:max_N], las=2)
mtext(text="Number of clusters", side = 3, line=3, cex=1.5)
dev.off()}

par(mar=c(5,5,1,10))
CalibrationPlot(stab, clustering = TRUE)
Argmax(stab, clustering=TRUE)

# Consensus matrix
par(mar=c(5,5,1,5))
mat=ConsensusMatrix(stab)
Heatmap(mat)

# Stable co-membership
adjacency=Adjacency(stab)
Heatmap(adjacency, colours = c("ivory", "darkred"), legend = FALSE)

# Network representation
set.seed(1)
mycolours=distinctColorPalette(k=length(n))
mygraph=Graph(adjacency = adjacency,
              node_colour = lighten(mycolours, amount=0.5)[simul$theta],
              satellites = TRUE)
V(mygraph)$size=3

# Stable clusters
myclusters=Clusters(adjacency=adjacency, implementation = igraph::cluster_walktrap)

set.seed(1)
mycolours_estimated=distinctColorPalette(k=length(unique(myclusters)))
par(mar=rep(0,4))
V(mygraph)$frame.color=mycolours_estimated[myclusters]
# V(mygraph)$frame.width=10
# plot(mygraph,
#      layout=layout_with_kk(mygraph, weights=ifelse(crossing(mycommunities, mygraph), 10, 1)))
# plot(mygraph,
#      layout=layout_with_fr(mygraph, weights=ifelse(crossing(mycommunities, mygraph), 1, 10)))


# Network representation
par(mar=rep(0,4))
mycommunities=cluster_louvain(mygraph)
mycommunities=edge.betweenness.community(mygraph)
mycommunities=cluster_leading_eigen(mygraph)
mycommunities=cluster_walktrap(mygraph)
plot(mycommunities, mygraph, 
     col=V(mygraph)$color,
     mark.border = "red", 
     mark.col=NA, 
     layout=layout_with_kk(mygraph, weights=ifelse(crossing(mycommunities, mygraph), 1, 0.5)))
plot(mycommunities, mygraph, 
     col=V(mygraph)$color,
     mark.border = "red", 
     mark.col=NA, 
     layout=layout_with_fr(mygraph, weights=ifelse(crossing(mycommunities, mygraph), 1, 100)))

# Performance as a function of stability score
implementations=c(igraph::components, igraph::cluster_louvain, igraph::cluster_walktrap)
metric="adj_rand"
par(mar=c(5,5,1,1), mfrow=c(1,length(implementations)), xpd=FALSE)
for (implementation in implementations){
  print(implementation)
  rand_index=matrix(NA, nrow=length(stab$nc), ncol=length(stab$params$pi_list))
  for (i in 1:length(stab$nc)){
    for (j in 1:length(stab$params$pi_list)){
      myclusters=Clusters(adjacency = stab, argmax_id = matrix(c(i,j), nrow=1), 
                          implementation=implementation)
      rand_index[i,j]=ClusteringPerformance(theta=myclusters, theta_star=simul)[,metric]
    }
  }
  
  plot(as.vector(stab$Sc_2d), 
       as.vector(rand_index),
       cex.lab=1.5, ylim=c(min(as.vector(rand_index)), 1), 
       pch=17, bty="n", las=1, col="slategrey",
       panel.first=c(abline(v=axisTicks(usr=c(0,max(as.vector(stab$Sc_2d), na.rm=TRUE)), log = FALSE), lty=3, col="grey"),
                     abline(h=axisTicks(usr=c(range(as.vector(rand_index), na.rm=TRUE)), log = FALSE), lty=3, col="grey"),
                     abline(v=max(as.vector(stab$Sc_2d), na.rm=TRUE), 
                            lty=3, col="darkred"),
                     abline(h=as.vector(rand_index)[which.max(as.vector(stab$Sc_2d))], 
                            lty=3, col="darkred")),
       xlab="Stability score", ylab="Adjusted Rand index")
  points(max(as.vector(stab$Sc_2d), na.rm=TRUE), 
         as.vector(rand_index)[which.max(as.vector(stab$Sc_2d))],
         pch=17, col="darkred")
}

# Comparison with hierarchical clustering
myhclust=hclust(d=dist(x), method = 'complete')
par(mar=c(5,5,2,1))
plot(myhclust, hang=-1, xlab="", main="", sub = "", 
     cex.lab=1.5, yaxt="n", labels=FALSE)
axis(side=2, at=seq(0,10,by=2), las=1)
for (i in 1:nrow(x)){
  axis(side=1, at=i, labels=myhclust$labels[myhclust$order[i]], 
       line=-1, las=2, 
       col.axis=darken(mycolours[simul$theta[myhclust$order[i]]], amount=0.5),
       col=mycolours[simul$theta[myhclust$order[i]]])
}

# Clustering performance
myclusters=cutree(myhclust, k=3)
ClusteringPerformance(theta=myclusters, theta_star=simul)$adj_rand

ClusteringPerformance(theta=stab, theta_star=simul, implementation=igraph::cluster_walktrap)$adj_rand


