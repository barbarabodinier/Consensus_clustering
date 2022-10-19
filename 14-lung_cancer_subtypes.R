rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

library(data.table)
library(readxl)
library(dendextend)
library(plotrix)
library(colorspace)
library(annotate)
library(hgu95a.db)
# library(sharp)
devtools::load_all("/Users/barbara/Dropbox/R_packages/Stability/sharp")

# Parameters
tau <- 0.5
noit <- 20
niter <- 10
max_nc <- 20
distance <- "euclidian"

# Loading and preparing the data
mydata <- fread("Data/Lung_cancer/fig1tree.cdt.tsv", data.table = FALSE)
mydata <- mydata[-c(1, 2), ]
rownames(mydata) <- mydata[, 2]
mydata <- mydata[, -c(1:4)]
for (k in 1:ncol(mydata)) {
  mydata[, k] <- as.numeric(mydata[, k])
}
x <- t(mydata)

# Checking median transformation of items
apply(mydata, 2, median)
apply(mydata, 2, FUN = function(z) {
  sum(z^2)
})

# Removing N=2 outliers
x <- x[-c(9, 166), ]

# Defining the classes
covars <- data.frame(
  normal = ifelse(grepl("^NL", rownames(x)), yes = 1, no = 0),
  adenocarcinoma = ifelse(grepl("^AD", rownames(x)), yes = 1, no = 0),
  carcinoid = ifelse(grepl("^COID", rownames(x)), yes = 1, no = 0),
  small_cell = ifelse(grepl("^SMCL", rownames(x)), yes = 1, no = 0),
  squamous = ifelse(grepl("^SQ", rownames(x)), yes = 1, no = 0)
)
subtypes <- sharp:::DummyToCategories(covars)
names(subtypes) <- rownames(x)
print(all(apply(covars, 1, sum) == 1))

# Excluding adenocarcinomas
x <- x[which(subtypes != 2), ]
subtypes <- subtypes[which(subtypes != 2)]

# Loading participant information
samples <- data.frame(read_excel("Data/Lung_cancer/sampledata.xls", sheet = 2))
rownames(samples) <- samples[, 4]
samples$Summary.Stage[which(samples$Summary.Stage == "_")] <- NA
samples$Summary.Stage[which(samples$Summary.Stage == "NA-")] <- NA
samples$Summary.Stage[which(samples$Summary.Stage == "I")] <- NA

# Hierarchical clustering
myhclust <- hclust(d = dist(scale(x), method = distance), method = "complete")
dend <- as.dendrogram(myhclust)
labels_colors(dend) <- subtypes[labels(dend)]

# Consensus hierarchical clustering (unweighted)
stab_unw <- Clustering(
  xdata = x,
  nc = 1:max_nc,
  distance = distance,
  implementation = HierarchicalClustering,
  tau = tau
)
saveRDS(stab_unw, paste0("Results/Application/Consensus_unweighted_hclust_", distance, "_", tau, ".rds"))

# Consensus hierarchical clustering (weighted)
system.time({
  stab_w <- Clustering(
    xdata = x,
    nc = 1:max_nc,
    distance = distance,
    implementation = HierarchicalClustering,
    tau = tau,
    Lambda = LambdaSequence(lmax = 10, lmin = 0.1, cardinal = 10),
    noit = noit,
    niter = niter,
    n_cores = 4
  )
})
table(subtypes, Clusters(stab_w))
saveRDS(stab_w, paste0("Results/Application/Consensus_weighted_hclust_", distance, "_", tau, ".rds"))

# Figure parameters
nc_star <- 4
subtype_colours <- lighten(c("darkolivegreen", NA, "darkred", "navy", "salmon4"), amount = 0.5)
cluster_colours <- lighten(c("tomato", "dodgerblue", "seagreen4", "tan"), amount = 0)

# Dendrogram
{
  pdf("Figures/Dendrogram_hierarchical_lung_subtypes.pdf",
    width = 14, height = 7
  )
  par(mar = c(7, 5, 1, 1))
  dend <- as.dendrogram(myhclust)
  labels_colors(dend) <- darken(subtype_colours[subtypes[labels(dend)]], amount = 0.5)
  labels_cex(dend) <- 0.8
  dend <- color_branches(dend, k = nc_star, col = cluster_colours)
  plot(dend, ylab = "Euclidian distance", cex.lab = 1.5)
  points(1:nrow(x), rep(0, nrow(x)),
    pch = 19, col = subtype_colours[subtypes[labels(dend)]]
  )
  myhclust$height[length(myhclust$height) - nc_star + c(0, 1)]
  abline(
    h = mean(myhclust$height[length(myhclust$height) - nc_star + c(1, 2)]),
    lty = 2, col = "darkred"
  )
  ordered <- cutree(myhclust, k = nc_star)[myhclust$order]
  for (k in 1:nc_star) {
    tmpx <- range(which(ordered == unique(ordered)[k])) + c(-0.25, 0.25)
    axis(
      side = 1, at = tmpx, labels = NA, line = 3,
      col = cluster_colours[k]
    )
    par(xpd = TRUE)
    text(
      x = mean(tmpx), y = par("usr")[3] - 20,
      labels = paste0("Cluster ", k),
      srt = 180, cex = 1.25,
      col = darken(cluster_colours[k], amount = 0.5)
    )
  }
  dev.off()
}

for (type in c("unweighted", "weighted")) {
  if (type == "unweighted") {
    stab <- stab_unw
  } else {
    stab <- stab_w
  }

  # Calibration curves
  {
    pdf(paste0("Figures/Calibration_", type, "_lung_subtypes.pdf"),
      width = 12, height = 12
    )
    par(mar = rep(9, 4))
    CalibrationPlot(stab, xlab = "G", cex.lab = 2, cex.legend = 1.5)
    dev.off()
  }

  # Consensus matrix
  {
    pdf(paste0("Figures/Consensus_", type, "_lung_subtypes.pdf"),
      width = 12, height = 12
    )
    par(mar = rep(9, 4))
    mat <- ConsensusMatrix(stab)
    myhclust_consensus <- hclust(as.dist(1 - mat))
    mat <- mat[myhclust_consensus$order, myhclust_consensus$order]
    Heatmap(mat,
      axes = FALSE,
      legend = ifelse(type == "weighted", yes = TRUE, no = FALSE)
    )
    for (k in 1:ncol(mat)) {
      axis(
        side = 1, at = k - 0.5, las = 2,
        labels = colnames(mat)[k],
        col = subtype_colours[subtypes[colnames(mat)[k]]],
        col.axis = darken(subtype_colours[subtypes[colnames(mat)[k]]], amount = 0.5),
        cex.axis = 0.8
      )
    }
    for (k in 1:ncol(mat)) {
      axis(
        side = 2, at = nrow(mat) - k + 0.5, las = 2,
        labels = colnames(mat)[k],
        col = subtype_colours[subtypes[colnames(mat)[k]]],
        col.axis = darken(subtype_colours[subtypes[colnames(mat)[k]]], amount = 0.5),
        cex.axis = 0.8
      )
    }
    box()
    ordered <- Clusters(stab)[colnames(mat)]
    for (k in 1:nc_star) {
      tmpx <- range(which(ordered == unique(ordered)[k]) - 0.5) + c(-0.25, 0.25)
      axis(
        side = 1, at = tmpx, labels = NA, line = 5,
        col = cluster_colours[k]
      )
      axis(
        side = 1, at = mean(tmpx),
        labels = paste0("Cluster ", k),
        line = 5, tick = FALSE, cex.axis = 1.25,
        col.axis = darken(cluster_colours[k], amount = 0.5)
      )
    }
    ordered <- Clusters(stab)[colnames(mat)]
    for (k in 1:nc_star) {
      tmpx <- range(nrow(mat) - which(ordered == unique(ordered)[k]) + 0.5) + c(-0.25, 0.25)
      axis(
        side = 2, at = tmpx, labels = NA, line = 5,
        col = cluster_colours[k]
      )
      axis(
        side = 2, at = mean(tmpx),
        labels = paste0("Cluster ", k),
        line = 5, tick = FALSE, cex.axis = 1.25,
        col.axis = darken(cluster_colours[k], amount = 0.5)
      )
    }
    abline(
      v = which(!duplicated(Clusters(stab)[myhclust_consensus$order]))[-1] - 1,
      lwd = 2, col = "blue"
    )
    abline(
      h = nrow(mat) - which(!duplicated(Clusters(stab)[myhclust_consensus$order]))[-1] + 1,
      lwd = 2, col = "blue"
    )
    dev.off()
  }
}


# # Ordered median of median weights
# stab=stab_w
# myclusters=Clusters(stab)
# ranked_list=rank(-apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median))
# WeightBoxplot(stab, at=ranked_list)
# myclusters=Clusters(stab)
# N=20
# {pdf("Figures/Boxplots_top_genes_lung.pdf",
#      width=20, height=5)
#   par(mfrow=c(1, 10), mar=c(5,3,5,1))
#   mycolours=lighten(c("navy", "darkolivegreen", "darkred", "darkmagenta"), amount=0.5)
#   for (k in 1:N){
#     tmpx=split(x[,names(sort(apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median), decreasing = TRUE))[k]],
#                f=paste0("C", myclusters))
#     boxplot(tmpx, range = 0, boxwex=0.35, frame=FALSE,
#             col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours,
#             medcol = darken(mycolours, amount = 0.4),
#             main=names(ranked_list)[k])
#   }
#   dev.off()}
# 
# {pdf("Figures/Boxplots_bottom_genes_lung.pdf",
#      width=20, height=5)
#   par(mfrow=c(1, 10), mar=c(5,3,5,1))
#   mycolours=lighten(c("navy", "darkolivegreen", "darkred", "darkmagenta"), amount=0.5)
#   for (k in 1:N){
#     tmpx=split(x[,names(sort(apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median), decreasing = TRUE))[ncol(stab$Beta)-k+1]],
#                f=paste0("C", myclusters))
#     boxplot(tmpx, range = 0, boxwex=0.35, frame=FALSE,
#             col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours,
#             medcol = darken(mycolours, amount = 0.4),
#             main=names(ranked_list)[k])
#   }
#   dev.off()}



# stab=stab_w
# WeightBoxplot(stability=stab)
# colnames(stab$Beta)[1:5]
#
# annot=data.frame(id=colnames(stab$Beta),
#                  chr=as.character(mget(colnames(stab$Beta), env=hgu95aCHR, ifnotfound=NA)),
#                  chr=as.character(mget(colnames(stab$Beta), env=hgu95aCHRLOC, ifnotfound=NA)),
#                  gene=as.character(mget(colnames(stab$Beta), env=hgu95aSYMBOL, ifnotfound=NA)))
#
# mget(colnames(stab$Beta), env=hgu95aPATH, ifnotfound=NA)
# kegg_pathways=do.call(c, mget(colnames(stab$Beta), env=hgu95aPATH, ifnotfound=NA))
# kegg_pathways=kegg_pathways[!is.na(kegg_pathways)]
# kegg_pathways=unique(kegg_pathways)
#
# disease_class_gad=fread("Data/Lung_cancer/chart_617FC072E0CF1659630322320.txt",
#                         sep="\t", data.table = FALSE)
# disease_class_gad$Term
# disease_class_gad$Genes
# cancer_genes=strsplit(disease_class_gad$Genes[which(disease_class_gad$Term=="CANCER")], split = ", ")[[1]]
#
# disease_gad=fread("Data/Lung_cancer/chart_617FC072E0CF1659630829212.txt",
#                   sep="\t", data.table = FALSE)
# lung_cancer_genes=unique(strsplit(disease_class_gad$Genes[which(disease_gad$Term=="lung cancer")], split = ", ")[[1]])
#
#
#
#
# rownames(annot)=annot$id
# gene_list=annot$gene
# gene_list=gene_list[!is.na(gene_list)]
# gene_list=unique(gene_list)
# write.table(cbind(gene_list), "Data/Lung_cancer/Gene_list.txt",
#             row.names=FALSE, col.names=FALSE, quote=FALSE)
#
# panther=data.table::fread("Data/Lung_cancer/pantherGeneList.txt",
#                           header = FALSE, sep="\t")
# table(panther[,6]) # Panther-GO Molecular Function
#
# mycolours=randomcoloR::distinctColorPalette(k=length(unique(annot$chr)))
# # colorRampPalette(c("navy","red"))(length(unique(annot$chr)))
# names(mycolours)=unique(annot$chr)
# x=apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median)
# # x=x[sort.list(annot[names(x),"chr"])]
# plot(x, col=mycolours[annot[names(x), "chr"]])
#
#
#
# x=apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median)
# plot(x, col=ifelse(annot[names(x),"gene"]%in%cancer_genes, yes="red", no="grey"))
#
# x=apply(stab$Beta[ArgmaxId(stab)[1],,], 1, median)
# plot(x, col=ifelse(annot[names(x),"gene"]%in%lung_cancer_genes, yes="red", no="grey"))
#
#
# boxplot(list(x[which(annot[names(x),"gene"]%in%lung_cancer_genes)],
#              x[which(!annot[names(x),"gene"]%in%lung_cancer_genes)]))
#
# boxplot(list(x[which(annot[names(x),"gene"]%in%cancer_genes)],
#              x[which(!annot[names(x),"gene"]%in%cancer_genes)]))
