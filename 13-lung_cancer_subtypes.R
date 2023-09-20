rm(list = ls())

library(sharp)
library(data.table)
library(readxl)
library(dendextend)
library(plotrix)
library(colorspace)

# Parameters
tau <- 0.5
noit <- 20
niter <- 10
max_nc <- 20
distance <- "euclidian"
linkage <- "complete"

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
myhclust <- hclust(d = dist(scale(x), method = distance), method = linkage)
dend <- as.dendrogram(myhclust)
labels_colors(dend) <- subtypes[labels(dend)]

# Consensus hierarchical clustering (unweighted)
stab_unw <- Clustering(
  xdata = x,
  nc = 1:max_nc,
  distance = distance,
  implementation = HierarchicalClustering,
  linkage = linkage,
  tau = tau
)
plot(stab_unw, theta_star = subtypes)
saveRDS(stab_unw, paste0("Results/Application/Consensus_unweighted_hclust_", distance, "_", tau, ".rds"))

# Consensus hierarchical clustering (weighted)
system.time({
  stab_w <- Clustering(
    xdata = x,
    nc = 1:max_nc,
    distance = distance,
    implementation = HierarchicalClustering,
    linkage = linkage,
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
cluster_colours <- lighten(c("tomato", "dodgerblue", "seagreen4", "tan", "darkmagenta"), amount = 0)

# Dendrogram
{
  pdf(paste0("Figures/Dendrogram_hierarchical_lung_subtypes_", linkage, ".pdf"),
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
    pdf(paste0("Figures/Calibration_", type, "_lung_subtypes_", linkage, ".pdf"),
      width = 12, height = 12
    )
    par(mar = rep(9, 4))
    CalibrationPlot(stab, xlab = "G", cex.lab = 2, cex.legend = 1.5)
    dev.off()
  }

  # Consensus matrix
  {
    pdf(paste0("Figures/Consensus_", type, "_lung_subtypes_", linkage, ".pdf"),
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
    nc=Argmax(stab)[1]
    for (k in 1:nc) {
      tmpx <- range(which(ordered == unique(ordered)[k]) - 0.5) + c(-0.25, 0.25)
      axis(
        side = 1, at = tmpx, labels = NA, line = 5,
        col = cluster_colours[k]
      )
      axis(
        side = 1, at = mean(tmpx),
        labels = paste0("Cluster ", k),
        line = 5, tick = FALSE, cex.axis = 1.1,
        col.axis = darken(cluster_colours[k], amount = 0.5)
      )
    }
    ordered <- Clusters(stab)[colnames(mat)]
    for (k in 1:nc) {
      tmpx <- range(nrow(mat) - which(ordered == unique(ordered)[k]) + 0.5) + c(-0.25, 0.25)
      axis(
        side = 2, at = tmpx, labels = NA, line = 5,
        col = cluster_colours[k]
      )
      axis(
        side = 2, at = mean(tmpx),
        labels = paste0("Cluster ", k),
        line = 5, tick = FALSE, cex.axis = 1.1,
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
