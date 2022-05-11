rm(list = ls())
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
for (name in r) eval(parse(text = paste0(name, "<-sharp:::", name)))

# Loading all additional functions
myfunctions <- list.files("Scripts/Functions/")
myfunctions <- myfunctions[myfunctions != "Former"]
for (k in 1:length(myfunctions)) {
  source(paste0("Scripts/Functions/", myfunctions[k]))
}

# Simulation of data with clusters
set.seed(1)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 3,
  nu_xc = 1,
  ev_xc = 0.8
)

# Making figure
mycolours <- darken(c("navy", "red", "gold"), amount = 0.2)
{
  pdf("Figures/Visualisation_simulation_3_variables.pdf",
    width = 12, height = 4
  )
  par(mfrow = c(1, 3), mar = c(5, 5, 1, 1))
  tmp <- matrix(c(1, 2, 1, 3, 2, 3), byrow = TRUE, ncol = 2)
  for (k in 1:nrow(tmp)) {
    xcomp <- tmp[k, 1]
    ycomp <- tmp[k, 2]
    plot(NA,
      xlim = range(simul$data[, 1:max(tmp)]),
      ylim = range(simul$data[, 1:max(tmp)]),
      las = 1, cex.lab = 1.5, cex.lab = 1.5,
      xlab = paste0("Variable ", xcomp),
      ylab = paste0("Variable ", ycomp),
    )
    abline(v = axTicks(1), lty = 3, col = "grey")
    abline(h = axTicks(2), lty = 3, col = "grey")
    text(simul$data[, xcomp], simul$data[, ycomp],
      pch = 19, cex = 0.7, las = 1,
      labels = gsub("obs", "", rownames(simul$data)),
      col = mycolours[simul$theta]
    )
  }
  dev.off()
}


# Simulation of data with clusters
set.seed(1)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 10,
  nu_xc = 1,
  ev_xc = 0.5
)

# Making figure
{
  pdf("Figures/Visualisation_simulation.pdf",
    height = 10, width = 11
  )
  # Heatmap of pairwise distances
  layout(
    mat = matrix(c(1, 1, 1, 2, 3, 4),
      nrow = 2, byrow = TRUE
    ),
    heights = c(2, 1)
  )
  par(mar = c(5, 19, 1, 19))
  x <- simul$data
  Heatmap(as.matrix(dist(x)))

  # Principal Component Analysis
  par(mar = c(5, 4.5, 2, 4.5))
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
  dev.off()
}

# Simulation of data with clusters
{
  pdf("Figures/Visualisation_simulation_different_p.pdf",
    height = 10, width = 15
  )
  par(mfrow = c(3, 4))
  for (p in c(5, 10, 30)) {
    set.seed(1)
    n <- c(20, 50, 30)
    simul <- SimulateClustering(
      n = n,
      pk = p,
      nu_xc = 1,
      ev_xc = 0.5
    )

    # Heatmap of pairwise distances
    par(mar = c(5, 5, 1.5, 4))
    x <- simul$data
    Heatmap(as.matrix(dist(x)))

    # Principal Component Analysis
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
  }
  dev.off()
}


# Simulation of data with clusters
mycolours <- darken(c("tan", "royalblue", "gold", "tomato", "darkolivegreen"), amount = 0.2)
{
  pdf("Figures/Visualisation_simulation_different_not_even.pdf",
    height = 10, width = 15
  )
  layout(
    mat = matrix(1:6, byrow = TRUE, nrow = 3),
    widths = c(1, 3)
  )
  for (seed in 1) {
    set.seed(seed)
    n <- c(20, 50, 30, 10, 40)
    p <- 5
    simul <- SimulateClustering(
      n = n,
      pk = p,
      nu_xc = 1,
      ev_xc = 0.8
    )

    # Heatmap of pairwise distances
    par(mar = c(5, 5, 1.5, 4))
    x <- simul$data
    Heatmap(as.matrix(dist(x)))

    # Silhouette plot
    mysilhouette <- silhouette(x = simul$theta, dist = dist(simul$data))
    plot(mysilhouette[, 3],
      col = mycolours[mysilhouette[, 1]],
      las = 1, bty = "n",
      # ylim = c(min_silhouette, max_silhouette),
      type = "h", lwd = 2, lend = 1,
      xaxt = "n", xlab = "",
      ylab = "Silhouette coefficient",
      cex.lab = 1.5,
      panel.first = abline(h = 0, lty = 2)
    )
  }
  dev.off()
}
