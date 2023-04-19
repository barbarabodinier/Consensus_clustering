rm(list = ls())

library(fake)
library(rCOSA)
library(sharp)
library(colorspace)

dir.create("Figures/Graphical_abstract")

# Simulation of data with clusters
set.seed(0)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 10,
  theta_xc = c(rep(1, 10), rep(0, 10)),
  nu_xc = 1,
  ev_xc = 0.7
)
x <- simul$data

# Consensus clustering
stab <- Clustering(
  xdata = x,
  implementation = HierarchicalClustering,
  nc = 1:10,
  Lambda = c(10, 1, 0.1)
)
CalibrationPlot(stab)
dev.off()

# Calculating distance matrices
Lambda <- c(10, 1, 0.1)
for (k in 1:3) {
  tmpdist <- rCOSA::cosa2(
    X = as.data.frame(x),
    lX = rep(1, ncol(x)),
    lambda = Lambda[k],
    stand = 0,
    pwr = 2
  )
  print(range(tmpdist$D, na.rm = TRUE))
  assign(paste0("dist_", k), tmpdist$D)
}

# Saving figures
for (k in 1:3) {
  for (id_j in 1:3) {
    # Subsampling
    set.seed(id_j)
    s <- Resample(data = x)
    mydist_full <- as.matrix(eval(parse(text = paste0("dist_", k))))
    mydist_full[s, s] <- 0
    rownames(mydist_full) <- colnames(mydist_full) <- paste0("obs", 1:ncol(mydist_full))
    print(range(mydist_full, na.rm = TRUE))

    # Distance matrix
    par(mar = rep(5, 4), bg = NA)
    Heatmap(mydist_full,
      legend_range = c(0, 6),
      legend = ifelse((id_j == 3) & (k == 3), yes = TRUE, no = FALSE)
    )
    dev.copy(png, paste0("Figures/Graphical_abstract/Distance_", id_j, "_", k, ".png"))
    dev.off()

    # Co-membership matrix
    for (id_i in 2:4) {
      mydist <- mydist_full[-s, -s]
      myhclust <- hclust(as.dist(mydist), method = "complete")
      members <- CoMembership(cutree(myhclust, k = id_i))
      members_full <- matrix(-1, ncol = nrow(x), nrow = nrow(x))
      rownames(members_full) <- colnames(members_full) <- rownames(x)
      for (i in 1:nrow(members_full)) {
        row_id <- rownames(members_full)[i]
        if (row_id %in% rownames(members)) {
          for (j in 1:ncol(members_full)) {
            col_id <- colnames(members_full)[j]
            if (col_id %in% colnames(members)) {
              members_full[row_id, col_id] <- members[row_id, col_id]
            }
          }
        }
      }
      Heatmap(members_full, col = c("white", "royalblue", "darkred"), legend = FALSE)
      dev.copy(png, paste0("Figures/Graphical_abstract/Comembership_matrix_", id_i, "_", id_j, "_", k, ".png"))
      dev.off()
    }
  }

  # Consensus matrix
  for (id_i in 2:4) {
    id <- which((stab$nc == id_i) & (stab$Lambda == Lambda[k]))
    Heatmap(ConsensusMatrix(stab, argmax_id = id),
      legend = ifelse(k == 3, yes = TRUE, no = FALSE)
    )
    dev.copy(png, paste0("Figures/Graphical_abstract/Consensus_matrix_", id_i, "_", k, ".png"))
    dev.off()
  }
}

pdf("Figures/Graphical_abstract/Calibration_plot.pdf",
  width = 14, height = 5
)
par(mar = c(5, 5, 5, 1))
mycolours <- c("red", "blue", "forestgreen")
plot(NULL,
  xlim = c(2, 10),
  ylim = c(100, 350),
  xlab = "", ylab = "",
  xaxt = "n", cex.axis = 1,
  panel.first = abline(v = 2:10, lty = 2)
)
for (k in 1:3) {
  ids <- which(stab$Lambda == Lambda[k])
  lines(stab$nc[ids], stab$Sc[ids],
    col = mycolours[k]
  )
}
points(stab$nc, stab$Sc,
  cex.lab = 1.5,
  pch = 19, cex = 3,
  col = adjustcolor(rep(mycolours, each = 10),
    alpha.f = 0.6
  )
)
axis(side = 3, at = 2:10, las = 2, cex.axis = 2)
dev.off()
