rm(list = ls())

library(fake)

# Consensus score
f <- expression(sqrt(n_w + n_b) * (x_w / n_w - x_b / n_b) * sqrt(n_w * n_b) / sqrt((x_w + x_b) * (n_w + n_b - x_w - x_b)))

# Parameters
n_w <- 10
n_b <- 20

fmat <- matrix(NA, nrow = n_b + 1, ncol = n_w + 1)
for (x_w in seq(0, n_w)) {
  for (x_b in seq(0, n_b)) {
    fmat[x_b + 1, x_w + 1] <- eval(f)
  }
}
rownames(fmat) <- 0:n_b
colnames(fmat) <- 0:n_w

pdf("Figures/Consensus_score_heatmap.pdf")
par(mar = rep(5, 4))
fake::Heatmap(fmat, xlas = 1)
dev.off()
