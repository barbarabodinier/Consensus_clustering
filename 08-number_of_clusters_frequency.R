rm(list = ls())

library(colorspace)
library(openxlsx)

# Simulation parameters
method <- "hclust"
simul_study_id <- "1"

# Defining the names and colours
full_names <- c(
  "G*",
  "Silhouette",
  "CH",
  "DB",
  "GAP statistic",
  "G*",
  "Silhouette",
  "CH",
  "DB",
  "GAP statistic",
  "Delta",
  "PAC",
  "PINS discrepancy",
  "RCSI (PAC)",
  "RCSI (entropy)",
  "sharp score"
)
names(full_names) <- c(
  "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap",
  "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db", "consensus_gap",
  "delta", "pac", "pins_discrepancy", "rcsi_pac", "rcsi_entropy", "consensus_score"
)
clustering_names <- c(
  rep("SR", 5),
  rep("C", 11)
)
mycolours <- lighten(
  c(
    "tan",
    "darkgreen",
    lighten("darkgreen", amount = 0.3),
    lighten("darkgreen", amount = 0.5),
    "darkorange",
    "tan",
    "darkgreen",
    lighten("darkgreen", amount = 0.3),
    lighten("darkgreen", amount = 0.5),
    "darkorange",
    darken("maroon4", amount = 0.3), "maroon4", lighten("maroon4", amount = 0.3),
    "navy", lighten("navy", amount = 0.3),
    "darkred"
  ),
  amount = 0.3
)

# Loop over the simulation IDs
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
  performances[, "G", ]

  # Calculating the counts per number of clusters
  g_counts <- NULL
  for (method_id in 1:nrow(performances)) {
    perf_tmp <- factor(performances[method_id, "G", ], levels = 2:20)
    table_tmp <- table(perf_tmp)
    table_tmp["10"] <- sum(table_tmp[as.character(10:20)])
    table_tmp <- table_tmp[as.character(2:10)]
    g_counts <- rbind(g_counts, table_tmp)
  }
  g_counts <- rbind(g_counts, NA)

  pdf(paste0("Figures/Simulations_consensus_", method, "/Barplot_G_", simul_study_id, "_", simul_id, ".pdf"),
    width = 14, height = 5
  )
  mylwd <- 5
  par(mar = c(5, 5, 1, 1))
  plot(as.vector(g_counts),
    type = "h",
    col = c(mycolours, "white"),
    xaxt = "n", las = 1, cex.lab = 1.5,
    xlab = "Number of clusters (G)",
    ylab = "Count",
    lend = 1, lwd = mylwd
  )
  xseq <- seq(0, length(g_counts) + 1, by = nrow(g_counts))
  abline(v = xseq, lty = 2, col = "grey")
  axis(1,
    at = apply(rbind(xseq[-1], xseq[-length(xseq)]), 2, mean),
    labels = c(2:4, "G* = 5", 6:9, "10+")
  )
  if (simul_id == 1) {
    legend("topright",
      ncol = 4, bg = "white", cex = 0.9,
      legend = paste0(clustering_names, ": ", full_names),
      col = mycolours, lty = 1, lwd = mylwd
    )
  }
  dev.off()
}
