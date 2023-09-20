rm(list = ls())

library(colorspace)

# Simulation parameters
method <- "hclust"
metric <- "ari"
n_lambda <- 10
ylim <- c(-0.05, 1.05)

if (method %in% c("cosa_hclust", "sparcl_hclust")) {
  mycolours <- lighten(
    c(
      "tan", "darkorange",
      "tan",
      "darkred"
    ),
    amount = 0.3
  )
}

dataset_list <- seq(1, 5)

dir.create(paste0("Figures/Real_data_consensus_", method), showWarnings = FALSE)
for (dataset_id in dataset_list) {
  print(dataset_id)

  # Loading the results
  performances <- readRDS(paste0("Results/HPC/Real_data_consensus_unweighted_", method, "/Performances_", dataset_id, ".rds"))

  # Defining the algorithm names
  algo_names <- c("Single run", "Consensus")

  # Calibration name
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
    "Delta",
    "PAC",
    "PINS discrepancy",
    "RCSI (PAC)",
    "RCSI (entropy)",
    "Consensus score"
  )
  full_names <- paste0("'", full_names, " (G = ", performances[, "G"], ")'")
  full_names[c(1, 6)] <- paste0("'G* = ", performances[1, "G"], "'")
  names(full_names) <- c(
    "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap",
    "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db",
    "delta", "pac", "pins_discrepancy", "rcsi_pac", "rcsi_entropy", "consensus_score"
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
      darken("maroon4", amount = 0.3), "maroon4", lighten("maroon4", amount = 0.3),
      "navy", lighten("navy", amount = 0.3),
      "darkred"
    ),
    amount = 0.3
  )
  zseq <- c(0.5, 5.5, 15.5)
  id_list <- c(1, 6, 15)

  pdf(paste0("Figures/Real_data_consensus_", method, "/Plot_unweighted_", metric, "_", dataset_id, ".pdf"),
    width = 10, height = 7
  )
  par(mar = c(11, 5, 5, 6))
  plot(as.numeric(performances[, metric]),
    xlab = "", ylab = toupper(metric),
    las = 1, cex.lab = 1.5,
    bty = "n", xaxt = "n",
    ylim = ylim,
    xlim = c(0.5, nrow(performances) + 0.5),
    type = "h", lend = 1, lwd = 20,
    col = mycolours
  )
  abline(h = axTicks(2), lty = 3, col = "grey")
  xseq <- (1:length(full_names))
  abline(v = zseq, lty = 2, col = "black")
  for (id in id_list) {
    abline(h = performances[id, metric], col = mycolours[id], lty = 2)
  }
  axis(side = 1, at = xseq, labels = NA)
  for (i in 1:length(full_names)) {
    axis(
      side = 1, at = xseq[i],
      labels = eval(parse(text = paste0("expression(", (full_names)[i], ")"))),
      las = 2, tick = FALSE
    )
  }
  for (i in 1:length(algo_names)) {
    axis(
      side = 3, at = c(zseq[i], zseq[i + 1]),
      labels = NA,
      cex.axis = 1, tick = TRUE
    )
    axis(
      side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean)[i],
      labels = algo_names[i],
      cex.axis = 1, tick = FALSE
    )
  }
  points(as.numeric(performances[, metric]),
    type = "h", lend = 1, lwd = 20,
    col = mycolours
  )
  dev.off()

  pdf(paste0("Figures/Real_data_consensus_", method, "/Plot_unweighted_calibration_", dataset_id, ".pdf"),
    width = 10, height = 7
  )
  # Loading the results
  stab <- readRDS(paste0("Results/HPC/Real_data_consensus_unweighted_", method, "/Stability_unweighted_", dataset_id, ".rds"))
  par(mar = c(11, 5, 5, 25))
  CalibrationPlot(stab)
  dev.off()
}
