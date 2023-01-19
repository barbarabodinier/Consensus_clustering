rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

library(fake)
library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
for (simul_study_id in 1:3) {
  print(paste0("Simulation study ", simul_study_id))
  params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )

  # Heatmap of pairwise Euclidian distances
  {
    pdf(paste0("Figures/Heatmaps_examples_simul_study_", simul_study_id, ".pdf"),
      width = 14, height = 4
    )
    par(mfrow = c(1, nrow(params_list)), mar = c(5, 5, 2, 6))
    max_silhouette <- 0
    min_silhouette <- 0
    for (params_id in 1:nrow(params_list)) {
      # Extracting simulation parameters
      nc <- params_list[params_id, "nc"]
      equal_size <- params_list[params_id, "equal_size"]
      n_tot <- params_list[params_id, "n_tot"]
      p <- params_list[params_id, "p"]
      ev_xc <- params_list[params_id, "ev_xc"]
      nu_xc <- params_list[params_id, "nu_xc"]
      v_min <- params_list[params_id, "v_min"]
      v_max <- params_list[params_id, "v_max"]

      # Data simulation
      set.seed(1)
      if (equal_size) {
        n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
      } else {
        n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
      }
      pk <- round(rep(0.2, 5) * p)
      sigma <- SimulateCorrelation(
        pk = pk,
        nu_within = 1,
        nu_between = 0,
        v_within = c(v_min, v_max),
        v_between = 0,
        v_sign = -1,
        pd_strategy = "min_eigenvalue"
      )$sigma
      simul <- SimulateClustering(
        n = n,
        pk = pk,
        sigma = sigma,
        ev_xc = ev_xc,
        nu_xc = nu_xc,
        output_matrices = TRUE
      )
      simul$data <- scale(simul$data)

      # Heatmap
      Heatmap(as.matrix(dist(simul$data)))
      mysilhouette <- silhouette(x = simul$theta, dist = dist(simul$data))
      max_silhouette <- max(c(max_silhouette, mysilhouette[, 3]))
      min_silhouette <- min(c(min_silhouette, mysilhouette[, 3]))
    }
    dev.off()
  }

  # Silhouette plot of simulated clusters
  {
    pdf(paste0("Figures/Silhouette_examples_simul_study_", simul_study_id, ".pdf"),
      width = 14, height = 4
    )
    par(mfrow = c(1, nrow(params_list)), mar = c(8, 5, 2, 6))
    for (params_id in 1:nrow(params_list)) {
      # Extracting simulation parameters
      nc <- params_list[params_id, "nc"]
      n_tot <- params_list[params_id, "n_tot"]
      p <- params_list[params_id, "p"]
      ev_xc <- params_list[params_id, "ev_xc"]
      nu_xc <- params_list[params_id, "nu_xc"]
      v_min <- params_list[params_id, "v_min"]
      v_max <- params_list[params_id, "v_max"]

      # Data simulation
      set.seed(1)
      if (nc == 5) {
        n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
        mycolours <- darken(c("tan", "royalblue", "gold", "tomato", "darkolivegreen"), amount = 0.2)
      } else {
        n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
        mycolours <- darken(randomcoloR::distinctColorPalette(k = nc), amount = 0.2)
      }
      pk <- round(rep(0.2, 5) * p)
      sigma <- SimulateCorrelation(
        pk = pk,
        nu_within = 1,
        nu_between = 0,
        v_within = c(v_min, v_max),
        v_between = 0,
        v_sign = -1,
        pd_strategy = "min_eigenvalue"
      )$sigma
      simul <- SimulateClustering(
        n = n,
        pk = pk,
        sigma = sigma,
        ev_xc = ev_xc,
        nu_xc = nu_xc,
        output_matrices = TRUE
      )
      simul$data <- scale(simul$data)

      # Silhouette plot
      mysilhouette <- silhouette(x = simul$theta, dist = dist(simul$data))
      plot(mysilhouette[, 3],
        col = mycolours[mysilhouette[, 1]],
        las = 1, bty = "n",
        ylim = c(min_silhouette, max_silhouette),
        type = "h", lwd = 1, lend = 1,
        xaxt = "n", xlab = "",
        ylab = ifelse(params_id == 1, yes = "Silhouette coefficient", no = ""),
        cex.lab = 1.5,
        panel.first = abline(h = 0, lty = 2)
      )
    }
    dev.off()
  }
}
