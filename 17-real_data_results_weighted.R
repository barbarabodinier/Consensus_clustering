rm(list = ls())

library(sharp)
library(colorspace)

# Simulation parameters
# method <- "cosa_hclust"
# method <- "sparcl_kmeans"
# method <- "sparcl_hclust"
method <- "impacc_hclust"
metric <- "ari"
n_lambda <- 10
ylim <- c(-0.05, 1.05)

if (method %in% c("cosa_hclust", "sparcl_hclust", "sparcl_kmeans")) {
  mycolours <- lighten(
    c(
      "tan", "darkorange",
      lighten("maroon4", amount = seq(0.1, 0.7, length.out = 10)),
      "tan",
      "darkred",
      lighten("navy", amount = seq(0.1, 0.7, length.out = 10)),
      "darkred"
    ),
    amount = 0.3
  )
}
if (method == "impacc_hclust") {
  mycolours <- lighten(
    c(
      "tan",
      "darkorange",
      "tan",
      lighten("navy", amount = seq(0.1, 0.7, length.out = 9))
    ),
    amount = 0.3
  )
}
dataset_list <- 6:10

dir.create(paste0("Figures/Real_data_consensus_", method), showWarnings = FALSE)
for (dataset_id in dataset_list) {
  print(dataset_id)

  # Loading the results
  performances <- readRDS(paste0("Results/HPC/Real_data_consensus_", method, "/Performances_", dataset_id, ".rds"))
  if (method=="impacc_hclust"){
    performances=performances[-which(rownames(performances)=="consensus"),]
  }

  # Defining the algorithm names
  if (method == "cosa_hclust") {
    algo_names <- c("Unweighted", "COSA", "Unweighted", "COSA")
  }
  if (method %in% c("sparcl_hclust", "sparcl_kmeans")) {
    algo_names <- c("Unweighted", "sparcl", "Unweighted", "sparcl")
  }
  if (method == "impacc_hclust") {
    algo_names <- c("Unweighted", "Unweighted", "IMPACC")
  }
  if (method == "sparse_gmm") {
    algo_names <- c("Unweighted", "clustvarsel", "VarSelLCM")
  }

  if (method == "cosa_hclust") {
    lambda_list <- formatC(as.numeric(performances[, "lambda"]), format = "f", digits = 2)
    full_names <- c(
      "G*", 
      "GAP",
      paste0("lambda[", 1:10, "]*'=", lambda_list[3:12], "'"),
      "G*",
      "sharp score",
      paste0("lambda[", 1:10, "]*'=", lambda_list[15:24], "'"),
      "sharp score"
    )
    names(full_names) <- c(
      "single_run_star", "single_run_gap",
      paste0("single_run_cosa_star_", 1:n_lambda),
      "consensus_star", "consensus",
      paste0("consensus_cosa_star_", 1:n_lambda),
      "consensus_cosa"
    )
    ids_g=c(2, 14, length(full_names))
    full_names[ids_g] <- paste0("'", full_names[ids_g], " (G = ", performances[ids_g, "G"], ")'")
    full_names[c(1, 13)] <- paste0("'G* = ", performances[1, "G"], "'")
  }
  if (method %in% c("sparcl_hclust", "sparcl_kmeans")) {
    lambda_list <- formatC(as.numeric(performances[, "lambda"]), format = "f", digits = 2)
    full_names <- c(
      "G*", "GAP",
      paste0("lambda[", 1:10, "]*'=", lambda_list[3:12], "'"),
      "G*",
      "sharp score",
      paste0("lambda[", 1:10, "]*'=", lambda_list[15:24], "'"),
      "sharp score"
    )
    names(full_names) <- c(
      "single_run_star", "single_run_gap",
      paste0("single_run_sparcl_star_", 1:n_lambda),
      "consensus_star", "consensus",
      paste0("consensus_sparcl_star_", 1:n_lambda),
      "consensus_sparcl"
    )
    ids_g=c(2, 14, length(full_names))
    full_names[ids_g] <- paste0("'", full_names[ids_g], " (G = ", performances[ids_g, "G"], ")'")
    full_names[c(1, 13)] <- paste0("'G* = ", performances[1, "G"], "'")
  }
  if (method == "impacc_hclust") {
    lambda_list <- formatC(as.numeric(performances[, "lambda"]), format = "f", digits = 2)
    full_names <- c(
      "'G*'", "GAP",
      "'G*'",
      seq(0.1, 0.9, by = 0.1)
    )
    names(full_names) <- c(
      "single_run_star", "single_run_gap",
      "consensus_star", 
      paste0("impacc_", seq(1, 9))
    )
    ids_g=c(2)
    full_names[ids_g] <- paste0("'", full_names[ids_g], " (G = ", performances[ids_g, "G"], ")'")
    full_names[c(1, 3)] <- paste0("'G* = ", performances[1, "G"], "'")
  }

  pdf(paste0("Figures/Real_data_consensus_", method, "/Plot_", metric, "_", dataset_id, ".pdf"),
    width = 14, height = 7
  )
  par(mar = c(11, 5, 5, 6))
  plot(as.numeric(performances[, metric]),
    xlab = "", ylab = toupper(metric),
    las = 1, cex.lab = 1.5,
    bty = "n", xaxt = "n", ylim = ylim,
    type = "h", lend = 1, lwd = 20,
    col = mycolours
  )
  abline(h = axTicks(2), lty = 3, col = "grey")
  xseq <- (1:length(full_names))
  if (method %in% c("cosa_hclust", "sparcl_hclust", "sparcl_kmeans")) {
    zseq <- c(0.5, 2.5, 12.5, 14.5, length(xseq) + 0.5)
    abline(v = zseq, lty = 2, col = "black")
  }
  if (method %in% c("impacc_hclust")) {
    zseq <- c(0.5, 2.5, 3.5, length(xseq) + 0.5)
    abline(v = zseq, lty = 2, col = "black")
  }
  if (method %in% c("cosa_hclust", "sparcl_hclust", "sparcl_kmeans")) {
    id_list <- c(1, 13, 14, 25)
    for (id in id_list) {
      abline(h = performances[id, metric], col = mycolours[id], lty = 2)
    }
  }
  if (method %in% c("impacc_hclust")) {
    id_list <- c(1, 3, 4)
    for (id in id_list) {
      abline(h = performances[id, metric], col = mycolours[id], lty = 2)
    }
  }
  axis(side = 1, at = xseq, labels = NA)
  for (i in 1:length(full_names)) {
    axis(
      side = 1, at = xseq[i],
      labels = eval(parse(text = paste0("expression(", (full_names)[i], ")"))),
      las = 2, tick = FALSE
    )
  }
  if (method %in% c("cosa_hclust", "sparcl_hclust", "sparcl_kmeans")) {
    myline <- 8
    axis(side = 1, at = xseq[c(3, 12)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(3, 12)]), 
         labels = paste0("G* = ", performances[1, "G"]), 
         line = myline - 0.5, tick = FALSE)
    axis(side = 1, at = xseq[c(15, 24)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(15, 24)]), 
         labels = paste0("G* = ", performances[1, "G"]), 
         line = myline - 0.5, tick = FALSE)
    axis(side = 3, at = zseq, labels = NA)
  }
  if (method %in% c("impacc_hclust")) {
    myline <- 3
    axis(side = 1, at = xseq[c(4, 12)] + c(-0.25, 0.25), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(4, 12)]), 
         labels = paste0("G* = ", performances[1, "G"]), 
         line = myline - 0.5, tick = FALSE)
    # axis(side = 1, at = xseq[c(15, 24)] + c(-0.5, 0.5), labels = NA, line = myline)
    # axis(side = 1, at = mean(xseq[c(15, 24)]), 
    #      labels = paste0("G* = ", performances[1, "G"]), 
    #      line = myline - 0.5, tick = FALSE)
    # axis(side = 3, at = zseq, labels = NA)
  }
  for (i in 1:length(algo_names)) {
    axis(
      side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean)[i],
      labels = algo_names[i],
      cex.axis = 1, tick = FALSE
    )
  }
  if (method != "sparse_gmm") {
    if (method == "impacc_hclust") {
      axis(side = 3, at = zseq[c(1, 2)] + c(0.1, -0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = mean(zseq[c(1, 2)]),
        labels = "Single run",
        cex.axis = 1, tick = FALSE, line = 2
      )
      zseq2 <- zseq[c(1, 2, length(zseq))]
      axis(side = 3, at = zseq2[2:3] + c(0.1, -0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = apply(rbind(zseq2[-1], zseq2[-length(zseq2)]), 2, mean)[2],
        labels = "Consensus clustering",
        cex.axis = 1, tick = FALSE, line = 2
      )
    } else {
      zseq2 <- zseq[c(1, 2, 3)]
      axis(side = 3, at = zseq2[c(1, 3)] + c(0.1, -0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = mean(zseq2[c(1, 3)]),
        labels = "Single run",
        cex.axis = 1, tick = FALSE, line = 2
      )
      zseq2 <- zseq[c(2, 3, length(zseq))]
      axis(side = 3, at = zseq2[2:3] + c(0.1, -0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = apply(rbind(zseq2[-1], zseq2[-length(zseq2)]), 2, mean)[2],
        labels = "Consensus clustering",
        cex.axis = 1, tick = FALSE, line = 2
      )
    }
  }
  points(as.numeric(performances[, metric]),
         type = "h", lend = 1, lwd = 20,
         col = mycolours
  )
  dev.off()

  if (method == "cosa_hclust") {
    pdf(paste0("Figures/Real_data_consensus_", method, "/Plot_weighted_calibration_", dataset_id, ".pdf"),
      width = 14, height = 7
    )
    # Loading the results
    stab <- readRDS(paste0("Results/HPC/Real_data_consensus_", method, "/Stability_weighted_", dataset_id, ".rds"))
    par(mar = c(11, 5, 5, 45))
    CalibrationPlot(stab, 
                    ylab="sharp score",
                    legend=(dataset_id=="6"),
                    ncol=2, 
                    cex.legend = 0.8)
    dev.off()
  }
}
