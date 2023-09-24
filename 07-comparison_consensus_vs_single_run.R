rm(list = ls())

library(colorspace)
library(openxlsx)

# Simulation parameters
method <- "hclust"
# method <- "pam"
# method <- "kmeans"
# method <- "gmm"

# Creating output folders
dir.create(paste0("Tables/Simulations_consensus_", method), showWarnings = FALSE)
dir.create(paste0("Figures/Simulations_consensus_", method), showWarnings = FALSE)

# Loop over the simulation study ID
for (simul_study_id in c("2")) {
  print(paste0("Simulation study ", simul_study_id))

  # Template design
  dimensionality <- c("", "", "")

  # Saving table
  continuous_metrics <- c("rand", "ari", "jaccard")
  integer_metrics <- c("G", "time")
  binary_metrics <- "signif"
  if (method=="gmm"){
    full_names <- c(
      "G*",
      "Silhouette",
      "CH",
      "DB",
      "GAP statistic",
      "BIC",
      "G*",
      "Silhouette",
      "CH",
      "DB",
      "Delta",
      "PAC",
      "PINS discrepancy",
      "Consensus score"
    )
    names(full_names) <- c(
      "single_run_star", "single_run_silhouette", "single_run_ch", "single_run_db", "single_run_gap", "single_run_bic",
      "consensus_star", "consensus_silhouette", "consensus_ch", "consensus_db",
      "delta", "pac", "pins_discrepancy", "consensus_score"
    )
    mycolours <- lighten(
      c(
        "tan",
        "darkgreen",
        lighten("darkgreen", amount = 0.3),
        lighten("darkgreen", amount = 0.5),
        "darkorange",
        lighten("darkorange", amount = 0.3),
        "tan",
        "darkgreen",
        lighten("darkgreen", amount = 0.3),
        lighten("darkgreen", amount = 0.5),
        darken("maroon4", amount = 0.3), "maroon4", lighten("maroon4", amount = 0.3),
        "darkred"
      ),
      amount = 0.3
    )
    zseq <- c(0.5, 6.5, 14.5)
    id_set=c(1, 6, 14)
  } else {
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
    zseq <- c(0.5, 5.5, 16.5)
    id_set=c(1, 6, 16)
  }
  if (method == "hclust") {
    algo_names <- c("Hierarchical", "Consensus hierarchical")
  }
  if (method == "kmeans") {
    algo_names <- c("K-means", "Consensus K-means")
  }
  if (method == "gmm") {
    algo_names <- c("GMM", "Consensus GMM")
  }
  if (method == "pam") {
    algo_names <- c("PAM", "Consensus PAM")
  }
  n_simul <- length(list.files(path = paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id)))
  mytable <- NULL
  for (simul_id in 1:n_simul) {
    performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))

    tmpmedian <- apply(performances, c(1, 2), FUN = function(x) {
      median(as.numeric(x))
    })
    mymedian <- tmpmedian
    mymedian[, c(continuous_metrics)] <- formatC(tmpmedian[, c(continuous_metrics)], format = "f", digits = 3)
    mymedian[, c(integer_metrics)] <- formatC(tmpmedian[, c(integer_metrics)], format = "f", digits = 0)
    mymedian <- mymedian[, c(integer_metrics[-length(integer_metrics)], continuous_metrics, "time")]

    tmpiqr <- apply(performances[, c(continuous_metrics, integer_metrics), ], c(1, 2), FUN = function(x) {
      IQR(as.numeric(x))
    })
    myiqr <- tmpiqr
    myiqr[, c(continuous_metrics)] <- formatC(tmpiqr[, c(continuous_metrics)], format = "f", digits = 3)
    myiqr[, c(integer_metrics)] <- formatC(tmpiqr[, c(integer_metrics)], format = "f", digits = 0)
    myiqr <- myiqr[, c(integer_metrics[-length(integer_metrics)], continuous_metrics, "time")]

    tmptable <- matrix(paste0(mymedian[, ], " [", myiqr, "]"), ncol = ncol(mymedian))
    colnames(tmptable) <- colnames(mymedian)

    mypercentage <- paste0(apply(performances[, binary_metrics, ], 1, sum) / dim(performances)[3] * 100, "%")
    mypercentage <- ifelse(mypercentage == "NA%", yes = "", no = mypercentage)

    tmptable <- cbind(tmptable[, -ncol(tmptable)], signif = mypercentage, tmptable[, ncol(tmptable), drop = FALSE])

    assign(paste0("tmptable", simul_id), tmptable)
    mytable <- rbind(mytable, tmptable)
  }
  # mytable <- rbind(mytable1, mytable2, mytable3)
  mytable <- cbind(rep(full_names[rownames(performances)], 3), mytable)
  mytable <- cbind(rep("", nrow(mytable)), mytable)

  write.xlsx(as.data.frame(mytable), paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".xlsx"),
    rowNames = FALSE, colNames = TRUE, overwrite = TRUE
  )
  write.table(mytable, paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".txt"),
    row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
  )

  # Saving figure
  for (metric in c("G", "ari")) {{ pdf(paste0("Figures/Simulations_consensus_", method, "/Boxplot_", metric, "_", simul_study_id, ".pdf"),
    width = 14, height = 4
  )
  par(mar = c(8, 5, 2, 6), mfrow = c(1, 3))
  for (simul_id in 1:n_simul) {
    performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
    mylist <- list()
    for (k in 1:nrow(performances)) {
      mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
      assign(paste0("median", k), median(as.numeric(performances[k, metric, ])))
    }
    xseq <- (1:length(full_names))
    if (metric == "ari") {
      ylim <- c(0, 1)
    } else {
      ylim <- c(0, 30)
    }
    myylab <- ifelse(simul_id == 1, yes = toupper(metric), no = "")
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, main = dimensionality[simul_id], cex.main = 1.5,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", ylim = ylim, frame = "F", boxwex = 0.35
    )
    abline(h = axTicks(2), lty = 3, col = "grey")
    abline(v = zseq, lty = 2, col = "black")
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, add = TRUE,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", frame = "F", boxwex = 0.35
    )
    for (id in id_set) {
      abline(h = eval(parse(text = paste0("median", id))), col = darken(mycolours[id], amount = 0.4), lty = 2)
    }
    axis(side = 1, at = xseq, labels = full_names[rownames(performances)], las = 2)
    axis(side = 3, at = zseq, labels = NA)
    axis(side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean), labels = algo_names, tick = FALSE)
  }
  dev.off() }}
}
