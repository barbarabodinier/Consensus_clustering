rm(list = ls())
setwd("~/Dropbox/Consensus_clustering")

library(openxlsx)

# Simulation parameters
method <- "weighted_hierarchical"
n_lambda <- 10
for (simul_study_id in "4") {
  print(paste0("Simulation study ", simul_study_id))

  # Template design
  mycolours <- lighten(c(
    "tan",
    "darkgreen", lighten("darkgreen", amount = 0.3),
    "orange", "darkred",
    rep(lighten("maroon4", amount = 0.2), 10), "darkred",
    rep(lighten("navy", amount = 0.2), 10), "darkred"
  ),
  amount = 0.3
  )
  dimensionality <- c("", "", "")

  # Saving table
  continuous_metrics <- c("rand", "ari", "jaccard")
  integer_metrics <- c("G", "q", "time")
  binary_metrics <- "signif"
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", 1, "_merged.rds"))
  lambda_list <- formatC(as.numeric(performances[, "lambda", 1]), format = "f", digits = 2)
  full_names <- c(
    "'G*'",
    "'Silhouette score'",
    "'GAP statistic'",
    "'G*'",
    "'Consensus score'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[6:15], "'"),
    "'Consensus score'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[17:26], "'"),
    "'Consensus score'"
  )
  names(full_names) <- c(
    "hclust_star", "hclust_silhouette", "hclust_gap",
    "unweighted_star", "unweighted",
    paste0("sparcl_star_", 1:n_lambda), "sparcl",
    paste0("cosa_star_", 1:n_lambda), "cosa"
  )
  algo_names <- c("Hierarchical", "Hierarchical", "sparcl", "COSA")
  n_simul <- length(list.files(path = paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id))) / 2
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
      IQR(as.numeric(x), na.rm = TRUE)
    })
    myiqr <- tmpiqr
    myiqr[, c(continuous_metrics)] <- formatC(tmpiqr[, c(continuous_metrics)], format = "f", digits = 3)
    myiqr[, c(integer_metrics)] <- formatC(tmpiqr[, c(integer_metrics)], format = "f", digits = 0)
    myiqr <- myiqr[, c(integer_metrics[-length(integer_metrics)], continuous_metrics, "time")]

    tmptable <- matrix(paste0(mymedian[, ], " [", myiqr, "]"), ncol = ncol(mymedian))
    colnames(tmptable) <- colnames(mymedian)

    assign(paste0("tmptable", simul_id), tmptable)
    mytable <- rbind(mytable, tmptable)
  }
  # mytable <- rbind(mytable1, mytable2, mytable3)
  mytable <- cbind(rep(full_names[rownames(performances)], n_simul), mytable)
  mytable <- cbind(rep("", nrow(mytable)), mytable)

  write.xlsx(as.data.frame(mytable), paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".xlsx"),
    rowNames = FALSE, colNames = TRUE, overwrite = TRUE
  )
  mytable[,2]=gsub("'", "", mytable[,2])
  ids=grep("lambda", mytable[,2])
  mytable[ids,2]=paste0("$", gsub("]*", "", gsub("lambda[", "\\lambda_", mytable[ids,2], fixed=TRUE), fixed=TRUE), "$")
  mytable[,4]=gsub("NA [NA]", "", mytable[,4], fixed=TRUE)
  write.table(mytable, paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".txt"),
    row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
  )

  # Saving figure
  for (metric in c("G", "ari")) {
    for (simul_id in 1:n_simul) {{ pdf(paste0("Figures/Simulations_consensus_", method, "/Boxplot_", metric, "_", simul_study_id, "_", simul_id, ".pdf"),
      width = 14, height = 7
    )
    par(mar = c(11, 5, 5, 6))
    performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
    mylist <- list()
    for (k in 1:nrow(performances)) {
      mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
      assign(paste0("median", k), median(as.numeric(performances[k, metric, ])))
    }
    tmpfullnames <- full_names
    tmpfullnames <- paste0(tmpfullnames, "*' (q=", apply(performances[, "q", ], 1, FUN = function(x) {
      formatC(median(as.numeric(x), na.rm = TRUE), format = "f", digits = 0)
    }), ")'")
    tmpfullnames <- gsub(" \\(q=NA\\)", "", tmpfullnames)
    names(tmpfullnames) <- names(full_names)
    xseq <- (1:length(full_names))
    if (metric == "ari") {
      ylim <- c(0, 1)
    } else {
      ylim <- c(0, 30)
    }
    myylab <- toupper(metric)
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, main = dimensionality[simul_id], cex.main = 1.5,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", ylim = ylim, frame = "F", boxwex = 0.25
    )
    abline(h = axTicks(2), lty = 3, col = "grey")
    zseq <- c(0.5, 3.5, 5.5, 16.5, length(xseq) + 0.5)
    abline(v = zseq, lty = 2, col = "black")
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, add = TRUE,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", frame = "F", boxwex = 0.25
    )
    for (id in c(1, 4, 5, 16, 27)) {
      abline(h = eval(parse(text = paste0("median", id))), col = darken(mycolours[id], amount = 0.4), lty = 2)
    }
    axis(side = 1, at = xseq, labels = NA)
    for (i in 1:nrow(performances)) {
      axis(
        side = 1, at = xseq[i],
        labels = eval(parse(text = paste0("expression(", (tmpfullnames[rownames(performances)])[i], ")"))), 
        las = 2, tick = FALSE
      )
    }
    axis(side = 1, at = xseq[c(6, 15)] + c(-0.5, 0.5), labels = NA, line = 8)
    axis(side = 1, at = mean(xseq[c(6, 15)]), labels = "G*", line = 8, tick = FALSE)
    axis(side = 1, at = xseq[c(17, 26)] + c(-0.5, 0.5), labels = NA, line = 8)
    axis(side = 1, at = mean(xseq[c(17, 26)]), labels = "G*", line = 8, tick = FALSE)
    axis(side = 3, at = zseq, labels = NA)
    for (i in 1:length(algo_names)) {
      axis(
        side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean)[i],
        labels = algo_names[i], 
        cex.axis = 1, tick = FALSE
      )
    }
    zseq2 <- zseq[c(1, 2, length(zseq))]
    axis(side = 3, at = zseq2[2:3], labels = NA, line = 2.5)
    axis(
      side = 3, at = apply(rbind(zseq2[-1], zseq2[-length(zseq2)]), 2, mean)[2],
      labels = "Consensus clustering", 
      cex.axis = 1, tick = FALSE, line = 2.5
    )
    dev.off() }}
  }
}