rm(list = ls())

library(colorspace)
library(openxlsx)

# Simulation parameters
type="clustering"
# type="selection"
# method <- "cosa_hclust"
# method <- "sparcl_hclust"
# method <- "sparcl_kmeans"
method <- "impacc_hclust"
# method <- "sparse_gmm"
n_lambda <- 10
simul_study_id="6"
simul_id=1

# Template design
print(paste0("Simulation study ", simul_study_id))
if (method == "cosa_hclust") {
  mycolours <- lighten(
    c(
      "tan",
      lighten("maroon4", amount = seq(0.1, 0.7, length.out = 10)),
      "tan",
      "darkred",
      lighten("navy", amount = seq(0.1, 0.7, length.out = 10)),
      "darkred"
      # lighten("darkgreen", amount = seq(0.1, 0.7, length.out = 9))
    ),
    amount = 0.3
  )
}
if (method %in% c("sparcl_hclust", "sparcl_kmeans")) {
  mycolours <- lighten(
    c(
      "tan",
      lighten("maroon4", amount = seq(0.1, 0.7, length.out = 10)),
      "darkorange",
      "tan",
      "darkred",
      lighten("navy", amount = seq(0.1, 0.7, length.out = 10)),
      "darkorange",
      "darkred"
    ),
    amount = 0.3
  )
}
if (method == "impacc_hclust") {
  mycolours <- lighten(
    c(
      "tan",
      "tan",
      "darkred",
      lighten("darkgreen", amount = seq(0.1, 0.7, length.out = 9))
    ),
    amount = 0.3
  )
}
if (method == "sparse_gmm") {
  mycolours <- lighten(
    c(
      "tan",
      lighten("darkorange", amount = 0.3),
      "tan",
      lighten("darkorange", amount = 0.3),
      "tan",
      lighten("darkorange", amount = 0.3)
    ),
    amount = 0.3
  )
}
dimensionality <- c("", "", "")

# Saving table
if (type=="selection"){
  continuous_metrics <- c("precision", "recall", "F1_score")
  integer_metrics=NULL
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Selection_performances_", simul_id, "_merged.rds"))
} else {
  continuous_metrics <- c("rand", "ari", "jaccard")
  if (method=="sparse_gmm"){
    integer_metrics <- c("G", "time")
  } else {
    integer_metrics <- c("G", "q", "time")
  }
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
}
clustering_performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
if (method == "cosa_hclust") {
  lambda_list <- formatC(as.numeric(clustering_performances[, "lambda", 1]), format = "f", digits = 2)
  full_names <- c(
    "'G*'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[2:11], "'"),
    "'G*'",
    "'Consensus score'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[14:23], "'"),
    "'Consensus score'"
  )
  names(full_names) <- c(
    "single_run_star",
    paste0("single_run_cosa_star_", 1:n_lambda),
    "consensus_star", "consensus",
    paste0("consensus_cosa_star_", 1:n_lambda),
    "consensus_cosa"
  )
}
if (method %in% c("sparcl_hclust", "sparcl_kmeans")) {
  lambda_list <- formatC(as.numeric(clustering_performances[, "lambda", 1]), format = "f", digits = 2)
  full_names <- c(
    "'G*'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[2:11], "'"),
    "'GAP statistic'",
    "'G*'",
    "'Consensus score'",
    paste0("lambda[", 1:10, "]*'=", lambda_list[15:24], "'"),
    "'GAP statistic'", 
    "'Consensus score'"
  )
  names(full_names) <- c(
    "single_run_star",
    paste0("single_run_sparcl_star_", 1:n_lambda),
    "single_run_sparcl_star_gap",
    "consensus_star", "consensus",
    paste0("consensus_sparcl_star_", 1:n_lambda),
    "consensus_sparcl_star_gap", # TODO: reorder
    "consensus_sparcl"
  )
}
if (method == "impacc_hclust") {
  lambda_list <- formatC(as.numeric(clustering_performances[, "lambda", 1]), format = "f", digits = 2)
  full_names <- c(
    "'G*'",
    "'G*'",
    "'Consensus score'",
    paste0(seq(0.1, 0.9, by = 0.1))
  )
  names(full_names) <- c(
    "single_run_star",
    "consensus_star", "consensus",
    paste0("impacc_", seq(1, 9))
  )
}
if (method == "sparse_gmm") {
  full_names <- c(
    "'G*'",
    "BIC",
    "'G*'",
    "BIC",
    "'G*'",
    "MICL"
  )
  names(full_names) <- c(
    "single_run_star",
    "single_run_bic",
    "clustvarsel_star",
    "clustvarsel_bic",
    "varsellcm_star",
    "varsellcm_micl"
  )
}
if (method=="cosa_hclust"){
  algo_names <- c("Unweighted", "COSA", "Unweighted", "COSA")
} 
if (method%in%c("sparcl_hclust", "sparcl_kmeans")){
  algo_names <- c("Unweighted", "sparcl", "Unweighted", "sparcl")
}
if (method=="impacc_hclust"){
  algo_names <- c("Unweighted", "Unweighted", "IMPACC")
} 
if (method=="sparse_gmm"){
  algo_names <- c("Unweighted", "clustvarsel", "VarSelLCM")
}
n_simul <- length(list.files(path = paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id))) / 2
mytable <- NULL
if (type=="selection"){
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Selection_performances_", simul_id, "_merged.rds"))
} else {
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
}

tmpmedian <- apply(performances, c(1, 2), FUN = function(x) {
  median(as.numeric(x))
})
mymedian <- tmpmedian
mymedian[, c(continuous_metrics)] <- formatC(tmpmedian[, c(continuous_metrics)], format = "f", digits = 3)
if (type=="clustering"){
  mymedian[, c(integer_metrics)] <- formatC(tmpmedian[, c(integer_metrics)], format = "f", digits = 0)
  mymedian <- mymedian[, c(integer_metrics[-length(integer_metrics)], continuous_metrics, "time")]
}

tmpiqr <- apply(performances[, c(continuous_metrics, integer_metrics), ], c(1, 2), FUN = function(x) {
  IQR(as.numeric(x), na.rm = TRUE)
})
myiqr <- tmpiqr
myiqr[, c(continuous_metrics)] <- formatC(tmpiqr[, c(continuous_metrics)], format = "f", digits = 3)
if (type=="clustering"){
  myiqr[, c(integer_metrics)] <- formatC(tmpiqr[, c(integer_metrics)], format = "f", digits = 0)
  myiqr <- myiqr[, c(integer_metrics[-length(integer_metrics)], continuous_metrics, "time")]
}

tmptable <- matrix(paste0(mymedian[, ], " [", myiqr, "]"), ncol = ncol(mymedian))
colnames(tmptable) <- colnames(mymedian)

assign(paste0("tmptable", simul_id), tmptable)
mytable <- tmptable
mytable <- cbind(full_names[rownames(performances)], mytable)
mytable <- cbind(rep("", nrow(mytable)), mytable)

dir.create(paste0("Tables/Simulations_consensus_", method), showWarnings = FALSE)
write.xlsx(as.data.frame(mytable), paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".xlsx"),
           rowNames = FALSE, colNames = TRUE, overwrite = TRUE
)
mytable[, 2] <- gsub("'", "", mytable[, 2])
ids <- grep("lambda", mytable[, 2])
mytable[ids, 2] <- paste0("$", gsub("]*", "", gsub("lambda[", "\\lambda_", mytable[ids, 2], fixed = TRUE), fixed = TRUE), "$")
mytable[, 4] <- gsub("NA [NA]", "", mytable[, 4], fixed = TRUE)
write.table(mytable, paste0("Tables/Simulations_consensus_", method, "/Table_performances_", simul_study_id, ".txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)

# Saving figure
if (type=="clustering"){
  metric_list=c("G", "ari")
}
if (type=="selection"){
  metric_list=c("F1_score")
}
dir.create(paste0("Figures/Simulations_consensus_", method), showWarnings = FALSE)
for (metric in metric_list) {
  pdf(paste0("Figures/Simulations_consensus_", method, "/Boxplot_", metric, "_", simul_study_id, "_", simul_id, ".pdf"),
      width = 14, height = 7
  )
  par(mar = c(11, 5, 5, 6))
  
  if (type=="selection"){
    performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Selection_performances_", simul_id, "_merged.rds"))
    if (method=="sparse_gmm"){
      mylist <- list(NA, NA)
    }
    if (method=="impacc_hclust"){
      mylist <- list(NA, NA, NA)
    }
    if (method%in%c("cosa_hclust", "sparcl_hclust", "sparcl_kmeans")) {
      mylist <- list(NA)
    }
  } else {
    performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
    mylist <- list()
  }
  clustering_performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
  
  for (k in 1:nrow(performances)) {
    mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
    if (type=="clustering"){
      assign(paste0("median", k), median(as.numeric(performances[k, metric, ])))
    }
    if (type=="selection"){
      assign(paste0("median", k+1), median(as.numeric(performances[k, metric, ])))
    }
  }
  tmpfullnames <- full_names
  if (method!="sparse_gmm"){
    tmpfullnames <- paste0(tmpfullnames, "*' (q=", apply(clustering_performances[, "q", ], 1, FUN = function(x) {
      formatC(median(as.numeric(x), na.rm = TRUE), format = "f", digits = 0)
    }), ")'")
  }
  tmpfullnames <- gsub(" \\(q=NA\\)", "", tmpfullnames)
  names(tmpfullnames) <- names(full_names)
  xseq <- (1:length(full_names))
  if (metric %in% c("ari", "F1_score")) {
    ylim <- c(0, 1)
  } else {
    ylim <- c(0, 30)
  }
  if (metric=="F1_score"){
    myylab <- expression(F[1]*"-score")
  } else {
    myylab <- toupper(metric)
  }
  boxplot(
    at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
    whisklty = 1, range = 0, las = 1, main = dimensionality[simul_id], cex.main = 1.5,
    ylab = myylab, cex.lab = 1.5, xaxt = "n", ylim = ylim, frame = "F", boxwex = 0.25
  )
  abline(h = axTicks(2), lty = 3, col = "grey")
  if (method=="cosa_hclust"){
    zseq <- c(0.5, 1.5, 11.5, 13.5, length(xseq) + 0.5)
  }
  if (method%in%c("sparcl_hclust", "sparcl_kmeans")){
    zseq <- c(0.5, 1.5, 12.5, 14.5, length(xseq) + 0.5)
  }
  if (method=="impacc_hclust"){
    zseq <- c(0.5, 1.5, 3.5, length(xseq) + 0.5)
  }
  if (method=="sparse_gmm"){
    zseq <- c(0.5, 2.5, 4.5, length(xseq) + 0.5)
  }
  abline(v = zseq, lty = 2, col = "black")
  boxplot(
    at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
    whisklty = 1, range = 0, las = 1, add = TRUE,
    ylab = myylab, cex.lab = 1.5, xaxt = "n", frame = "F", boxwex = 0.25
  )
  if (method=="cosa_hclust"){
    id_list=c(1, 12, 13, 24)
  }
  if (method%in%c("sparcl_hclust", "sparcl_kmeans")){
    id_list=c(1, 13, 14, 26)
  }
  if (method=="impacc_hclust"){
    id_list=c(1, 2, 3)
  }
  if (method=="sparse_gmm"){
    id_list=c(1,3,5)
  }
  if (type=="selection"){
    id_list=id_list[-1]
  }
  for (id in id_list) {
    abline(h = eval(parse(text = paste0("median", id))), col = darken(mycolours[id], amount = 0.4), lty = 2)
  }
  axis(side = 1, at = xseq, labels = NA)
  for (i in 1:length(tmpfullnames)) {
    axis(
      side = 1, at = xseq[i],
      # labels = eval(parse(text = paste0("expression(", (tmpfullnames[rownames(performances)])[i], ")"))),
      labels = eval(parse(text = paste0("expression(", (tmpfullnames)[i], ")"))),
      las = 2, tick = FALSE
    )
  }
  if (method=="cosa_hclust"){
    myline=8
    axis(side = 1, at = xseq[c(2, 12)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(2, 12)]), labels = "G*", line = myline-0.5, tick = FALSE)
    axis(side = 1, at = xseq[c(14, 23)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(14, 23)]), labels = "G*", line = myline-0.5, tick = FALSE)
    axis(side = 3, at = zseq, labels = NA)
  }
  if (method%in%c("sparcl_hclust", "sparcl_kmeans")){
    myline=9
    axis(side = 1, at = xseq[c(2, 12)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(2, 12)]), labels = "G*", line = myline-0.5, tick = FALSE)
    axis(side = 1, at = xseq[c(15, 25)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(15, 25)]), labels = "G*", line = myline-0.5, tick = FALSE)
    axis(side = 3, at = zseq, labels = NA)
  }
  if (method=="impacc_hclust"){
    myline=8
    axis(side = 1, at = xseq[c(4, 12)] + c(-0.5, 0.5), labels = NA, line = myline)
    axis(side = 1, at = mean(xseq[c(4, 12)]), labels = "G*", line = myline-0.5, tick = FALSE)
    axis(side = 3, at = zseq, labels = NA)
  }
  if (method=="sparse_gmm"){
    axis(
      side = 3, at = zseq,
      labels = NA,
      tick = TRUE
    )
  }
  for (i in 1:length(algo_names)) {
    axis(
      side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean)[i],
      labels = algo_names[i],
      cex.axis = 1, tick = FALSE
    )
  }
  if (method!="sparse_gmm"){
    if (method=="impacc_hclust"){
      axis(side = 3, at = zseq[c(1,2)]+c(0.1,-0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at=mean(zseq[c(1,2)]),
        labels = "Single run",
        cex.axis = 1, tick = FALSE, line = 2
      )
      zseq2 <- zseq[c(1, 2, length(zseq))]
      axis(side = 3, at = zseq2[2:3]+c(0.1,-0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = apply(rbind(zseq2[-1], zseq2[-length(zseq2)]), 2, mean)[2],
        labels = "Consensus clustering",
        cex.axis = 1, tick = FALSE, line = 2
      )
    } else {
      zseq2 <- zseq[c(1, 2, 3)]
      axis(side = 3, at = zseq2[c(1,3)]+c(0.1,-0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at=mean(zseq2[c(1,3)]),
        labels = "Single run",
        cex.axis = 1, tick = FALSE, line = 2
      )
      zseq2 <- zseq[c(2, 3, length(zseq))]
      axis(side = 3, at = zseq2[2:3]+c(0.1,-0.1), labels = NA, line = 2.5)
      axis(
        side = 3, at = apply(rbind(zseq2[-1], zseq2[-length(zseq2)]), 2, mean)[2],
        labels = "Consensus clustering",
        cex.axis = 1, tick = FALSE, line = 2
      )
    }
  }
  dev.off() 
}
