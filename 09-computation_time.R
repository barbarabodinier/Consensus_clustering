rm(list = ls())

library(colorspace)
library(openxlsx)

# Simulation parameters
method <- "hclust"
simul_study_id <- "2"

method_ids <- c(13:15)

all_times <- NULL
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/HPC/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
  all_times <- cbind(all_times, t(performances[method_ids, "time", ]))
  all_times <- cbind(all_times, NA)
}

pdf(paste0("Figures/Simulations_consensus_", method, "/Boxplot_time_", simul_study_id, ".pdf"),
  width = 6, height = 6
)
mycolours <- c("navy", lighten("navy", amount = 0.3), "darkred", "white")
par(mar = c(5, 5, 1, 1))
boxplot(
  all_times,
  col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
  whisklty = 1, range = 0, las = 1, cex.main = 1.5, xlim = c(0, ncol(all_times)),
  xlab = "Number of items", ylab = "Time (s)", cex.lab = 1.5, xaxt = "n", frame = "F", boxwex = 0.35
)
abline(v = seq(0, ncol(all_times), by = 4), lty = 2, col = "grey")
axis(
  side = 1, at = seq(2, ncol(all_times), by = 4),
  labels = c("150", "300", "600")
)
dev.off()
