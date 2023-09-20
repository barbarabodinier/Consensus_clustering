library(abind)

args <- commandArgs(trailingOnly = TRUE)
simul_study_id <- as.numeric(args[1])
algo <- as.character(args[2])

print(paste("ID of the simulation study:", simul_study_id))
cat("\n")

setwd(paste0("../../Results/Simulations_consensus_", algo, "/Simulations_", simul_study_id))

myfiles <- list.files(pattern = "Selection_performances")
params_id_list <- unique(gsub("_.*", "", gsub("Selection_performances_", "", myfiles)))

for (params_id in params_id_list) {
  print(paste0("Simulation ID: ", params_id))

  myfiles <- list.files(pattern = "Selection_performances")
  myfiles <- myfiles[grep(paste0("Selection_performances_", params_id, "_.*", ".rds"), myfiles)]
  myfiles <- myfiles[!grepl("merged", myfiles)]
  print(length(myfiles))

  results <- NULL
  for (j in 1:min(1000, length(myfiles))) {
    tmp <- readRDS(myfiles[j])
    results <- abind(results, tmp, along = 3)
  }

  saveRDS(results, paste0("Selection_performances_", params_id, "_merged.rds"))
}
