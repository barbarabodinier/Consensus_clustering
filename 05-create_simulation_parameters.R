rm(list = ls())

simul_study_id <- 4

nc <- 5
equal_size <- FALSE
n_tot <- 150
p <- 100
ev_xc <- 0.6
nu_xc <- 0.2
v_min <- 0
v_max <- 0

mylist <- expand.grid(nc = nc, equal_size = equal_size, n_tot = n_tot, p = p, ev_xc = ev_xc, nu_xc = nu_xc, v_min = v_min, v_max = v_max)
mylist <- cbind(simulation_id = 1:nrow(mylist), mylist)

dir.create("Simulation_parameters", showWarnings = FALSE)
write.table(mylist, paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"), sep = "\t", row.names = FALSE)
