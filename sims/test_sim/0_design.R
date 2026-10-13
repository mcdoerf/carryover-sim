##############################Michael Doerfler##################################
##############################2026-10-13########################################
#This script is part of a test simulation study to make sure the simulation functions in 
# 09_sim_engine.R are working correctly.
#This is the design script in which the design grid is specified for the test_sim
#simulation study.

design <- expand.grid(N1 = c(500, 2000, 8000), psi = c(0.05, 0.1), KEEP.OUT.ATTRS = FALSE)
design$scenario_id <- seq_len(nrow(design))
design$mechanism   <- "multinomial"
design$K           <- 3
design$N2          <- 0.8 * design$N1
design$N_tot1      <- 10 * design$N1
design$N_tot2      <- 20 * design$N2
design$theta1      <- I(rep(list(c(.3, .6, .25, .5, .4, .3, .2)), nrow(design)))
design$theta2      <- I(rep(list(c(.4, .5, .20, .6, .5, .4, .2)), nrow(design)))
design$carry       <- I(rep(list(c(0, 0, 1)), nrow(design)))
saveRDS(design, here::here("sims", "test_sim", "design.rds"))

