##############################Michael Doerfler##################################
##############################2026-10-20########################################
#This script is for testing to make sure that the new full MLE inference functions are working
#as intended.

design <- expand.grid(N1 = c(500, 2000, 8000), psi = c(0.1, 0.3, 0.5, 0.7), KEEP.OUT.ATTRS = FALSE)
design$scenario_id <- seq_len(nrow(design))
design$mechanism   <- "multinomial"
design$K           <- 3
design$N2          <- 0.8 * design$N1
design$N_tot1      <- 10 * design$N1
design$N_tot2      <- 20 * design$N2
design$theta1      <- I(rep(list(c(.3, .6, .25, .5, .4, .3, .2)), nrow(design)))
design$theta2      <- I(rep(list(c(.4, .5, .20, .6, .5, .4, .2)), nrow(design)))
design$carry       <- I(rep(list(c(0, 0, 1)), nrow(design)))
saveRDS(design, here::here("sims", "test_sim_full_MLE", "design.rds"))

