##############################Michael Doerfler##################################
##############################2026-10-20########################################
#This script is for testing to make sure that the new full MLE inference functions are working
#as intended.

design <- expand.grid(N_tot1 = c(1000, 2000, 5000), psi = c(0.1, 0.4, 0.7), KEEP.OUT.ATTRS = FALSE)
design$scenario_id <- seq_len(nrow(design))
design$mechanism   <- "multinomial"
design$K           <- 3
design$N1          <- 500
design$N2          <- 400
design$N_tot2      <- 2000
design$theta1      <- I(rep(list(c(.3, .6, .25, .5, .4, .3, .2)), nrow(design)))
design$theta2      <- I(rep(list(c(.4, .5, .20, .6, .5, .4, .2)), nrow(design)))
design$carry       <- I(rep(list(c(0, 0, 1)), nrow(design)))
saveRDS(design, here::here("sims", "test_sim_full_MLE_2", "design.rds"))

