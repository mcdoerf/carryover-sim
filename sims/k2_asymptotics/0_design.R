############################Michael Doerfler####################################
############################2026-10-14##########################################
#This is the R script for putting together the grid in the k2_asymptotics
#simulation study.
N1_vals <- round(exp(seq(log(50), log(20000), length.out = 19)))
N1_vals <- sort(unique(c(N1_vals, 159)))    # 20 values including the CRISP case
N1_vals

design <- expand.grid(N1 = N1_vals, psi = c(200/1029), KEEP.OUT.ATTRS = FALSE)
design$scenario_id <- seq_len(nrow(design))
design$mechanism   <- "multinomial"
design$K           <- 2
design$N2          <- round(0.8 * design$N1)
design$N_tot1      <- round(design$N1 / (159/1029))
design$N_tot2      <- round(design$N2 / (0.18))
design$theta1      <- I(rep(list(c(.4, 0.6, 0.3)), nrow(design)))
design$theta2      <- I(rep(list(c(.5, 0.55, 0.3)), nrow(design)))
design$carry       <- I(rep(list(c(0, 1)), nrow(design)))
saveRDS(design, here::here("sims", "k2_asymptotics", "design.rds"))


#A few checks:

# the carry-over really holds in every scenario
stopifnot(all(vapply(seq_len(nrow(design)), function(i)
  design$theta1[[i]][3] == design$theta2[[i]][3], logical(1))))

# how likely are zero counts? smallest expected first-capture count per scenario
design$min_expected <- vapply(seq_len(nrow(design)), function(i) {
  t1 <- design$theta1[[i]][c(1, 3)]; t2 <- design$theta2[[i]][c(1, 3)]
  p1 <- c(t1[1], (1 - t1[1]) * t1[2], (1 - t1[1]) * (1 - t1[2]) * design$psi[i])   # b1, b2, anchor-only
  p2 <- c(t2[1], (1 - t2[1]) * t2[2])                                               # a1, a2
  min(design$N1[i] * p1, design$N2[i] * p2)
}, numeric(1))
print(design[, c("scenario_id", "N1", "psi", "min_expected")])
