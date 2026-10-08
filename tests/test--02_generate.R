###################Michael Doerfler####################
##################2026-10-08###########################
#Script for testing 02_generate.R.

library(here)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

#1) Testing sim.multinomial()-------------------------------------
#The idea here is that for large N^{(t)}, the simulated proportions should compare
#similarly to the true underlying probabilities if my code is working correctly.
set.seed(1)

d <- sim_multinomial(c(.3, .6, .2), c(.4, .5, .2), 1e6, 1e6, 2e6, 2e6, psi = .1)

p1 <- cell_probs(c(.3, .6, .2), 2, psi = .1)

p2 <- cell_probs(c(.4, .5, .2), 2)

round(rbind(sim = d$counts[[1]] / 1e6, true = p1[-length(p1)]), 3)   # should agree

round(rbind(sim = d$counts[[2]] / 1e6, true = p2[-length(p2)]), 3)   # should agree

#2) Testing sim_stratified()---------------------------------------
#Two strata with clearly different parameters. For large N, each stratum's
#simulated proportions should match that stratum's own cell probabilities.
set.seed(2)

theta1_mat <- rbind(c(.3, .6, .2),     # stratum 1, cohort 1
                    c(.7, .2, .5))     # stratum 2, cohort 1
theta2_mat <- rbind(c(.4, .5, .2),     # stratum 1, cohort 2
                    c(.6, .3, .4))     # stratum 2, cohort 2

ds <- sim_stratified(theta1_mat, theta2_mat,
                     N1_vec = c(1e6, 1e6), N2_vec = c(1e6, 1e6),
                     N_tot1_vec = c(2e6, 2e6), N_tot2_vec = c(2e6, 2e6),
                     psi = .1)

length(ds)          # should be 2

for (s in 1:2) {
  p1 <- cell_probs(theta1_mat[s, ], 2, psi = .1)
  p2 <- cell_probs(theta2_mat[s, ], 2)
  cat("\nStratum", s, "\n")
  print(round(rbind(sim = ds[[s]]$counts[[1]] / 1e6, true = p1[-length(p1)]), 3))   # should agree
  print(round(rbind(sim = ds[[s]]$counts[[2]] / 1e6, true = p2[-length(p2)]), 3))   # should agree
}

#Bad input should give an error:
try(sim_stratified(theta1_mat, theta2_mat, c(1e6, 1e6), c(1e6, 1e6),
                   c(2e6, 2e6), c(2e6, 2e6), psi = 1.5))   # psi out of range
