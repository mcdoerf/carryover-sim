###################Michael Doerfler#########################
##################2026-10-12################################
#This is an R script for testing the 05_bayes script.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

##Simulate toy data
#K = 3. theta vectors are in the paper's order:
#theta_1, theta_2|1, theta_2|0, theta_3|11, theta_3|10, theta_3|01, theta_3|00.
#The chain parameters are positions 1, 3, 7 (chain_pos(3)).
set.seed(2026)
K   <- 3
psi <- 0.1
theta1 <- c(.3, .6, .25, .5, .4, .3, .2)
theta2 <- c(.4, .5, .20, .6, .5, .4, .2)
toy <- sim_multinomial(theta1, theta2, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
toy
fc <- first_capture_counts(toy)
fc

##First, we will run the frequentist inference,  THe Bayes and frequentist inference should
#look fairy similar.

fit_frequentist <- fit_carryover(data = toy, carry = c(0,0,1), method = "optim")
fit_frequentist_summary <- tidy(fit_frequentist)

#Now we run the Bayesian inference.

fit_bayesian <- fit_bayes(data = toy, carry = c(0,0,1))
fit_bayesian_summary <- tidy(fit_bayesian)

#Compare the results.  If things are working correctly, these should look pretty similar.

fit_frequentist_summary
fit_bayesian_summary
