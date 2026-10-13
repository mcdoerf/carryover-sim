##############################Michael Doerfler##################################
##############################2026-10-12########################################
#This is for testing the 08_stratified R script.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

# Testing pool_strata() ---------------------------------------------------------

#First, I simulate some toy data.

set.seed(2026)

theta1_mat <- matrix(c(0.4, 0.55, 0.45, 0.3, 0.65, 0.25), nrow = 2, ncol = 3, byrow = TRUE)
theta2_mat <- matrix(c(0.35, 0.65, 0.45, 0.3, 0.75, 0.25), nrow = 2, ncol = 3, byrow = TRUE)

N1_vec <- c(500, 2000)
N2_vec <- c(3000, 500)

N_tot1_vec <- c(10000, 10000)
N_tot2_vec <- c(10000, 10000)

psi <- 0.1

d_list <- sim_stratified(theta1_mat, theta2_mat, N1_vec, N2_vec, N_tot1_vec, N_tot2_vec, psi)

d_list[[1]]
d_list[[2]]

#According to my hand calculations, if pool_strata() is working correctly then it should return
# Cohort 1: nc = 1401, Ntot = 20000
# 111   110   101   100   011   010   001
#  59   438    24   252    50   464   114
# Cohort 2: nc = 2122 Ntot = 20000
#  11   10   01 
# 753  385  984
#ndbar: 1810

pool_strata(d_list)

# Testing fit_stratified()-------------------------------------------------------
#Next, I test fit_stratified().  Let us first do fits for each strata seperately using 
# carryover_fit().

fit1 <- fit_carryover(data = d_list[[1]], carry = c(0, 1), method = "optim")
fit2 <- fit_carryover(data = d_list[[2]], carry = c(1,1), method = "optim")

tidy(fit1)
tidy(fit2)

#So, if things are working correctly, once we pool the estimate of N1 should be 503 + 1924 = 2427
#and the estimate of N2 should be 2859.777 + 457.8071 = 3317.58.
#Morevoer, the se for N1 should be sqrt(39.12^2 + 93.43^2) = 101.3, and the se for N2 should
#be sqrt(250.38^2 + 32.72^2) = 252.5.
#Additionally, according to my hand calculations, the Wald interval for N1 should be (2228.5, 2625.5),
#and the Wald interval for N2 should be (2822.7, 3812.5).

fit_strat <- fit_stratified(data_list = d_list, carry_list = list(c(0,1), c(1,1)), method = "optim")
tidy(fit_strat)

# Testing bayes_stratified()----------------------------------------------------
fit_bayes1 <- fit_bayes(data = d_list[[1]], carry = c(0,1))
fit_bayes2 <- fit_bayes(data = d_list[[2]], carry = c(1,1))

tidy(fit_bayes1)
tidy(fit_bayes2)
#If everything is working correclty, the output from fit_strat_bayes should look roughly similar
#to that from fit_strat.

fit_strat_bayes <- bayes_stratified(data_list = d_list, carry_list = list(c(0,1), c(1,1)))
tidy(fit_strat_bayes)
