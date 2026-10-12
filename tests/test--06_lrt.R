##########################Michael Doerfler#################################
##########################2026-10-12#######################################
#This script is for testing whether the 06_lrt script is performing correctly.


library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

# 1. Every single-carry model is saturated -> same loglik
fit_carryover(toy, carry = c(0, 0, 1), method = "optim")$loglik
fit_carryover(toy, carry = c(0, 1, 0), method = "optim")$loglik   # should agree
fit_carryover(toy, carry = c(1, 0, 0), method = "optim")$loglik   # should agree

# 2. Now, we simulate many datasets so that the null is true in each one. I simulate so that
# (1,0,1) is true.

results <- vector(mode = "list", length = 10000)

set.seed(2026)
K   <- 3
psi <- 0.1
theta1 <- c(.3, .6, .25, .5, .4, .3, .2)
theta2 <- c(.3, .5, .20, .6, .5, .4, .2)


for (b in 1:10000){
  data <- sim_multinomial(theta1, theta2, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  results[[b]] <- lrt_carryover(data = data, carry_null = c(1,0,1), carry_alt = c(0,0,1))
}

#Calculate the rejection rate.  If everything is working correctly, this should be around 0.05.
rej <- 0
for (b in 1:10000){
  rej <- rej + as.numeric(results[[b]]$p_value < 0.05)
}

rej_rate <- rej/10000

rej_rate

##3. Now, I do the previous step again, just with different simulation settings.  I simulate
#now so the null is (1,1,0) and the alternative is (1,0,0).

results <- vector(mode = "list", length = 10000)

set.seed(2027)
K   <- 3
psi <- 0.1
theta1 <- c(.3, .6, .20, .5, .4, .3, .5)
theta2 <- c(.3, .5, .20, .6, .5, .4, .2)


for (b in 1:10000){
  data <- sim_multinomial(theta1, theta2, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  results[[b]] <- lrt_carryover(data = data, carry_null = c(1,1, 0), carry_alt = c(1,0,0))
}

#Calculate the rejection rate.  If everything is working correctly, this should be around 0.05.
rej <- 0
for (b in 1:10000){
  rej <- rej + as.numeric(results[[b]]$p_value < 0.05)
}

rej_rate <- rej/10000

rej_rate

##4. Another set of simulations.  This time I simulate so the null is (1,1,1) and the alternative
## is (0,1,0).

results <- vector(mode = "list", length = 10000)

set.seed(2028)
K   <- 3
psi <- 0.1
theta1 <- c(.3, .6, .20, .5, .4, .3, .2)
theta2 <- c(.3, .5, .20, .6, .5, .4, .2)


for (b in 1:10000){
  data <- sim_multinomial(theta1, theta2, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  results[[b]] <- lrt_carryover(data = data, carry_null = c(1,1, 1), carry_alt = c(0,1,0))
}

#Calculate the rejection rate.  If everything is working correctly, this should be around 0.05.
rej <- 0
for (b in 1:10000){
  rej <- rej + as.numeric(results[[b]]$p_value < 0.05)
}

rej_rate <- rej/10000

rej_rate

##5. Now I do some power simulations.  I simulate so that the null is (1,0,1), and this is false.
##The alternative is (0,0,1), which is true.  I am going to do some simulations where I let the discrepancy
##in theta1 between the two populations grow.  We should see the rejection rate go up.

results1 <- vector(mode = "list", length = 10000)
results2 <- vector(mode = "list", length = 10000)
results3 <- vector(mode = "list", length = 10000)

set.seed(2029)
K   <- 3
psi <- 0.1
theta1 <- c(.3, .6, .20, .5, .4, .3, .2)
theta21 <- c(.35, .5, .40, .6, .5, .4, .2)
theta22 <- c(.55, .5, .40, .6, .5, .4, .2)
theta23 <- c(.75, .5, .40, .6, .5, .4, .2)


for (b in 1:10000){
  data1 <- sim_multinomial(theta1, theta21, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  data2 <- sim_multinomial(theta1, theta22, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  data3 <- sim_multinomial(theta1, theta23, N1 = 2000, N2 = 3000, N_tot1 = 10000, N_tot2 = 12000, psi = psi)
  
  results1[[b]] <- lrt_carryover(data = data1, carry_null = c(1,0,1), carry_alt = c(0,0,1))
  
  results2[[b]] <- lrt_carryover(data = data2, carry_null = c(1,0,1), carry_alt = c(0,0,1))
  
  results3[[b]] <- lrt_carryover(data = data3, carry_null = c(1,0,1), carry_alt = c(0,0,1))
}

#Calculate the rejection rate.  These should be increasing..
rej1 <- 0
rej2 <- 0
rej3 <- 0

for (b in 1:10000){
  rej1 <- rej1 + as.numeric(results1[[b]]$p_value < 0.05)
  rej2 <- rej2 + as.numeric(results2[[b]]$p_value < 0.05)
  rej3 <- rej3 + as.numeric(results3[[b]]$p_value < 0.05)
}

rej_rate1 <- rej1/10000; rej_rate2 <- rej2/10000; rej_rate3 <- rej3/10000

rej_rate1;rej_rate2;rej_rate3
