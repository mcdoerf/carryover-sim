##################################Michael Doerfler##############################
##################################2026-10-20####################################
#This script is for testing the functions in the 04b_full_mle.R script.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

toy3 <- structure(
  list(
    K   = 3,
    psi = 0.1,
    counts = list(
      c("1111" = 6,  "1110" = 24, "1101" = 10, "1100" = 40,     # cohort 1: 2^(K+1) - 1 = 15 cells
        "1011" = 8,  "1010" = 32, "1001" = 12, "1000" = 48,
        "0111" = 5,  "0110" = 25, "0101" = 10, "0100" = 40,
        "0011" = 10, "0010" = 40, "0001" = 20),
      c("111" = 20, "110" = 40, "101" = 15, "100" = 75,          # cohort 2: 2^K - 1 = 7 cells
        "011" = 10, "010" = 50, "001" = 40)
    ),
    N_tot  = c(5000, 8000),
    n_dbar = 450,
    stream_names = c("pathology", "hospital", "death_cert")
  ),
  class = "cr_data"
)

carry = c(0,0,1)
# Testing chain_theta_given_N---------------------------------------------------
#According to my hand calculations, for N1 = 500 and N2 = 400 and carry = c(0,0,1), 
#chain_theta_given_N should return chain1 = (180/500, 80/(500-180), (50+40)/((500+400) - (180 + 80) - (150 + 60)))
# = (0.36, 0.25, 0.209).  Meanwhile, chain2 = (150/400, 60 / (400 - 150), 0.209) = 
#(0.375, 0.24, 0.209).

fc<- first_capture_counts(toy3)

chain_theta_given_N(N1 = 500, N2 = 400, b = fc[[1]], a = fc[[2]], carry = c(0,0,1), K = 3)

#I test again, but with carry switched to c(1,0,1).  In this case, we shoud have that
#the first entry in chain1 and chain2 is now (180 + 150) / (500 + 400) = 0.367.
chain_theta_given_N(N1 = 500, N2 = 400, b = fc[[1]], a = fc[[2]], carry = c(1,0,1), K = 3)


# Testing full_mle --------------------------------------------------------------
#Next, we are going to text the full_mle function.  We are going to make sure that 
#this indeed maximizes the likelihood.

res <- full_mle(data = toy3, carry = c(0,0,1))
res_cond <- fit_carryover(data = toy3, carry = c(0,0,1), method = "optim")


full_loglik(res$N_cont, res$chain1, res$chain2, toy3) #This should be greater
full_loglik(res_cond$N_hat, res_cond$chain1, res_cond$chain2, toy3) #This should be smaller

#On the other hand,

p_full <- c(res$chain1,      res$chain2[carry == 0])
p_cond <- c(res_cond$chain1, res_cond$chain2[carry == 0])

#This should be smaller since res maximizes the full likelihood
chain_loglik(qlogis(p_full), fc, psi = 0.1, K = 3, carry = carry)
#This should be greater since res_cond maximizes the conditional likelihood/
chain_loglik(qlogis(p_cond), fc, psi = 0.1, K = 3, carry = carry)

# Testing that profile is correct ----------------------------------------------

#At a fixed N, we are going to maximize the lieklihood over theta with optim and check
#that the result match chain_theta_given_N().

#I am writing this function that evaluates the full loglikelihood where par is
# betaU.

full_loglik_chain <- function(par, N, data, carry){
  K   <- data$K
  psi <- data$psi
  
  chain1 <- par[1:K]
  chain2 <- rep(NA, K)
  chain2[carry == 0] <- par[-(1:K)]
  chain2[carry == 1] <- chain1[carry == 1]
  
  stopifnot(
    "N must be c(N1, N2)" = is.numeric(N) && length(N) == 2,
    "chain1 and chain2 must have length K" = length(chain1) == K && length(chain2) == K,
    "chain parameters must lie in [0, 1]" = all(c(chain1, chain2) >= 0 & c(chain1, chain2) <= 1)
  )
  N1 <- N[[1]]; N2 <- N[[2]]
  
  fc  <- first_capture_counts(data)
  b   <- fc$b                      # length K + 1 (last entry: first caught by the anchor)
  a   <- fc$a                      # length K
  bK  <- b[1:K]
  nc1 <- sum(b); nc2 <- sum(a)
  N_tot1 <- data$N_tot[1]; n_dbar <- data$n_dbar
  
  # Outside the feasible region the likelihood is 0
  if (N1 < nc1 || N2 < nc2 || N1 > N_tot1 - n_dbar || N2 > data$N_tot[2]) return(-Inf)
  
  xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))      # 0 * log(0) = 0
  
  # At-risk counts: not caught by streams 1, ..., j - 1
  R1 <- N1 - head(c(0, cumsum(bK)), K)
  R2 <- N2 - head(c(0, cumsum(a)),  K)
  
  # (A) multinomial coefficients N! / (N - n_c)!
  A <- lgamma(N1 + 1) - lgamma(N1 - nc1 + 1) + lgamma(N2 + 1) - lgamma(N2 - nc2 + 1)
  
  # (B) first-capture (chain) terms for both cohorts
  B <- sum(xlogy(bK, chain1) + xlogy(R1 - bK, 1 - chain1)) +
    sum(xlogy(a,  chain2) + xlogy(R2 - a,  1 - chain2))
  
  # (B) anchor stage in cohort 1: b_{K+1} first caught by the anchor, N1 - n_c1 missed by it
  B_anchor <- xlogy(b[K + 1], psi) + (N1 - nc1) * log(1 - psi)
  
  # (C) anchor non-cases: n_dbar ~ Bin(N_tot1 - N1, psi), lgamma form so N1 can be non-integer
  C <- lgamma(N_tot1 - N1 + 1) - lgamma(N_tot1 - N1 - n_dbar + 1) - lgamma(n_dbar + 1) +
    n_dbar * log(psi) + (N_tot1 - N1 - n_dbar) * log(1 - psi)
  
  A + B + B_anchor + C
}

K = 3

carry = c(0,0,1)

start = runif(K + sum(carry == 0), 0.2, 0.7)

control <- list()

opt <- optim(
  par     = qlogis(start),
  fn      = function(beta) full_loglik_chain(plogis(beta), N = c(500, 400), data = toy3, carry = carry),
  method  = "BFGS",
  control = list(fnscale = -1, reltol = 1e-14, maxit = 1000)
)

#If everything is working correctly, the following should be the same:
plogis(opt$par)
chain_theta_given_N(500, 400, b = fc[[1]], a = fc[[2]], K = 3, carry = c(0,0,1))
