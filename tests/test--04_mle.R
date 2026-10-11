###############################Michael Doerfler#############################
##############################2026-10-09####################################
#This script is for testing the 04_mle R script.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

# Testing chain_loglik --------------------------------------------------------

logit <- function(x){
  log(x / (1-x))
}

K = 3

carry <- c(1, 0, 1)

betaU <- c(0.4, 0.45, 0.35, 0.25)

par <- logit(betaU)

#toy data

# Toy data: K = 3 streams; cohort 1's last digit is the anchor stream
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

#Extract first capture counts

fc <- first_capture_counts(toy3)

####By hand, I calculated that the joint log-likelihood should be -584.966 for this example.
###So, if the chain_loglik function is working correctly, that is what it should return.

#Plug into function.

chain_loglik(par = par, fc = fc, psi = 0.2, K = 3, carry = c(1,0,1))

#Next, testing to make sure it catches incorrect input.

chain_loglik(par = par, fc = fc, psi = 1.2, K = 3, carry = c(1,0,1))
chain_loglik(par = par, fc = fc, psi = 0.2, K = 5, carry = c(1,0,1))
chain_loglik(par = par, fc = fc, psi = 0.2, K = 3, carry = c(-2,0,1))


# Testing chain_mle -----------------------------------------------------

#Helper: N-tilde from fitted chains (n_c / p_c).
N_from_chains <- function(fit, fc, psi) {
  c(sum(fc$b) / (1 - (1 - psi) * prod(1 - fit$chain1)),
    sum(fc$a) / (1 - prod(1 - fit$chain2)))
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
toy <- sim_multinomial(theta1, theta2, N1 = 500, N2 = 400, N_tot1 = 5000, N_tot2 = 8000, psi = psi)
toy
fc <- first_capture_counts(toy)
fc

##1) Closed form (J = K) matches the formulas computed by hand 
fit_closed <- chain_mle(fc, psi, K, carry = c(0, 0, 1), method = "closed")
b <- fc$b; a <- fc$a
N1_hand    <- sum(b) - b[K + 1] + b[K + 1] / psi
theta_hand <- b[K] / (b[K] + b[K + 1] / psi)
N2_hand    <- sum(a) - a[K] + a[K] / theta_hand
rbind(chain_mle = N_from_chains(fit_closed, fc, psi), by_hand = c(N1_hand, N2_hand))   # rows should agree

##2) optim agrees with the closed form when J = K
fit_optim <- chain_mle(fc, psi, K, carry = c(0, 0, 1), method = "optim")
fit_optim$converged                                                     # should be TRUE
rbind(closed = N_from_chains(fit_closed, fc, psi),
      optim  = N_from_chains(fit_optim,  fc, psi))                       # rows should agree
fit_closed$loglik - fit_optim$loglik                                    # should be ~0

##3) optim agrees with the full conditional likelihood for other carry-overs 
#Independent reference: maximize the full conditional likelihood over ALL
#parameters (identifiable ones too), using cell_probs(). Slow, but only used here.
full_mle_N <- function(data, carry) {
  K  <- data$K; m <- 2^K - 1
  sp <- chain_pos(K)[carry == 1]
  n1 <- data$counts[[1]]; n2 <- data$counts[[2]]
  neg_loglik <- function(eta) {
    b  <- plogis(eta)
    t1 <- b[1:m]
    t2 <- numeric(m); t2[-sp] <- b[-(1:m)]; t2[sp] <- t1[sp]
    p1 <- cell_probs(t1, K, psi = data$psi); p1 <- p1[-length(p1)]
    p2 <- cell_probs(t2, K);                 p2 <- p2[-length(p2)]
    -(sum(n1 * log(p1 / sum(p1))) + sum(n2 * log(p2 / sum(p2))))
  }
  fit <- optim(rep(0, 2 * m - length(sp)), neg_loglik, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-14))
  b  <- plogis(fit$par)
  t1 <- b[1:m]
  t2 <- numeric(m); t2[-sp] <- b[-(1:m)]; t2[sp] <- t1[sp]
  c(sum(n1) / (1 - (1 - data$psi) * prod(1 - t1[chain_pos(K)])),
    sum(n2) / (1 - prod(1 - t2[chain_pos(K)])))
}
for (carry in list(c(0, 1, 0), c(1, 0, 0), c(1, 0, 1))) {
  fit <- chain_mle(fc, psi, K, carry = carry, method = "optim")
  cat("carry =", carry,
      "  chain_mle =", round(N_from_chains(fit, fc, psi), 3),
      "  full likelihood =", round(full_mle_N(toy, carry), 3),
      "  converged =", fit$converged, "\n")
}
#the two sets of numbers should agree

##4) Shared parameters are equal in both chains-
fit <- chain_mle(fc, psi, K, carry = c(1, 0, 1), method = "optim")
fit$chain1[c(1, 3)] - fit$chain2[c(1, 3)]    # should be exactly 0

##5) The answer does not depend on the starting value
set.seed(1)
starts <- replicate(5, runif(K + 1, 0.05, 0.95), simplify = FALSE)
sapply(starts, function(s)
  round(N_from_chains(chain_mle(fc, psi, K, c(1, 0, 1), "optim", start = s), fc, psi), 2))
#all columns should (essentially) agree

##6) Analytic gradient matches a numerical gradient 
num_grad <- function(f, x, ..., h = 1e-6) {
  sapply(seq_along(x), function(i) {
    up <- x; dn <- x
    up[i] <- up[i] + h; dn[i] <- dn[i] - h
    (f(up, ...) - f(dn, ...)) / (2 * h)
  })
}
set.seed(2)
for (carry in list(c(0, 0, 1), c(0, 1, 0), c(1, 0, 1), c(1, 1, 1))) {
  par <- rnorm(K + sum(carry == 0))
  analytic <- gradient_chain_loglik(par, fc, psi, K, carry)
  numeric  <- num_grad(chain_loglik, par, fc = fc, psi = psi, K = K, carry = carry)
  cat("carry =", carry, "  max relative difference =",
      signif(max(abs(analytic - numeric) / pmax(1, abs(numeric))), 2), "\n")
}
#differences should be ~1e-7 or smaller

##7) Closed form refuses carry-overs other than J = K 
try(chain_mle(fc, psi, K, carry = c(0, 1, 0), method = "closed"))   # should give an error

##8) With large N, the estimates are close to the truth 
set.seed(3)
big     <- sim_multinomial(theta1, theta2, 1e6, 1e6, 3e6, 3e6, psi = psi)
fc_big  <- first_capture_counts(big)
fit_big <- chain_mle(fc_big, psi, K, carry = c(0, 0, 1), method = "optim")
round(rbind(estimate = fit_big$chain1, truth = theta1[chain_pos(K)]), 3)   # should agree
round(rbind(estimate = fit_big$chain2, truth = theta2[chain_pos(K)]), 3)   # should agree
N_from_chains(fit_big, fc_big, psi) / 1e6                                   # both should be close to 1


# Testing fit_carryover ---------------------------------------------------------------
#By hand calculations, N1_hat should be 413, N2_hat should be 323.2, and V(N1_hat) should be
# 1,170 for the toy data.

fit_carryover(data = toy, carry = c(0,0,1), method = "closed")
fit_carryover(data = toy, carry = c(0,0,1), method = "optim")

# Testing confint.carryover_fit()-----------------------------------------------------

#According to my hand calculations for the toy data, the wald interval for N1 should be 
#(346, 480) and for N2 it should be (254.1, 392.3).
object <- fit_carryover(data = toy, carry = c(0,0,1), method = "optim")

confint.carryover_fit(object = object, type = "wald")

#According to my hand calculations for the toy data, the logit interval for N1 should be
#(350.7, 485.2) and for N2 it should be (260.7, 399.84).

confint.carryover_fit(object = object, type = "logit")

##Let us also check to make sure the truncation is being implemented correctly.

#Testing confint.carryover_fit() truncation
c1 <- c("111" = 12, "110" = 48, "101" = 20, "100" = 80, "011" = 15, "010" = 60, "001" = 25)
c2 <- c("11" = 30, "10" = 90, "01" = 45)
#With these counts and carry = c(0, 1): N1_tilde = 485, N2_tilde = 315, n_c = (260, 165)

#A) No truncation: plenty of room in both cohorts
dA <- cr_data(c1, c2, psi = 0.1, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)
fA <- fit_carryover(dA, c(0, 1))
fA$truncated                        # FALSE FALSE
confint(fA); confint(fA, type = "logit")

#B) N2 truncated: N_tot2 = 300 is below N2_tilde = 315
dB <- cr_data(c1, c2, psi = 0.1, N_tot1 = 5000, N_tot2 = 300, n_dbar = 470)
fB <- fit_carryover(dB, c(0, 1))
fB$N_hat; fB$truncated              # N2 = 300; FALSE TRUE
confint(fB); confint(fB, type = "logit")
#N2: logit falls back to Wald, so both rows for N2 should be identical, upper = 300

#C) N1's interval clipped above: U1 = 990 - 470 = 520, inside the raw Wald interval
dC <- cr_data(c1, c2, psi = 0.1, N_tot1 = 990, N_tot2 = 8000, n_dbar = 470)
fC <- fit_carryover(dC, c(0, 1))
fC$truncated                        # FALSE FALSE (estimate 485 is below 520)
confint(fC); confint(fC, type = "logit")
#N1 upper limit should be 520 for both types

#D) Lower limits clipped at n_c: few people caught only by the last stream
c1s <- c1; c1s["001"] <- 1
c2s <- c2; c2s["01"]  <- 2
dD <- cr_data(c1s, c2s, psi = 0.1, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)
fD <- fit_carryover(dD, c(0, 1))
fD$N_hat                            # 245, 122.27
confint(fD); confint(fD, type = "logit")
#N1 lower limit should be n_c = 236 (raw Wald lower is 226.4); N2 lower should be 122


# Testing tidy.carryover_fit()-----------------------------------------------------

tidy(x = fit_carryover(data = toy, carry = c(0,0,1), method = "optim"))
tidy(x = fit_carryover(data = dA, carry = c(0,1), method = "optim"))
tidy(x = fit_carryover(data = dB, carry = c(0,1), method = "optim"))
