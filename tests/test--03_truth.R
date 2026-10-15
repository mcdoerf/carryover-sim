##############Michael Doerfler###############
#############2026-10-08######################
#This R script is for testing the functions in the 03_truth.R file.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}


#1) Exact case ---------------------------------------------------------------
#When only theta_{K|0...0} is carried over, cohort 1's estimator is
#N1_tilde = n_c - b_{K+1} + b_{K+1}/psi, whose variance is known exactly:
#Sigma[1,1] = (1 - p_caught)(1 - psi)/psi, p_caught = P(caught by one of the K streams).
set.seed(1)
K <- 3; psi <- 0.1
chain1 <- runif(K, .2, .7)
chain2 <- runif(K, .2, .7)
carry  <- c(0, 0, 1)
chain1[carry == 1] <- chain2[carry == 1]

S <- sigma_NN(chain1, chain2, 2e4, 3e4, K, psi, carry)
p_caught <- 1 - prod(1 - chain1)
c(function_value = S[1, 1], exact = (1 - p_caught) * (1 - psi) / psi)   # should agree

#2) Compare with the full information matrix --------------------------------
#Build Sigma^{-1} from Theorem 2 using ALL parameters (identifiable ones too),
#with derivatives of cell_probs() taken numerically, invert it, and take the
#N block. This doesn't use the Schur-complement shortcut, so it's an independent check.
full_sigma_NN <- function(theta1, theta2, N1, N2, K, psi, carry) {
  lambda <- c(N1, N2) / (N1 + N2)
  shared_pos <- chain_pos(K)[carry == 1]
  beta <- c(theta1, theta2[-shared_pos])
  d <- length(beta)
  probs <- function(beta) {
    t1 <- beta[1:(2^K - 1)]
    t2 <- numeric(2^K - 1)
    t2[-shared_pos] <- beta[-(1:(2^K - 1))]
    t2[shared_pos]  <- t1[shared_pos]
    list(cell_probs(t1, K, psi = psi), cell_probs(t2, K))
  }
  p <- probs(beta)
  eps <- 1e-6
  grad <- lapply(1:2, function(t) sapply(1:d, function(r) {
    up <- beta; up[r] <- up[r] + eps
    dn <- beta; dn[r] <- dn[r] - eps
    (probs(up)[[t]] - probs(dn)[[t]]) / (2 * eps)
  }))
  Sinv <- matrix(0, d + 2, d + 2)
  for (t in 1:2) {
    G <- grad[[t]]; pt <- p[[t]]; z <- length(pt)      # z = the all-zero history
    Sinv[1:d, 1:d]   <- Sinv[1:d, 1:d] + lambda[t] * t(G / pt) %*% G
    Sinv[1:d, d + t] <- Sinv[d + t, 1:d] <- -sqrt(lambda[t]) * G[z, ] / pt[z]
    Sinv[d + t, d + t] <- (1 - pt[z]) / pt[z]
  }
  solve(Sinv)[d + 1:2, d + 1:2]
}

for (carry in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1))) {
  theta1 <- runif(2^K-1, .2, .7)
  th2 <- runif(2^K - 1, .2, .7)
  th2[chain_pos(K)][carry == 1] <- theta1[chain_pos(K)][carry == 1]
  ch1 <- theta1[chain_pos(K)]
  ch2 <- th2[chain_pos(K)]
  diff <- max(abs(sigma_NN(ch1, ch2, 2e4, 3e4, K, psi, carry) -
                    full_sigma_NN(theta1, th2, 2e4, 3e4, K, psi, carry)))
  cat("carry =", carry, "  max difference =", signif(diff, 2), "\n")   # should be ~1e-9
}

#3) Basic properties ---------------------------------------------------------
isSymmetric(S)              # should be TRUE
all(eigen(S)$values > 0)    # positive definite: should be TRUE


#4) Time saving --------------------------------------------------------------
#We look at how much time the Schur complement shortcut saves in calculating compared to 
#inverting the full matrix.

carry = c(1, 0, 1)


tic()
for(i in 1:1000){
  theta1 <- runif(2^K-1, .2, .7)
  th2 <- runif(2^K - 1, .2, .7)
  th2[chain_pos(K)][carry == 1] <- theta1[chain_pos(K)][carry == 1]
  ch1 <- theta1[chain_pos(K)]
  ch2 <- th2[chain_pos(K)]
  sigma_NN(ch1, ch2, 2e4, 3e4, K, psi, carry)
}
toc()

tic()
for(i in 1:1000){
  theta1 <- runif(2^K-1, .2, .7)
  th2 <- runif(2^K - 1, .2, .7)
  th2[chain_pos(K)][carry == 1] <- theta1[chain_pos(K)][carry == 1]
  full_sigma_NN(theta1, th2, 2e4, 3e4, K, psi, carry)
}
toc()

##Schur complement shortcut speeds things up by about 800!