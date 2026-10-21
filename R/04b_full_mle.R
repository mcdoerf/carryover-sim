################################Michael Doerfler################################
################################2026-10-20######################################
#This R script contains code for doing inference with the full likelihood (not just conditional).

# chain_theta_given_N ----------------------------------------------------------
#This function returns the MLEs of chain 1 and chain 2 for KNOWN N^{(1)}, N^{(2)}.
#This function is important for creating the porfile likelihood.

#b is the first capture counts for cohort 1, a is the first capture counts for cohort 2,
# and carry is a vector of length K of 0's and 1's indicating which chain parameters to
# carry-over.

chain_theta_given_N <- function(N1, N2, b, a, K, carry){
  stopifnot(
    "b must have length K + 1 (streams 1..K plus the anchor)" = length(b) == K + 1,
    "a must have length K" = length(a) == K,
    "carry must have length K" = length(carry) == K
  )
  
  bK <- b[1:K]                                   # drop the anchor's first-capture count
  
  # At-risk counts: individuals not caught by streams 1, ..., j - 1
  R1 <- N1 - head(c(0, cumsum(bK)), K)
  R2 <- N2 - head(c(0, cumsum(a)),  K)
  
  # Carried parameters pool both populations; the rest use their own population
  num1 <- bK + a * carry
  num2 <- a + bK * carry
  chain1 <- num1 / (R1 + R2 * carry)
  chain2 <- num2 / (R2 + R1 * carry)
  
  # 0/0 at the lower bound of N: the MLE limit is 0
  chain1[num1 == 0] <- 0
  chain2[num2 == 0] <- 0
  
  list(
    chain1 = chain1,
    chain2 = chain2
  )
}


# 0 * log(0) = 0 convention
xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))

# full_profile_loglik ----------------------------------------------------------
# Profile log-likelihood in N = c(N1, N2), with the chain parameters profiled out.
#Note that par = c(N1, N2).
# Constants (factorials of cell counts, non-chain theta terms, psi powers for observed
# individuals) are dropped.

full_profile_loglik <- function(par, b, a, psi, K, carry, N_tot1, n_dbar){
  N1 <- par[1]; N2 <- par[2]
  
  chain  <- chain_theta_given_N(N1, N2, b, a, K, carry)
  chain1 <- chain$chain1
  chain2 <- chain$chain2
  
  bK  <- b[1:K]                       # streams 1..K only (drop the anchor)
  nc1 <- sum(b)                       # all observed cases in cohort 1 (incl. anchor-only)
  nc2 <- sum(a)
  
  R1 <- N1 - head(c(0, cumsum(bK)), K)
  R2 <- N2 - head(c(0, cumsum(a)),  K)
  
  A1 <- lgamma(N1 + 1) - lgamma(N1 - nc1 + 1)                       # (A), cohort 1
  A2 <- lgamma(N2 + 1) - lgamma(N2 - nc2 + 1)                       # (A), cohort 2
  A3 <- sum(xlogy(bK, chain1) + xlogy(R1 - bK, 1 - chain1))         # (B) chain, cohort 1
  A4 <- sum(xlogy(a,  chain2) + xlogy(R2 - a,  1 - chain2))         # (B) chain, cohort 2
  A5 <- (N1 - nc1) * log(1 - psi)                                   # (B) never caught: missed by anchor
  A6 <- lgamma(N_tot1 - N1 + 1) - lgamma(N_tot1 - N1 - n_dbar + 1) +
    (N_tot1 - N1 - n_dbar) * log(1 - psi)                       # (C) anchor non-cases
  
  A1 + A2 + A3 + A4 + A5 + A6
}


# full_mle ---------------------------------------------------------------------
#Takes data, a cr_data object, as input, along with carry (a K length vector of 0's and 1's
#indicating which chain parameters to carry over) along with option start values.

#This function returns N_hat, N_cont (the real maximizer), chain1, chain2 
#loglik, converged, and at_bound.

full_mle <- function(data, carry, start = NULL, control = list()){
  
  # Inputs from data
  fc     <- first_capture_counts(data)
  b      <- fc$b
  a      <- fc$a
  K      <- data$K
  psi    <- data$psi
  N_tot1 <- data$N_tot[1]
  N_tot2 <- data$N_tot[2]
  n_dbar <- data$n_dbar
  
  # Feasible box: N1 in [n_c1, N_tot1 - n_dbar], N2 in [n_c2, N_tot2]
  lower <- c(sum(b), sum(a))
  upper <- c(N_tot1 - n_dbar, N_tot2)
  
  # Starting values: the conditional estimates (explicitly conditional, to avoid recursion)
  if (is.null(start)) {
    start <- fit_carryover(data, carry, method = "optim")$N_tilde
  }
  start <- pmin(pmax(start, lower), upper)
  
  # Maximize the profile log-likelihood over (N1, N2)
  opt <- optim(
    par     = start,
    fn      = full_profile_loglik,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(fnscale = -1, factr = 1e3, parscale = c(100, 100), maxit = 1000),
                         control),
    b = b, a = a, psi = psi, K = K, carry = carry, N_tot1 = N_tot1, n_dbar = n_dbar
  )
  
  N_cont <- setNames(opt$par, c("N1", "N2"))
  
  # (Optional) integer step goes here; for now N_hat = N_cont
  N_hat <- N_cont
  
  # Chain parameters at the maximizer
  chain <- chain_theta_given_N(N_hat[1], N_hat[2], b, a, K, carry)
  
  list(
    N_hat     = N_hat,
    N_cont    = N_cont,
    chain1    = chain$chain1,
    chain2    = chain$chain2,
    loglik    = opt$value,
    converged = opt$convergence == 0,
    at_bound  = abs(N_hat - lower) < 1e-6 | abs(N_hat - upper) < 1e-6
  )
}


# vcov_from_sigma --------------------------------------------------------------
# Shared post-processing for the conditional and full fits: computes Sigma_NN with
# sigma_NN() and rescales it to the covariance of (N1_hat, N2_hat).
# If a chain parameter is on the boundary, or sigma_NN() fails, both are returned as NA.

vcov_from_sigma <- function(chain1, chain2, N_hat, K, psi, carry, e = 0, boundary = FALSE){
  S_NN <- if (isTRUE(boundary)) {
    matrix(NA_real_, 2, 2)
  } else {
    tryCatch(
      sigma_NN(chain1 = chain1, chain2 = chain2, N1 = N_hat[1], N2 = N_hat[2],
               K = K, psi = psi, carry = carry, e = e),
      error = function(err) matrix(NA_real_, 2, 2)
    )
  }
  dimnames(S_NN) <- list(c("N1", "N2"), c("N1", "N2"))
  
  list(
    Sigma_NN = S_NN,
    vcov_N   = outer(sqrt(N_hat), sqrt(N_hat)) * S_NN      # elementwise, not %*%
  )
}


# fit_carryover_full -----------------------------------------------------------
# Full-likelihood MLE of (N1, N2): uses the anchor's non-diseased count n_dbar as well
# as the case counts. Takes data (a cr_data object), carry (0/1 vector of length K),
# optional start = c(N1, N2), and optim controls passed to full_mle().
#
# Returns an object of class c("carryover_fit_full", "carryover_fit"), so confint()
# works unchanged. Extra fields: N_cont (real-valued maximizer), e (the added
# information for N1), at_bound.

fit_carryover_full <- function(data, carry, start = NULL, control = list(), tol_boundary = 1e-4){
  stopifnot("data must be a cr_data object" = inherits(data, "cr_data"))
  K   <- data$K
  psi <- data$psi
  
  # Same identifiability check as the conditional fit
  fc <- first_capture_counts(data)
  J  <- which(carry == 1)
  if (length(J) == 1 && fc$b[J] == 0 && fc$a[J] == 0) {
    stop("N2 not identified: no one in either cohort was first caught by the carried-over stream ", J)
  }
  
  # Maximize the full profile likelihood over (N1, N2)
  fit    <- full_mle(data, carry, start = start, control = control)
  N_hat  <- fit$N_hat
  chain1 <- unname(fit$chain1)      # unname: names from N_hat can trip all.equal() in sigma_NN()
  chain2 <- unname(fit$chain2)
  
  # Chain parameters on the boundary -> no Wald-type variance (same rule as chain_mle)
  betaU    <- c(chain1, chain2[carry == 0])
  boundary <- any(betaU < tol_boundary | betaU > 1 - tol_boundary)
  
  # Extra information for N1 from the anchor's non-diseased (N1-standardized scale)
  pi_hat <- N_hat[[1]] / data$N_tot[1]
  e_hat   <- pi_hat * psi / ((1 - pi_hat) * (1 - psi))
  
  V <- vcov_from_sigma(chain1 = chain1, chain2 = chain2, N_hat = N_hat,
                       K = K, psi = psi, carry = carry, e = e_hat, boundary = boundary)
  
  structure(list(
    N_hat       = N_hat,
    N_cont      = fit$N_cont,
    truncated   = fit$at_bound,        # confint() uses this; for the full fit it means "on the box edge"
    at_bound    = fit$at_bound,
    chain1      = chain1,
    chain2      = chain2,
    Sigma_NN    = V$Sigma_NN,
    vcov_N      = V$vcov_N,
    e           = e_hat,
    loglik      = fit$loglik,
    converged   = fit$converged,
    boundary    = boundary,
    method_used = "full_profile",
    carry       = carry,
    data        = data
  ), class = c("carryover_fit_full", "carryover_fit"))
}


# full_loglik ------------------------------------------------------------------
# Evaluates the full log-likelihood log L_T at given (N1, N2) and chain parameters,
# written in first-capture form. Used for testing only: the full MLE should give a
# value at least as large as the conditional MLE.
#
# Terms that are the same for every fit are dropped: the cell-count factorials
# prod(n_h!) and the non-chain theta terms (their MLEs are the same empirical
# proportions in both fits). Everything that depends on N or the chain parameters is kept.
#
# N      : c(N1, N2)
# chain1 : chain parameters theta_{j|0...0}^{(1)}, j = 1..K
# chain2 : chain parameters theta_{j|0...0}^{(2)}, j = 1..K (carried ones equal to chain1)
# data   : a cr_data object

full_loglik <- function(N, chain1, chain2, data){
  K   <- data$K
  psi <- data$psi
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