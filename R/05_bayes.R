###############################Michael Doerfler#################################
###############################2026-10-11#######################################
#This script contains the functions for the Bayesian inference.  

# fit_bayes() ------------------------------------------------------------------
#This function takes a cr_data object as input, as well as a carry vector of 0's and 1's
#of length K indicating which parameters in the non-identifiable chain carry-over.
#init is a vector of length 2^{K+1} - (# of shared parameters).  The initial values should
#be in the order theta1, theta21, theta2bar1, ... shared parameters, pi1, pi2.

fit_bayes <- function(data, carry, n_keep = 2000, burn = 1000, thin = 5, prior = c(0.5, 0.5)){
  #Extract important variables from data.
  K <- data$K
  
  psi <- data$psi
  
  counts1 <- data$counts[[1]]
  
  counts2 <- data$counts[[2]]
  
  n0bullet1 <- data$N_tot[1] - sum(data$counts[[1]]) - data$n_dbar
  
  n0bullet2 <- data$N_tot[2] - sum(data$counts[[2]])
  
  nc1 <- sum(counts1)
  
  nc2 <- sum(counts2)
  
  fc <- first_capture_counts(data = data)
  
  zero_vector <- strrep("0", 1:K)
  
  zero_history_counts1 <- rep(NA, K)
  
  zero_history_counts2 <- rep(NA, K)
  
  for (j in 1:K){
    zero_history_counts1[j] <- sum(counts1[startsWith(names(counts1), zero_vector[j])])
    
    zero_history_counts2[j] <- sum(counts2[startsWith(names(counts2), zero_vector[j])])
  }
  
  fc1 <-fc$b[-length(fc$b)]
  
  fc2 <- fc$a
  
  shape1_1 <- fc1 + fc2 * carry + prior[1]
  
  fail_1   <- zero_history_counts1 + zero_history_counts2 * carry + prior[2]
  
  shape1_2 <- fc2 + prior[1]
  
  fail_2   <- zero_history_counts2 + prior[2]
  
  #Calculate number of total runs we will need
  P <- burn
  
  k <- thin
  
  M <- n_keep
  
  B <- P + (M * k)
  
  cond_names <- cond_hist_names(K)
  
  #Create matrices to save output
  
  chain1_draws <- matrix(NA, B+1, K)
  
  chain2_draws <- matrix(NA, B+1, K)
  
  pi_draws <- matrix(NA, B+1, 2)
  
  N_draws <- matrix(NA, B+1, 2)
  
  #Initializing values -------------------------------------------------------------------
  
  # Try the MLE as a starting point; fall back to simple values if it fails
  fit <- tryCatch(fit_carryover(data = data, carry = carry, method = "optim"),
                  error = function(e) NULL)
  
  use_mle <- !is.null(fit) && isTRUE(fit$converged) && !isTRUE(fit$boundary) &&
    all(is.finite(c(fit$N_hat, fit$chain1, fit$chain2)))
  
  if (use_mle) {
    chain1_init <- fit$chain1
    chain2_init <- fit$chain2
    N_init      <- fit$N_hat
  } else {
    chain1_init <- rep(0.5, K)
    chain2_init <- rep(0.5, K)          # shared entries are automatically equal
    p_c1 <- 1 - (1 - psi) * prod(1 - chain1_init)
    p_c2 <- 1 - prod(1 - chain2_init)
    N_init <- pmin(c(nc1 / p_c1, nc2 / p_c2),
                   c(data$N_tot[1] - data$n_dbar, data$N_tot[2]))
  }
  
  chain1_draws[1, ] <- chain1_init
  chain2_draws[1, ] <- chain2_init
  N_draws[1, ]      <- N_init
  # keep pi strictly inside (0, 1), even if N_init equals N_tot
  pi_draws[1, ]     <- pmin(pmax(N_init / data$N_tot, 0.001), 0.999)
    
    #Main MCMC loop --------------------------------------------------------------------
    
    for (b in 2:(B+1)){
      
      p_miss1 <- (1 - psi) * prod(1 - chain1_draws[b-1, ])
      
      p_miss2 <- prod(1 - chain2_draws[b-1, ])
      
      #Update n01 and n02
      
      n01 <- rbinom(1, n0bullet1, pi_draws[b-1, 1] * p_miss1 / (pi_draws[b-1, 1] * p_miss1 + (1 - pi_draws[b-1,1]) * (1 - psi)))
      
      n02 <- rbinom(1, n0bullet2, pi_draws[b-1, 2] * p_miss2 / (pi_draws[b-1,2] * p_miss2 + (1 - pi_draws[b-1,2])))
      
      #Update N1 and N2
      
      N_draws[b, ] <- c(nc1 + n01, nc2 + n02)
      
      #Update pi1 and pi2
      
      pi1 <- rbeta(1, N_draws[b, 1] + prior[1], data$N_tot[1] - N_draws[b, 1] + prior[2])
      
      pi2 <- rbeta(1, N_draws[b, 2] + prior[1], data$N_tot[2] - N_draws[b, 2] + prior[2])
      
      pi_draws[b, ] <- c(pi1, pi2)
      
      
     #Update chain parameters in population 1
      
      chain1_draws[b,] <- rbeta(K, shape1_1, fail_1 + n01 + n02 * carry)
      
      #Update chain parameters in population 2
      
      chain2_draws[b,] <- rbeta(K, shape1_2, fail_2 + n02)
      
      chain2_draws[b, ][carry == 1] <- chain1_draws[b, ][carry == 1]
      
      
    }
    
    #Apply thinning and burnin
    
    colnames(chain1_draws) <- cond_hist_names(K)[chain_pos(K)]
    
    colnames(chain2_draws) <- cond_hist_names(K)[chain_pos(K)]
    
    colnames(pi_draws) <- c("pi1", "pi2")
    
    colnames(N_draws) <- c("N1", "N2")
    
    
    keep <- (seq_along(1:(B+1)) > (P + 1) & seq_along(1:(B + 1))%%k == 0)
    
    structure(list(chain1 = chain1_draws[keep ,], chain2 = chain2_draws[keep,],
                   pi_draws = pi_draws[keep,], N_draws = N_draws[keep, ],
                   carry = carry, data = data),
              class = "carryover_bayes")
  
}

# tidy.carryover_bayes -------------------------------------------------------------

#This function returns posterior median/mean and equal tailed credible intervals in the
#same columns as tidy.carryover_fit.  Takes a carryover_bayes object as input.

tidy.carryover_bayes <- function(x, level = 0.95, ...){
  res <- data.frame(
    target = c("N1", "N2"),
    method = c(rep("bayes_median", 2)),
    est = rep(0, 2),
    se  = rep(0,2),
    lower = rep(0,2),
    upper = rep(0,2),
    converged = rep(NA, 2),
    boundary = rep(NA, 2),
    truncated = rep(NA, 2)
  )

  ci1 <- unname(quantile(x$N_draws[,1], probs = c((1 - level)/2, level + (1 - level)/2)))
  
  ci2 <- unname(quantile(x$N_draws[,2], probs = c((1 - level)/2, level + (1 - level)/2)))
  
  
  res$est       <- unname(c(apply(x$N_draws, 2, median)))
  res$se        <- unname(c(sd(x$N_draws[,1]), sd(x$N_draws[,2])))
  res$lower     <- c(ci1[1], ci2[1])
  res$upper     <- c(ci1[2], ci2[2])
  
  return(res)
  
}