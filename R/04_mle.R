#########################Michael Doerfler###############################
########################2026-10-09######################################
#This script contains code for calculating conditional MLEs.

# chain_loglik()------------------------------
#This function takes par: free chain parameters on logit scale, where it is ordered as
#theta1^{(1)}, theta2bar1^{(1)},...,thetaKbar12...K-1^{(1)}, cohort-2-only parameters.
#carry is a vector of 0's and 1's of length K indicating which chain parameters to carry-over.

chain_loglik <- function(par, fc, psi, K, carry){
  stopifnot(
    "K must be a single positive integer" =
      is.numeric(K) && length(K) == 1 && K >= 1 && K == round(K),
    "carry must be a 0/1 vector of length K with at least one 1" =
      is.numeric(carry) && length(carry) == K && all(carry %in% c(0, 1)) && any(carry == 1),
    "par must be a numeric vector of length K + sum(carry == 0), with no NAs" =
      is.numeric(par) && length(par) == K + sum(carry == 0) && !anyNA(par),
    "psi must be a single number in (0, 1)" =
      is.numeric(psi) && length(psi) == 1 && psi > 0 && psi < 1,
    "fc must be a list with elements b (length K + 1) and a (length K)" =
      is.list(fc) && length(fc$b) == K + 1 && length(fc$a) == K,
    "b and a in fc must be non-negative integer counts" =
      is.numeric(fc$b) && is.numeric(fc$a) &&
      !anyNA(c(fc$b, fc$a)) && all(c(fc$b, fc$a) >= 0) &&
      all(c(fc$b, fc$a) == round(c(fc$b, fc$a)))
  )
  
  betaU <- plogis(par)
  
  betaU1 <- betaU[1:K]
  
  betaU2 <- betaU[-(1:K)]
  
  #Extracting b's and a's from first counts.
  
  b <- fc$b
  
  a <- fc$a
  
  #Extracting sums (sum_{i=j+1}^{K+1} b_i) and (sum_{i=j+1}^{K} a_i)
  
  b_partial_sums <- sum(b) - cumsum(b)
  
  b_partial_sums <- head(b_partial_sums, -1)
  
  a_partial_sums <- sum(a) - cumsum(a)
  
  #Extracting numbers captured
  
  nc1 <- sum(b)
  
  nc2 <- sum(a)

  #Extracting capture probability
  
  pc1 <- 1 - (1 - psi) * prod(1 - betaU1)
  
  pc2 <- 1 - prod(1 - betaU1[carry == 1]) * prod(1 - betaU2)
  
  #Putting together joint log-likelihood
  
  lU1 <- sum(head(b, -1) * log(betaU1) + b_partial_sums * log(1 - betaU1)) +
         sum( (a * log(betaU1) + a_partial_sums * log(1 - betaU1)) * carry) -
         nc1 * log(pc1)
  
  lU2 <- sum(a[carry == 0] * log(betaU2) + a_partial_sums[carry == 0] * log(1 - betaU2))- 
         nc2 * log(pc2)
  
  #Return the joint conditional log-likelihood
  
  return(lU1 + lU2)
  
}


##gradient_chain_loglik()------------------------------------------------------------
#This function takes par: free chain parameters on logit scale, where it is ordered as
#theta1^{(1)}, theta2bar1^{(1)},...,thetaKbar12...K-1^{(1)}, cohort-2-only parameters.
#carry is a vector of 0's and 1's of length K indicating which chain parameters to carry-over.
#Returns the gradient of the joint log likelihood.

gradient_chain_loglik <- function(par, fc, psi, K, carry){
  
  betaU <- plogis(par)
  
  betaU1 <- betaU[1:K]
  
  betaU2 <- betaU[-(1:K)]
  
  #Extracting b's and a's from first counts.
  
  b <- fc$b
  
  a <- fc$a
  
  #Extracting sums (sum_{i=j+1}^{K+1} b_i) and (sum_{i=j+1}^{K} a_i)
  
  b_partial_sums <- sum(b) - cumsum(b)
  
  b_partial_sums <- head(b_partial_sums, -1)
  
  a_partial_sums <- sum(a) - cumsum(a)
  
  #Extracting numbers captured
  
  nc1 <- sum(b)
  
  nc2 <- sum(a)
  
  #Extracting capture probability
  
  pc1 <- 1 - (1 - psi) * prod(1 - betaU1)
  
  pc2 <- 1 - prod(1 - betaU1[carry == 1]) * prod(1 - betaU2)
  
  #Putting together gradient
  
  gradU1 <- head(b, -1) / betaU1 - b_partial_sums / (1 - betaU1) - nc1/pc1 * (1 - pc1) / (1 - betaU1) +
            (a / betaU1 - (a_partial_sums / (1 - betaU1)) - nc2 / pc2 * (1 - pc2) / (1 - betaU1)) * carry
  
  gradU2 <- a[carry == 0] / betaU2 - (a_partial_sums[carry == 0]) / (1 - betaU2) -
            nc2/pc2 * (1 - pc2) / (1 - betaU2)
  
  #Return the gradient.  Do not forget to multiply by betaU * (1 - betaU) for the chain rule.
  
  return(c(gradU1, gradU2) * betaU * (1 - betaU))
  
}



# unpack_chain(par, K, carry)--------------------------------------------------------------

#This function maps the free parameter chain vector betaU to list(chain1, chain2), inserting
#shared values in both cohorts.

unpack_chain<- function(betaU, K, carry){
  
  stopifnot("betaU must have length K + sum(carry == 0)" =
              length(betaU) == K + sum(carry == 0))
  
  chain1 <- betaU[1:K]
  chain2 <- rep(0, K)
  
  chain2[carry == 1] <- chain1[carry == 1]
  chain2[carry == 0] <- betaU[-(1:K)]
  
  return(list(chain1 = chain1, chain2 = chain2))
}

# chain_mle ----------------------------------------------------------------------------

#This function maximizes the chain log-likelihood.  "closed" uses the closed form (valid when J = K),
# and "optim" uses BFGS on the logit scale.  
#Returns list: chain1, chain2, loglik, converged, method_used, boundary (logical: any estmiate
#with 10^{-6} of 0 or 1).

chain_mle <- function(fc, psi, K, carry, method = c("closed", "optim"), start = NULL,
                      control = list()) {
  method <- match.arg(method)
  if (method == "closed" && !identical(as.numeric(carry), c(rep(0, K - 1), 1)))
    stop('method = "closed" is only available when carrying over theta_{K|0...0} alone; use method = "optim"')
  
  if (method == "closed") {
    b <- fc$b
    a <- fc$a
    
    N1_tilde   <- sum(b) - b[K + 1] + b[K + 1] / psi
    chain1_hat <- b[1:K] / (N1_tilde - c(0, head(cumsum(b), -2)))
    
    N2_tilde   <- sum(a) - a[K] + a[K] / chain1_hat[K]
    chain2_hat <- a / (N2_tilde - c(0, head(cumsum(a), -1)))
    
    betaU     <- c(chain1_hat, chain2_hat[carry == 0])
    converged <- TRUE
  }
  
  if (method == "optim") {
    n_par <- K + sum(carry == 0)
    if (is.null(start)) start <- rep(0.5, n_par)
    stopifnot("start must have length K + sum(carry == 0), with values in (0, 1)" =
                length(start) == n_par && all(start > 0 & start < 1))
    
    opt_result <- optim(
      par     = qlogis(start),                 # probability -> logit scale
      fn      = chain_loglik,
      gr      = gradient_chain_loglik,
      method  = "BFGS",
      control = modifyList(list(fnscale = -1, maxit = 1000, reltol = 1e-12), control),
      fc = fc, psi = psi, K = K, carry = carry  # passed on to fn and gr
    )
    betaU     <- plogis(opt_result$par)
    converged <- opt_result$convergence == 0
  }
  
  par_hat <- qlogis(betaU)
  c(unpack_chain(betaU, K, carry),
    list(loglik      = chain_loglik(par_hat, fc, psi, K, carry),
         converged   = converged,
         method_used = method,
         boundary    = any(betaU < 1e-6 | betaU > 1 - 1e-6)))
}


# ident_mle------------------------------------------------------
##This function provides the closed-form estimates for the identifiable theta parameters.
#Used for reporting and for the Bayesian starting values, but not needed for \hat{N^{(t)}}.

ident_mle <- function(data){
  K <- data$K
  
  id1 <- rep(0, 2^K - 1)
  
  id2 <- rep(0, 2^K - 1)
  
  counts1 <- data$counts[[1]]
  
  counts2 <- data$counts[[2]]
  
  for (j in 2:K){
    for (f in hist_names(j - 1, observable = TRUE)){
      
    num1 <- sum(counts1[startsWith(names(counts1), paste0(f, 1))])
    
    den1 <- num1 + sum(counts1[startsWith(names(counts1), paste0(f, 0))])
    
    num2 <- sum(counts2[startsWith(names(counts2), paste0(f, 1))])
    
    den2 <- num2 + sum(counts2[startsWith(names(counts2), paste0(f, 0))])
    
    id1[theta_index(j, f)] <- num1 / den1
    
    id2[theta_index(j, f)] <- num2 / den2
    
    
      
    }
  }
  
  return(list(
    hat_thetaI1 =  id1[-(chain_pos(K))],
    hat_thetaI2 = id2[-(chain_pos(K))]
  ))
  
}


# fit_carryover() --------------------------------------------------------------
#This function takes a cr_data object as input, along with a carry vector of length K of 
#0's and 1's that indicated which parameters in the non-identifiable chain to carryover.

##The functoin validates input, computes first-capture counts, the chain MLE, and the
##asymptotic covariance.

fit_carryover <- function(data, carry, method = c("closed", "optim"), start = NULL, control = list()){
  #extract relevant parameters from the data
  
  psi <- data$psi
  
  K <- data$K
  
  #compute first capture counts
  
  fc <- first_capture_counts(data)
  
  #extract MLEs of nonidentifiable chain parameters
  
  chain_est <- chain_mle(fc = fc, psi = psi, K = K, carry = carry, method = method,
                         start = start, control = control)
  
  #Compute N_tilde
  
  N1_tilde <- sum(data$counts[[1]]) / (1 - (1 - psi) * prod(1 - chain_est$chain1))
  
  N2_tilde <- sum(data$counts[[2]]) / (1 - prod(1 - chain_est$chain2))
  
  N_tilde <- c(N1_tilde, N2_tilde)
  
  #Compute the truncated N_hat
  
  U <- c(data$N_tot[1] - data$n_dbar, data$N_tot[2])
  
  N_hat <- pmin(N_tilde, U)
  
  #Compute whether truncation occured
  
  truncated <- N_tilde > U
  
  #Compute the asymptotic covariance matrix
  
  Sigma_NN <- sigma_NN(chain1 = chain_est$chain1, chain2 = chain_est$chain2,
                       N1 = N_hat[1], N2 = N_hat[2], K = K, psi = psi, carry = carry)
  
  vcov_N <- outer(sqrt(N_hat), sqrt(N_hat)) * Sigma_NN
  
  structure(list(
    N_hat = N_hat,
    N_tilde = N_tilde,
    truncated = truncated,
    chain1 = chain_est$chain1,
    chain2 = chain_est$chain2,
    Sigma_NN = Sigma_NN,
    vcov_N = vcov_N,
    loglik = chain_est$loglik,
    converged = chain_est$converged,
    boundary = chain_est$boundary,
    method_used = chain_est$method_used,
    carry = carry,
    data = data
  ), class = "carryover_fit")
  
  
}


# confint.carryover_fit() ------------------------------------------------------
#This function takes a carryover_fit object and returns a Wald interval or logit-transformed
#interval as defined in the paper.  Lower limits are truncated at nc^{(t)} and upper limits
#at N^{(t)}.

confint.carryover_fit <- function(object, parm = c("N1", "N2"), level = 0.95, type = c("wald", "logit"), ...){
  type  <- match.arg(type)
  z     <- qnorm(1 - (1 - level) / 2)
  N_hat <- object$N_hat
  se    <- sqrt(diag(object$vcov_N))
  n_c   <- vapply(object$data$counts, sum, numeric(1))
  U     <- c(object$data$N_tot[1] - object$data$n_dbar, object$data$N_tot[2])
  N_tot <- object$data$N_tot
  
  if(type == "wald"){
    lower <- N_hat - z * se
    upper <- N_hat + z * se
  } else{
    var_logit <- 1 / ((N_hat / N_tot)^2 * (1 - N_hat / N_tot)^2) * 1 / (N_tot^2) * se^2
    lower <- N_tot * (plogis(qlogis(N_hat / N_tot) - z * sqrt(var_logit)))
    upper <- N_tot * (plogis(qlogis(N_hat / N_tot) + z * sqrt(var_logit)))
    if(object$truncated[2] == TRUE){
      #If hat(N2) = Ntot2, the logit interval falls apart for N2.  Fall back on the 
      #Wald interval in that case.
      lower[2] <- N_hat[2] - z * se[2]
      upper[2] <- N_hat[2] + z * se[2]
    }
  }
  
  # truncate to the plausible range
  lower <- pmax(lower, n_c)
  upper <- pmin(upper, U)
  
  ci <- cbind(lower, upper)
  dimnames(ci) <- list(c("N1", "N2"),
                       paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), "%"))
  ci[parm, , drop = FALSE]
  
}


# tidy.carryover_fit()----------------------------------------------------------
#This function will take a fit_carryover object x and output a dataframe for each target
# (N1 and N2) along with a point est, se, ci, (by Wald and Logit).

tidy <- function(x, ...) UseMethod("tidy")

tidy.carryover_fit <- function(x, level = 0.95,...){
  res <- data.frame(
    target = c("N1", "N2", "N1", "N2"),
    method = c(rep("mle_wald", 2), rep("mle_logit", 2)),
    est = rep(0, 4),
    se  = rep(0,4),
    lower = rep(0,4),
    upper = rep(0,4),
    converged = rep(FALSE, 4),
    boundary = rep(FALSE, 4),
    truncated = rep(FALSE, 4)
  )
  
  ci_wald <- confint(object = x, level = level, type = "wald")
  
  ci_logit <- confint(object = x, level = level, type = "logit")
  
  
  res$est       <- rep(unname(x$N_hat), 2)
  res$se        <- rep(unname(sqrt(diag(x$vcov_N))), 2)
  res$lower     <- c(ci_wald[, 1], ci_logit[, 1])
  res$upper     <- c(ci_wald[, 2], ci_logit[, 2])
  res$converged <- x$converged
  res$boundary  <- x$boundary
  res$truncated <- rep(unname(x$truncated), 2)
  
  return(res)
  

}



