########################Michael Doerfler###################
######################2026-10-08###########################
#This script contains a function for calculating the theoretical \Sigma_{NN} at the true parameters.


#1) true_sigma_NN()--------------------------------------------------------------
#This function takes true or estimated population level chain1, chain2, N1, N2, N_tot1,
#N_tot2, and psi as input and outputs the asymptotic covariance matrix for 1/sqrt{N_1} (\hat{N}_1 - N_1),
#1/\sqrt{N_2} (\hat{N}_2 - N_2).
# carry is a vector of length K indicating which non-identifiable parameters carry-over from population 1.

#2026-10-20 update: added an e argument as well to sigma_NN.  The e is for computing
#sigma_NN for the full likleihood MLE, 
# e: extra information for N1 from the anchor's non-cases (full likelihood only),
#    e = rho * psi / ((1 - rho) * (1 - psi)), rho = N1 / N_tot1.
#    In the N1^{-1/2}(N1_hat - N1) scale used by D. Use e = 0 for the conditional fit.
sigma_NN <- function(chain1, chain2, N1, N2, K, psi, carry, e = 0){
  
  stopifnot(
    "K must be a single positive integer" =
      is.numeric(K) && length(K) == 1 && K >= 1 && K == round(K),
    "chain1 and chain2 must be numeric vectors of length K" =
      is.numeric(chain1) && is.numeric(chain2) &&
      length(chain1) == K && length(chain2) == K,
    "chain1 and chain2 must lie strictly between 0 and 1 (no NAs)" =
      !anyNA(c(chain1, chain2)) && all(chain1 > 0 & chain1 < 1) && all(chain2 > 0 & chain2 < 1),
    "N1 and N2 must be single positive numbers" =
      is.numeric(N1) && is.numeric(N2) && length(N1) == 1 && length(N2) == 1 && N1 > 0 && N2 > 0,
    "psi must be a single number in (0, 1)" =
      is.numeric(psi) && length(psi) == 1 && psi > 0 && psi < 1,
    "carry must be a 0/1 vector of length K with at least one 1" =
      length(carry) == K && all(carry %in% c(0, 1)) && any(carry == 1),
    "e must be a single non-negative number" =
      is.numeric(e) && length(e) == 1 && is.finite(e) && e >= 0
  )
  
  shared <- carry == 1
  stopifnot(
    "carried-over chain parameters must be equal in theta1 and theta2" =
      isTRUE(all.equal(chain1[shared], chain2[shared]))
  )
  
  lambda1 <- N1 / (N1 + N2)
  
  lambda2 <- 1 - lambda1
  
  #Extract the non-identifiable chain parameters
  
  betaU <- c(chain1, chain2[carry == 0])
  
  dU    <- length(betaU)
  
  #Create cumulative products 
  cum_chain1 <- head(c(1, cumprod(1 - chain1)), -1)
  
  cum_chain2 <- head(c(1, cumprod(1 - chain2)), -1)
  
  # Lengthen carry vector to length of dU
  
  carry <- c(carry, rep(0, dU - K))
  
  #Create a vector of indicators for whether parameter in betaU belongs only to pop2
  
  pop2_only <- c(rep(0, K), rep(1, dU - K))
  
  #Calculate probability of being missed by all K+1 streams in cohort 1, and all K streams in cohort 2
  
  p_miss1 <- (1 - psi) * prod(1 - chain1)
  
  p_miss2 <- prod(1 - chain2)
  
  #Calculate the D matrix
  
  D <- diag(c((1 - p_miss1) / p_miss1, (1 - p_miss2) / p_miss2))
  
  #Calculate the A matrix

  A <- diag(rep(0, dU))
  
  denominator <- c(lambda1 * cum_chain1 + lambda2 * cum_chain2 * carry[1:K],   # cohort-1 and shared entries
                   lambda2 * cum_chain2[carry[1:K] == 0])                      # cohort-2-only entries
  
  diag_values <- (betaU * (1 - betaU)) / denominator

  diag(A) <- diag_values
  
  #Calculate the B matrix
  
  B1 <-(lambda1^{1/2}/(1 - betaU)) * (1 - pop2_only)
  
  B2 <- (lambda2^{1/2} / (1 - betaU)) * (1 - (1 - carry) * (1 - pop2_only))
  
  B <- cbind(B1, B2)
  
  #Calculate and return SigmaNN
  D[1,1] <- D[1,1] + e
  
  Sigma_NN <- solve(D - t(B) %*% A %*% B)
  
  return(Sigma_NN)


}
  
  
  


