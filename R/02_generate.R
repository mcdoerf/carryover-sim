#################Michael Doerfler#################
##################2026-10-08######################
#This script contains functions for simulating R data.

# Simulating multinomial -----------------------------

#sim_multinomial(theta, N, N_tot, psi).  theta = list(theta1, theta2), N = list(N1, N2),
#N_tot = list(N_tot1, N_tot2)
#Output a cr_data object.

sim_multinomial <- function(theta1, theta2, N1, N2, N_tot1, N_tot2, psi){
  if (N1 > N_tot1){
    stop("N1 cannot be greater than N_tot1")
  }
  
  if (N2 > N_tot2){
    stop("N2 cannot be greater than N_tot2")
  }
  
  if (!(psi >0 && psi <1)){
    stop("psi must be in (0,1)")
  }
  
  if(!(length(theta1) == length(theta2))){
    stop("theta1 and theta2 must be of same length")
  }
  
  #Extract K: length(theta1) = 2^{K}-1, so K = log_2(length(theta1) + 1)
  
  K = log(length(theta1) + 1, base = 2)
  
  #Construct cell probabilities for cohort 1 (ph1) and cohort 2 (ph2).
  
  ph1 <- cell_probs(theta1, K, psi = psi)
  
  ph2 <- cell_probs(theta2, K, psi = NULL)
  
  #Simulate from multinomials
  
  counts1 <- setNames(rmultinom(1, N1, prob = ph1)[,1], names(ph1))
  
  counts2 <- setNames(rmultinom(1, N2, prob = ph2)[,1], names(ph2))
  
  counts1 <- counts1[-length(counts1)]     # drop "0...0" (unobserved)
  
  counts2 <- counts2[-length(counts2)]
  
  n_dbar  <- rbinom(1, N_tot1 - N1, prob = psi)
  
  #Return cr_data object

  return(cr_data(counts1, counts2, psi, N_tot1, N_tot2, n_dbar))
  
}


# Simulating stratified --------------------------

#sim_stratified(theta1_mat, theta2_mat, N1_vec, N2_vec, N_tot1_vec, N_tot2_vec, psi) 
#theta1_mat and theta2_mat are S by 2^{K} -1, where S it the number of strata.
#N1_vec is of length S.  N2_vec is of length S. 
#N_tot1_vec is of length S, and N_tot2_vec is of length S.
#Returns a list of cr_data objects, one per strata.

sim_stratified <- function(theta1_mat, theta2_mat, N1_vec, N2_vec, N_tot1_vec, N_tot2_vec, psi){

  lengths <- c(nrow(theta1_mat), nrow(theta2_mat), length(N1_vec), 
               length(N2_vec), length(N_tot1_vec), length(N_tot2_vec))
  
  if (!all(lengths == lengths[1])) {
    stop("theta1_mat, theta2_mat, N1_vec, N2_vec, N_tot1_vec, and N_tot2_vec must have same number of rows")
  }
  
  if(!(ncol(theta1_mat) == ncol(theta2_mat))){
    stop("theta1_mat and theta2_mat must have same number of columns")
  }
  
  if(!(psi > 0 && psi < 1)){
    stop("psi must be in (0,1)")
  }
  
  if(any(N1_vec > N_tot1_vec) || any(N2_vec > N_tot2_vec)){
    stop("Number of diseased in a strata cannot exceed total size of the strata")
  }
  
  S <- nrow(theta1_mat)
  
  cr_data_strata <- vector(mode = "list", length = S)
  
  for (s in 1:S){
    cr_data_strata[[s]] <- sim_multinomial(theta1_mat[s,], theta2_mat[s,], N1_vec[s], N2_vec[s], N_tot1_vec[s], N_tot2_vec[s], psi)
  }
  
  return(cr_data_strata)
}
