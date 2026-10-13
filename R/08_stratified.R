#########################Michael Doerfler##################################
#########################2026-10-12########################################
##This R script contains functions for stratified inerences.

# pool_strata() ----------------------------------------------------------------
#This function pools simulated stratified cr_data into one cr_data object.

pool_strata <- function(data_list){
  stopifnot(
    "data_list must be a list of cr_data objects" =
      is.list(data_list) && all(vapply(data_list, inherits, logical(1), "cr_data")),
    "all strata must have the same K" =
      length(unique(vapply(data_list, function(d) d$K, numeric(1)))) == 1,
    "all strata must have the same psi" =
      length(unique(vapply(data_list, function(d) d$psi, numeric(1)))) == 1
  )
  counts1 <- Reduce("+", lapply(data_list, function(d) d$counts[[1]]))
  counts2 <- Reduce("+", lapply(data_list, function(d) d$counts[[2]]))
  N_tot1   <- Reduce("+", lapply(data_list, function(d) d$N_tot[[1]]))
  N_tot2 <- Reduce("+", lapply(data_list, function(d) d$N_tot[[2]]))
  n_dbar <- Reduce("+", lapply(data_list, function(d) d$n_dbar))
  
  return(cr_data(counts1 = counts1, counts2 = counts2, psi = data_list[[1]]$psi,
                 N_tot1 = N_tot1, N_tot2 = N_tot2, n_dbar = n_dbar))
}


# fit_stratified() ---------------------------------------------------------------
# data_list is a list of cr_data objects, one per stratum.  carry is a list of  0,1 vectors of length
# K for each stratum, indicating which parameters in the non-identifiable chain to carry-over
# within that stratum.
fit_stratified <- function(data_list, carry_list, method = c("auto", "closed", "optim")) {
  method <- match.arg(method)
  stopifnot("data_list must be a list of cr_data objects" =
              is.list(data_list) && all(vapply(data_list, inherits, logical(1), "cr_data")))
  S <- length(data_list)
  if (!is.list(carry_list)) carry_list <- rep(list(carry_list), S)
  stopifnot("carry_list must have one entry per stratum" = length(carry_list) == S)
  
  data <- pool_strata(data_list)    # first, so its K/psi checks run before fitting
  
  # fit each stratum
  fits <- Map(fit_carryover, data_list, carry_list, MoreArgs = list(method = method))
  
  # combine across strata
  structure(list(
    N_hat     = Reduce("+", lapply(fits, function(f) f$N_hat)),
    N_tilde   = Reduce("+", lapply(fits, function(f) f$N_tilde)),
    vcov_N    = Reduce("+", lapply(fits, function(f) f$vcov_N)),
    truncated = Reduce("|", lapply(fits, function(f) f$truncated)),
    loglik    = Reduce("+", lapply(fits, function(f) f$loglik)),
    converged = all(vapply(fits, function(f) f$converged, logical(1))),
    boundary  = any(vapply(fits, function(f) f$boundary, logical(1))),
    data      = data,
    fits      = fits,
    carry     = carry_list
  ), class = c("carryover_strat", "carryover_fit"))
}


# bayes_stratified()------------------------------------------------------------
#This gets posterior draws within each stratum.  Then, we sum over stratum.

bayes_stratified <- function(data_list, carry_list, n_keep = 2000, burn = 1000, thin = 5, prior = c(0.5, 0.5)){

  stopifnot("data_list must be a list of cr_data objects" =
              is.list(data_list) && all(vapply(data_list, inherits, logical(1), "cr_data")))
  S <- length(data_list)
  if (!is.list(carry_list)) carry_list <- rep(list(carry_list), S)
  stopifnot("carry_list must have one entry per stratum" = length(carry_list) == S)
  
  data <- pool_strata(data_list)
  
  fits <- Map(fit_bayes, data_list, carry_list, MoreArgs = list(n_keep = n_keep, burn = burn, thin = thin, prior = prior))
  
  structure(list(
    N_draws = Reduce("+", lapply(fits, function(f) f$N_draws)),
    carry   = carry_list,
    data = data,
    fits = fits
  ), class = c("carryover_bayes_strat", "carryover_bayes"))
}

