############################Michael Doerfler####################################
###########################2026-10-12###########################################
#This contains code for conducting the conditional likelihood ratio test of additional
#carry-over constraints. carry_alt is a vector of length K of 0's and 1's indicating which
#parameters in the non-identifiable chain to carry over.  carry_null is a vector of length
#K of 0's and 1's indicating which parameters in the non-identifiable chain carry-over.
#carry_null must have 1's at least in every place carry_alt has 1's.

lrt_carryover <- function(data, carry_null, carry_alt) {
  stopifnot(
    "carry vectors must have the same length" = length(carry_null) == length(carry_alt),
    "carry_null must contain carry_alt (nested)" = all(carry_null >= carry_alt),
    "carry_alt must carry over at least one parameter" = sum(carry_alt) >= 1
  )
  df <- sum(carry_null) - sum(carry_alt)
  stopifnot("models are identical (df = 0)" = df > 0)
  
  fit_null <- fit_carryover(data = data, carry = carry_null, method = "optim")
  fit_alt  <- fit_carryover(data = data, carry = carry_alt,  method = "optim")
  
  if (!fit_null$converged || !fit_alt$converged) warning("a fit did not converge")
  if (any(fit_null$boundary, fit_alt$boundary))  warning("estimate on boundary; chi-square reference may be poor")
  
  stat <- max(2 * (fit_alt$loglik - fit_null$loglik), 0)
  
  list(stat = stat, df = df,
       p_value = pchisq(stat, df = df, lower.tail = FALSE),
       carry_null = carry_null, carry_alt = carry_alt)
}


