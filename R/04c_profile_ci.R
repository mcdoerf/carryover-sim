# profile_ci_full --------------------------------------------------------------
# Profile-likelihood interval for N1 or N2 from a carryover_fit_full object.
# For each candidate value x of the target N, maximize the profile log-likelihood
# over the other N; the interval is {x : 2 * (loglik_max - l_p(x)) <= qchisq(level, 1)}.
# If the likelihood has not dropped far enough at a bound, that bound is the endpoint.

profile_ci_full <- function(object, which = c("N1", "N2"), level = 0.95){
  which <- match.arg(which)
  j     <- match(which, c("N1", "N2"))
  other <- 3 - j
  
  data  <- object$data
  fc    <- first_capture_counts(data)
  b <- fc$b; a <- fc$a
  lower <- c(sum(b), sum(a))
  upper <- c(data$N_tot[1] - data$n_dbar, data$N_tot[2])
  
  ll <- function(N) full_profile_loglik(N, b = b, a = a, psi = data$psi, K = data$K,
                                        carry = object$carry, N_tot1 = data$N_tot[1],
                                        n_dbar = data$n_dbar)
  ll_max <- ll(object$N_cont)           # maximum of the (continuous) profile
  crit   <- qchisq(level, 1) / 2
  
  # profile out the other N at a fixed value x of the target
  lp <- function(x) {
    optimize(function(y) { N <- numeric(2); N[j] <- x; N[other] <- y; ll(N) },
             interval = c(lower[other], upper[other]), maximum = TRUE, tol = 1e-8)$objective
  }
  g <- function(x) ll_max - lp(x) - crit        # 0 at the endpoints, < 0 inside the interval
  
  N_j <- object$N_cont[[j]]
  lo <- if (g(lower[j]) > 0) uniroot(g, c(lower[j], N_j), tol = 1e-6)$root else lower[j]
  hi <- if (g(upper[j]) > 0) uniroot(g, c(N_j, upper[j]), tol = 1e-6)$root else upper[j]
  c(lower = lo, upper = hi)
}


# confint.carryover_fit_full ---------------------------------------------------
# type = "wald" or "logit" -> same as the conditional fit (inherited method);
# type = "profile"         -> profile-likelihood interval.

confint.carryover_fit_full <- function(object, parm = c("N1", "N2"), level = 0.95,
                                       type = c("wald", "logit", "profile"), ...){
  type <- match.arg(type)
  if (type != "profile") return(confint.carryover_fit(object, parm = parm, level = level, type = type, ...))
  
  ci <- t(vapply(c("N1", "N2"), function(w)
    tryCatch(profile_ci_full(object, which = w, level = level),
             error = function(e) c(lower = NA_real_, upper = NA_real_)),
    numeric(2)))
  dimnames(ci) <- list(c("N1", "N2"), paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), "%"))
  ci[parm, , drop = FALSE]
}


# tidy.carryover_fit_full ------------------------------------------------------
# Rows: full_wald, full_logit (relabelled from the conditional tidy method) and,
# if profile = TRUE, full_profile. The se column of the profile rows repeats the
# Wald SE, for reference only.

tidy.carryover_fit_full <- function(x, level = 0.95, profile = TRUE, ...){
  out <- tidy.carryover_fit(x, level = level, ...)
  out$method <- sub("^mle_", "full_", out$method)
  if (!profile) return(out)
  
  prof <- out[out$method == "full_wald", ]
  ci   <- confint(x, level = level, type = "profile")
  prof$method <- "full_profile"
  prof$lower  <- ci[prof$target, 1]
  prof$upper  <- ci[prof$target, 2]
  rbind(out, prof)
}
