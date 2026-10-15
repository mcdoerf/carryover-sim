##############################Michael Doerfler##################################
#############################2026-10-12#########################################
#This R script contains functions for competing inferential methods to compare to 
#the carry-over methods.

# fit_chapman ------------------------------------------------------------------
#This function returns the Chapman estimator and its usual variance.  The function takes
# a cr_data function as input, as well as  cohort (1 or 2) indicating which cohort we are 
# estimating N for.  Finally, streams is a list of two integer vectors, where streams within
# a set are merged into one list.

fit_chapman <- function(data, cohort = 2, streams = list(1, 2)) {
  K <- data$K
  L <- if (cohort == 1) K + 1 else K
  stopifnot(
    "data must be a cr_data object" = inherits(data, "cr_data"),
    "cohort must be 1 or 2"         = cohort %in% c(1, 2),
    "streams must be a list of two sets" = is.list(streams) && length(streams) == 2,
    "stream positions out of range" = all(unlist(streams) %in% seq_len(L)),
    "stream sets must not overlap"  = length(intersect(streams[[1]], streams[[2]])) == 0
  )
  
  counts <- data$counts[[cohort]]
  nc     <- sum(counts)
  U      <- if (cohort == 1) data$N_tot[1] - data$n_dbar else data$N_tot[2]
  
  H   <- do.call(rbind, strsplit(names(counts), "")) == "1"
  inA <- rowSums(H[, streams[[1]], drop = FALSE]) > 0
  inB <- rowSums(H[, streams[[2]], drop = FALSE]) > 0
  
  n11 <- sum(counts[inA & inB])
  n10 <- sum(counts[inA & !inB])
  n01 <- sum(counts[!inA & inB])
  n1. <- n11 + n10
  n.1 <- n11 + n01
  
  N_tilde <- (n1. + 1) * (n.1 + 1) / (n11 + 1) - 1
  N_hat <- min(max(N_tilde, nc), U)
  var_chap <- (n1. + 1) * (n.1 + 1) * n10 * n01 / ((n11 + 1)^2 * (n11 + 2))
  
  truncated_lower = N_tilde < nc
  truncated_upper = N_tilde > U
  truncated       = N_tilde < nc || N_tilde > U
  
  structure(list(
    N_hat = N_hat, N_tilde = N_tilde, truncated = truncated, var = var_chap,
    n1. = n1., n.1 = n.1, n11 = n11, n10 = n10, n01 = n01,
    nc = nc, U = U, cohort = cohort, streams = streams,
    label = paste0(paste(streams[[1]], collapse = ""), "v", paste(streams[[2]], collapse = "")),
    data = data
  ), class = "chapman_fit")
}

# confint.chapman_fit() ---------------------------------------------------------------
#This function takes a chapman_fit object and returns either a Sadinle or a Waring interval.

confint.chapman_fit <- function(object, parm = NULL, level = 0.95,
                                type = c("sadinle", "waring"), ...) {
  type  <- match.arg(type)
  alpha <- 1 - level
  n11 <- object$n11; n10 <- object$n10; n01 <- object$n01
  n2  <- n11 + n10 + n01            # two-list total (formulas)
  nc  <- object$nc                  # all observed (clipping)
  U   <- object$U
  
  if (type == "sadinle") {
    h_n      <- (n10 + 0.5) * (n01 + 0.5) / (n11 + 0.5)
    sigma_tl <- sqrt(1/(n11 + 0.5) + 1/(n10 + 0.5) + 1/(n01 + 0.5) +
                       (n11 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
    z   <- qnorm(1 - alpha/2)
    lwr <- (n2 - 0.5) + h_n * exp(-z * sigma_tl)
    upr <- (n2 - 0.5) + h_n * exp( z * sigma_tl)
  } else {
    lwr <- n2 + qwar(alpha/2,     n10 + 1, n01 + 1, n2 + 3)
    upr <- n2 + qwar(1 - alpha/2, n10 + 1, n01 + 1, n2 + 3)
  }
  
  lwr <- max(lwr, nc); upr <- min(upr, U)
  if (upr <= lwr) { lwr <- nc; upr <- U }   # rare degenerate case
  
  ci <- cbind(lwr, upr)
  dimnames(ci) <- list(paste0("N", object$cohort),
                       paste(format(100 * c(alpha/2, 1 - alpha/2)), "%"))
  ci
}


# Helper functions for Waring interval ------------------------------------------------------
## R Code for Rivest CI to accompany Chapman's estimator ----
##The code for pwar and qwar are adapted from the supplementary materials of Rivest and Yuack 2025.
##Rivest, L.P., Yauck, M., 2025. Small sample inference for two-way capture-recapture experiments. International Statistical Review 93, 62–72. doi:https://doi.org/10.1111/insr.12574.

### R code for the generalized Waring cumulative distribution function ----

pwar <- function(q, a, b, c, lower.tail = TRUE) {
  
  f <- function(x, q) {
    pnbinom(q, size = a, prob = x, lower.tail = lower.tail) *
      dbeta(x, shape1 = c - b - a, shape2 = b)
  }
  
  if (length(q) == 1) {
    integrate(f, lower = 0, upper = 1,
              q = q,
              rel.tol = 1e-10,
              abs.tol = 0)$value
  } else {
    sapply(q, function(qi)
      integrate(f, lower = 0, upper = 1,
                q = qi,
                rel.tol = 1e-10,
                abs.tol = 0)$value)
  }
}




### R code for the generalized Waring quantile function ----    
# p is a probability in (0,1)
# a b and c are positive real parameters satisfying c>a+b


qwar <- function(p, a, b, c) {
  
  if (p <= 0) return(0)
  
  # Step 1: Find an upper bound
  upper <- 1
  while (pwar(upper, a, b, c) < p) {
    upper <- upper * 2
  }
  
  lower <- 0
  
  # Step 2: Bisection search
  while (upper - lower > 1) {
    mid <- floor((lower + upper) / 2)
    if (pwar(mid, a, b, c) >= p) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  return(upper)
}


# tidy.chapman_fit() -----------------------------------------------------------
#This function takes an object x of type chapman_fit and returns a summary.

tidy.chapman_fit <- function(x, level = 0.95, ...) {
  ci_s <- confint(x, level = level, type = "sadinle")
  ci_w <- confint(x, level = level, type = "waring")
  
  data.frame(
    target    = paste0("N", x$cohort),
    method    = paste0("chapman_", x$label, c("_sadinle", "_waring")),
    est       = unname(x$N_hat),
    se        = sqrt(unname(x$var)),
    lower     = c(ci_s[, 1], ci_w[, 1]),
    upper     = c(ci_s[, 2], ci_w[, 2]),
    converged = NA,
    boundary  = x$n11 == 0,
    truncated = unname(x$truncated)
  )
}
