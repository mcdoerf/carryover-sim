#######################Michael Doerfler####################
#####################2026-10-08############################
# 01_data.R.  Contains functions for creating cr_data object.

#cr_data(counts1, counts2, psi, N_tot, n_dbar, stream_names = NULL)
#counts1: named counts for cohort 1's observable histories.
#counts2: named counts for cohort 2's observable histories.
#psi, N_tot, n_dbar and optional stream labels
#Returns cr_data object.

cr_data <- function(counts1, counts2, psi, N_tot1, N_tot2, n_dbar, stream_names = NULL) {
  # --- scalars ---
  if (!(length(psi) == 1 && psi > 0 && psi < 1)) stop("psi must be a single number in (0, 1)")
  is_count <- function(x) all(!is.na(x) & x >= 0 & x == round(x))
  if (!is_count(c(N_tot1, N_tot2, n_dbar))) stop("N_tot1, N_tot2 and n_dbar must be non-negative integers")
  
  # --- counts: names, lengths, values ---
  if (is.null(names(counts1)) || is.null(names(counts2))) stop("counts1 and counts2 must be named")
  K <- nchar(names(counts2)[1])
  if (!setequal(names(counts1), hist_names(K + 1)) || anyDuplicated(names(counts1)))
    stop("names(counts1) must be exactly the 2^(K+1) - 1 observable histories of length K + 1")
  if (!setequal(names(counts2), hist_names(K)) || anyDuplicated(names(counts2)))
    stop("names(counts2) must be exactly the 2^K - 1 observable histories of length K")
  if (!is_count(counts1) || !is_count(counts2)) stop("counts must be non-negative integers")
  
  # --- totals ---
  if (sum(counts1) + n_dbar > N_tot1) stop("n_c^(1) + n_dbar cannot exceed N_tot1")
  if (sum(counts2) > N_tot2)          stop("n_c^(2) cannot exceed N_tot2")
  if (!is.null(stream_names) && length(stream_names) != K) stop("stream_names must have length K")
  
  # --- canonical order ---
  counts1 <- counts1[hist_names(K + 1)]
  counts2 <- counts2[hist_names(K)]
  
  structure(
    list(K = K, psi = psi, counts = list(counts1, counts2),
         N_tot = c(N_tot1, N_tot2), n_dbar = n_dbar, stream_names = stream_names),
    class = "cr_data"
  )
}


#print.cr_data(x) Takes a cr_data object x as input and outputs a readable summary.

print.cr_data <- function(x, ...) {
  cat("<cr_data>  K =", x$K, "  psi =", x$psi, "\n")
  for (t in 1:2) {
    cat("\nCohort", t, " (n_c =", sum(x$counts[[t]]), ", N_tot =", x$N_tot[t], ")\n")
    print(x$counts[[t]])
  }
  cat("\nNon-diseased sampled by anchor (n_dbar):", x$n_dbar, "\n")
  if (!is.null(x$stream_names)) cat("Streams:", paste(x$stream_names, collapse = ", "), "\n")
  invisible(x)
}
