###########################Michael Doerfler#####################################
###########################2026-10-12###########################################
#In this R script, I test the chapman estimator and its corresponding intervals from the
#07_comparators code.
library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

#I create a toy data set.  The toy dataset is designed so that streams 1 and 2 show the 
#same counts as in the 5-cell paper.

toy <- structure(
  list(
    K   = 2,
    psi = 0.1,
    counts = list(
      c("111" = 6, "110" = 6, "101" = 26, "100" = 26,   # cohort 1: 2^(K+1) - 1 = 7 cells
        "011" = 4, "010" = 15, "001" = 25),
      c("11"  = 12, "10"  = 52, "01"  = 19)                # cohort 2: 2^K - 1 = 3 cells
    ),
    N_tot  = c(1029, 1029),
    n_dbar = 169,
    stream_names = c("pathology", "hospital")
  ),
  class = "cr_data"
)

# Testing fit_chapman() --------------------------------------------------------
#If these are working correctly, they should match the results form Table 7 of the 5-cell
#paper.
res1 <- fit_chapman(data = toy, cohort = 1, streams = list(1,2))
tidy(res1)

res2 <- fit_chapman(data = toy, cohort = 2, streams = list(1,2))
tidy(res2)

