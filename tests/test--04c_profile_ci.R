#####################################Michael Doerfler###########################
#####################################2026-20-10#################################
#This script is for testing the 04c_profile_ci script.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}


toy3 <- structure(
  list(
    K   = 3,
    psi = 0.1,
    counts = list(
      c("1111" = 6,  "1110" = 24, "1101" = 10, "1100" = 40,     # cohort 1: 2^(K+1) - 1 = 15 cells
        "1011" = 8,  "1010" = 32, "1001" = 12, "1000" = 48,
        "0111" = 5,  "0110" = 25, "0101" = 10, "0100" = 40,
        "0011" = 10, "0010" = 40, "0001" = 20),
      c("111" = 20, "110" = 40, "101" = 15, "100" = 75,          # cohort 2: 2^K - 1 = 7 cells
        "011" = 10, "010" = 50, "001" = 40)
    ),
    N_tot  = c(5000, 8000),
    n_dbar = 450,
    stream_names = c("pathology", "hospital", "death_cert")
  ),
  class = "cr_data"
)

# Testing profile_ci_full-------------------------------------------------------
object = fit_carryover_full(data = toy3, carry = c(0,0,1))
object_bayes = fit_bayes(data = toy3, carry = c(0,0,1))

#If things are working correctly, this should give a similar result to the Bayes interval.
profile_ci_full(object = object, which = "N1", level = 0.95)
profile_ci_full(object = object, which = "N2", level = 0.95)
tidy(object_bayes)

# Testing confint---------------------------------------------------------------

#Should also be similar, but the profile interval should be moving rightward to account
#for skewness of estimator.

#Based on hand calcs, the wald should be [423.7, 583.87]
confint(object = object, parm = "N1", carry = c(0,0,1), type = "wald")
confint(object = object, parm = "N1", carry = c(0,0,1), type = "logit")
confint(object = object, parm = "N1", carry = c(0,0,1), type = "profile")

#Based on hand calcs, the wald should be [307.1 ,497.6]
confint(object = object, parm = "N2", carry = c(0,0,1), type = "wald")
confint(object = object, parm = "N2", carry = c(0,0,1), type = "logit")
confint(object = object, parm = "N2", carry = c(0,0,1), type = "profile")

