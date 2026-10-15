###################################Michael Doerfler#############################
###################################2026-10-14###################################
#This is the 2_fit R script in the k2_asymptotics simulation study.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

methods <- list(
  
  # Carry-over MLE: Wald and logit intervals (4 rows: N1/N2 x wald/logit)
  mle = function(d, row) {
    tidy(fit_carryover(d, carry = row$carry[[1]], method = "closed"))
  },
  
  # Carry-over Bayes (2 rows: N1, N2)
  bayes = function(d, row) {
    tidy(fit_bayes(d, carry = row$carry[[1]], n_keep = 4000, burn = 1000, thin = 1))
  },
  
  # Chapman in cohort 2, streams 1 vs 2 (2 rows: sadinle, waring)
  chapman_c2 = function(d, row) {
    tidy(fit_chapman(d, cohort = 2, streams = list(1, 2)))
  },
  
  # Chapman in cohort 1, all K streams merged vs the anchor (2 rows)
  chapman_c1_anchor = function(d, row) {
    tidy(fit_chapman(d, cohort = 1, streams = list(1:row$K, row$K + 1)))
  }
)


run_study("./sims/k2_asymptotics", stage = "fit", methods = methods, cores = 10)
