########################Michael Doerfler########################################
########################2026-10-14##############################################
#Here I include some tests for the 09_sim_engine script.  In particular, what we are testing
#is whether the seeding makes the bayes results reproducible.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

sim <- readRDS(file.path("sims", "k2_asymptotics_pilot", "data", "scenario_001.rds"))
sim$datasets <- sim$datasets[1:5]  # small, for speed

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

a <- fit_scenario(sim, methods)
b <- fit_scenario(sim, methods)
identical(a$est[a$fit_name == "bayes"], b$est[b$fit_name == "bayes"])     # TRUE: reproducible

c2 <- fit_scenario(sim, methods[c("bayes", "mle")])                        # different method order
identical(a$est[a$fit_name == "bayes"], c2$est[c2$fit_name == "bayes"])   # TRUE: order doesn't matter
