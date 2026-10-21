########################################Michael Doerfler########################
########################################2026-10-13##############################
#This is the 1_generate script for the test_sim_full_MLE simulation study, in which the 
#data are generated.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

scenarios <- readRDS("./sims/test_sim_full_MLE/design.rds")

run_study("./sims/test_sim_full_MLE", stage = "generate", scenarios = scenarios$scenario_id, R = 1000,
          base_seed = 1)
