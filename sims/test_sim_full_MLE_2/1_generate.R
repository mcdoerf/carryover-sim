########################################Michael Doerfler########################
########################################2026-10-13##############################
#This is the 1_generate script for the test_sim_full_MLE_2 simulation study, in which the 
#data are generated.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

scenarios <- readRDS("./sims/test_sim_full_MLE_2/design.rds")

run_study("./sims/test_sim_full_MLE_2", stage = "generate", scenarios = scenarios$scenario_id, R = 200,
          base_seed = 2)
