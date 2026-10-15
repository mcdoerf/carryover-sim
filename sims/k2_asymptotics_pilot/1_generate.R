########################################Michael Doerfler########################
########################################2026-10-13##############################
#This is the 1_generate script for the k2_asymptotics_pilot simulation study, in which the 
#data are generated.

library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

scenarios <- readRDS("./sims/k2_asymptotics_pilot/design.rds")

run_study("./sims/k2_asymptotics_pilot", stage = "generate", scenarios = scenarios$scenario_id, R = 200,
          base_seed = 1, cores = 10)
