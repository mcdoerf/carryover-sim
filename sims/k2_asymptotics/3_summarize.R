#############################Michael Doerfler###################################
#############################2026-10-14#########################################
#This script contains R code for putting together performance measures for the k2_asymptotics
# simulation study.

library(tidyverse)
library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}


#Load the design grid corresponding to this simulation study.-------------------------------
design <- readRDS("./sims/k2_asymptotics/design.rds")

#Create the truth d.f. ---------------------------------------------------------------------
truth_res <- truth_generator(design = design)

#Extract the fits objects and row bind them into one large object.--------------------------

# 1. Get a vector of all .rds file paths in the 'fits' directory
file_paths <- list.files(path = "./sims/k2_asymptotics/fits", pattern = "\\.rds$", full.names = TRUE)

# 2. Load and row-bind them in one step
combined_df <- map(file_paths, readRDS) |> 
  list_rbind()

#Load the combined fits object into the performance function along with the truth object-------------------

performance_res <- performance(fits = combined_df, truth_res = truth_res)

#Save the results in the test_sim folder as "summary.rds"

saveRDS(performance_res, "./sims/k2_asymptotics/summary.rds")
