#############################Michael Doerfler###################################
#############################2026-10-20#########################################
#This script contains R code for putting together performance measures for the test_sim_full_MLE
#test simulation study.

library(tidyverse)
library(here)
library(tictoc)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}


#Load the design grid corresponding to this simulation study.-------------------------------
design <- readRDS("./sims/test_sim_full_MLE_2/design.rds")

#Create the truth d.f. ---------------------------------------------------------------------
truth_res <- truth_generator(design = design)

#Extract the fits objects and row bind them into one large object.--------------------------

# 1. Get a vector of all .rds file paths in the 'fits' directory
file_paths <- list.files(path = "./sims/test_sim_full_MLE_2/fits", pattern = "\\.rds$", full.names = TRUE)

# 2. Load and row-bind them in one step
combined_df <- map(file_paths, readRDS) |> 
  list_rbind()

#Load the combined fits object into the performance function along with the truth object-------------------

performance_res <- performance(fits = combined_df, truth_res = truth_res)

#Save the results in the test_sim folder as "summary.rds"

saveRDS(performance_res, "./sims/test_sim_full_MLE_2/summary.rds")
