#############################Michael Doerfler###################################
#############################2026-10-13#########################################
#This script contains R code for putting together performance measures for a simulation study.

# truth() ----------------------------------------------------------------------
#This function takes the design grid for a simulation study and builds a small d.f.
#of the true values, one row per scenario per target.

truth_generator <- function(design){
  design_long <- design %>%
    select(scenario_id, N1, N2) %>%
    pivot_longer(
      cols = c(N1, N2),
      names_to = "target",
      values_to = "truth"
    )
  
  return(design_long)
}

# performance() ----------------------------------------------------------------
#This function takes as input all the raw fit results for a study, stacked into on df
# (named fits).  It also takes truth_res, the output of the truth_generator() function.
#Returns: d.f. with one row per scenario \times fit_names \times method \times target, with
#performance summary measures.

performance <- function(fits, truth_res){
  join_df <- left_join(fits, truth_res, by = c("scenario_id", "target"))
  
  failures <- fits %>%
    distinct(scenario_id, rep, fit_names, error) %>%
    group_by(scenario_id, fit_names) %>%
    summarise(n_rep = n(), n_fail = sum(!is.na(error)), .groups = "drop")
  
  perf_summary <- join_df %>%
    filter(is.na(error)) %>%
    group_by(scenario_id, fit_names, method, target) %>%
    summarise(
      n_ok = n(),
      bias = mean(est - truth, na.rm = TRUE),
      bias_mcse = sd(est) / sqrt(n()),
      rel_bias = mean(est - truth, na.rm = TRUE) / mean(truth, na.rm = TRUE),
      emp_sd = sd(est, na.rm = TRUE),
      emp_sd_mcse = emp_sd / sqrt(2 * (n() - 1)),
      mean_se = sqrt(mean(se^2, na.rm = TRUE)),
      rmse = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      coverage = mean(lower <= truth & upper >= truth, na.rm = TRUE),
      coverage_mcse = sqrt(coverage * (1 - coverage) / n()),
      width = mean(upper - lower, na.rm = TRUE),
      width_mcse = sd(upper - lower) / sqrt(n()),
      rel_width = mean(upper - lower, na.rm = TRUE) / mean(truth, na.rm = TRUE),
      .groups = "drop"
      )
  
  left_join(perf_summary, failures, by = c("scenario_id", "fit_names")) %>%
    mutate(fail_rate = n_fail / n_rep)
}