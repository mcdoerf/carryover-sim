#############################Michael Doerfler###################################
#############################2026-10-13#########################################
#This script contains code for writing the simulation studies in R.


# scenario_seed() --------------------------------------------------------------

##Creates a seed for a simulation scenario.
scenario_seed <- function(base_seed, scenario_id){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed + scenario_id)
  .Random.seed
}

# rep_seeds()-------------------------------------------------------------------
##From the scenario's seed, step through streams so that replication r gets stream r.
rep_seeds <- function(seed, R) {
  seeds <- vector("list", R)
  for (r in seq_len(R)) {
    seed <- parallel::nextRNGStream(seed)
    seeds[[r]] <- seed
  }
  seeds
}
# generate_scenario() ----------------------------------------------------------
#This functions simulates R datasets for one scenario using the specifications in design_row.
#It returns a list of R cr_data objects.

generate_scenario <- function(design_row, R, base_seed) {
  stopifnot("design_row must be a single row" = nrow(design_row) == 1)
  seed  <- scenario_seed(base_seed, design_row$scenario_id)
  seeds <- rep_seeds(seed, R)
  
  sim_one <- switch(design_row$mechanism,
                    multinomial = function() sim_multinomial(design_row$theta1[[1]], design_row$theta2[[1]],
                                                             design_row$N1, design_row$N2,
                                                             design_row$N_tot1, design_row$N_tot2, design_row$psi),
                    stratified  = function() sim_stratified(design_row$theta1_mat[[1]], design_row$theta2_mat[[1]],
                                                            design_row$N1_vec[[1]], design_row$N2_vec[[1]],
                                                            design_row$N_tot1_vec[[1]], design_row$N_tot2_vec[[1]],
                                                            design_row$psi),
                    stop("unknown mechanism: ", design_row$mechanism)
  )
  
  datasets <- vector("list", R)
  for (r in seq_len(R)) {
    assign(".Random.seed", seeds[[r]], envir = globalenv())
    datasets[[r]] <- sim_one()
  }
  
  list(datasets = datasets, design_row = design_row, seed = seed, R = R)
}



# fit_scenario()----------------------------------------------------------------
#This functions takes sim as input, where sim is the output from generate_scenario,
#along with methods, which is a named list of functions that each take a cr_data object.

fit_scenario <- function(sim, methods){
  datasets <- sim$datasets
  row <- sim$design_row
  
  results <- vector(mode = "list", length = length(datasets) * length(methods))

  for (i in 1:length(datasets)){
    for(j in 1:length(methods)){
      tic("Method time")
      
      method_result <- tryCatch(
        cbind(methods[[j]](datasets[[i]], row), error = NA_character_),
        error = function(e) data.frame(target = NA, method = NA,
                                       est = NA, se = NA, lower = NA, upper = NA,
                                       converged = NA, boundary = NA, truncated = NA, 
                                       error = conditionMessage(e))
      )
      
      time_result <- toc(quiet = TRUE)
      
      elapsed_time <- time_result$toc - time_result$tic
      
      method_result$scenario_id <- row$scenario_id
      
      method_result$rep <- i
      
      method_result$time <- elapsed_time
      
      method_result$fit_names <- names(methods)[j]
      
      results[[(i - 1) * length(methods) + j]] <- method_result
      
    }
  }
  
 do.call(rbind, results)
  
}



# run_study() ------------------------------------------------------------------
#This function loops (in parallel) over scenarios for one stage, writing one file per 
#scenario and skipping scenarios whose output already exists.

run_study <- function(study_dir, stage = c("generate", "fit"), scenarios = NULL, R = NULL, base_seed = NULL, methods = NULL, cores = 1){
  
  stage <- match.arg(stage)
  
  # ---- Setup shared by both stages -------------------------------------------
  design_grid <- readRDS(paste0(study_dir, "/design.rds"))
  
  if (is.null(scenarios)) scenarios <- design_grid$scenario_id
  
  stopifnot("scenarios must lie in scenario_id" = all(scenarios %in% design_grid$scenario_id))

  
  if (stage == "generate"){
    stopifnot("R must be an positive integer" = is.numeric(R) && length(R) == 1 &&  R >= 1 && round(R) == R)
    stopifnot("Base seed must be specified" = is.numeric(base_seed))
  }
  
  if (stage == "fit"){
    stopifnot("Methods must be specified" = !is.null(methods))
  }
  
  data_dir <- paste0(study_dir, "/data")
  fits_dir <- paste0(study_dir, "/fits")
  
  # --- One scenario, generate stage

  generate_one <- function(id){
    
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
    
    out_file <- paste0(data_dir, sprintf("/scenario_%03d.rds", id))
    
    if (file.exists(out_file)) { message("scenario ", id, ": exists, skipped"); return(invisible("skipped"))}
    row <- design_grid[design_grid$scenario_id == id, ]
    
    sim <- generate_scenario(row, R = R, base_seed = base_seed)
    
    saveRDS(sim, file = out_file)
    
    message("scenario ", id, ": generated")
    
    invisible("done")
  }
  
  

  # --- One scenario, fit stage
  
  fit_one <- function(id){
    
    dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
    
    in_file <- paste0(data_dir, sprintf("/scenario_%03d.rds", id))
    out_file <- paste0(fits_dir, sprintf("/scenario_%03d.rds", id))
    
    if (file.exists(out_file)) { message("scenario ", id, ": exists, skipped"); return(invisible("skipped"))}
    
    if (!file.exists(in_file)) {stop("generate scenario ", id, "first")}
    
    sim <- readRDS(in_file)
    
    fits <- fit_scenario(sim = sim, methods = methods)
    
    saveRDS(fits, out_file)
    
    message("scenario ", id, ": fitted")
    
    invisible("done")
  }
  
  # --- Run every requested scenario
  
  do_one <- if (stage == "generate") generate_one else fit_one
  
  if (cores == 1){
    status <- lapply(scenarios, do_one)
  } else{
    future::plan(future::multisession, workers = cores)
    status <- future.apply::future_lapply(scenarios, do_one, future.seed = TRUE)
  }
  
  # --- Small summary of what happened
  invisible(data.frame(scenario_id = scenarios, status = unlist(status)))
  
}