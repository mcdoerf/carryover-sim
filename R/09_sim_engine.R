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



# run_study --------------------------------------------------------------------
# Runs one stage ("generate" or "fit") of a simulation study over a set of
# scenarios, writing one file per scenario and skipping scenarios whose output
# file already exists (so interrupted runs can simply be restarted).
#
#   study_dir : folder containing design.rds; data/ and fits/ are created inside it
#   stage     : "generate" (needs R, base_seed) or "fit" (needs methods)
#   scenarios : scenario_ids to run (default: all in the design)
#   cores     : 1 = sequential; > 1 = parallel via future (works on Windows)
#
# Returns (invisibly) a data frame: scenario_id, status, seconds.

run_study <- function(study_dir, stage = c("generate", "fit"), scenarios = NULL,
                      R = NULL, base_seed = NULL, methods = NULL, cores = 1) {
  
  stage <- match.arg(stage)
  
  # ---- Setup shared by both stages ------------------------------------------
  design_grid <- readRDS(file.path(study_dir, "design.rds"))
  
  if (is.null(scenarios)) scenarios <- design_grid$scenario_id
  stopifnot(
    "scenarios must all be scenario_ids in the design" =
      all(scenarios %in% design_grid$scenario_id),
    "cores must be a positive integer" =
      is.numeric(cores) && length(cores) == 1 && cores >= 1 && cores == round(cores)
  )
  
  if (stage == "generate") {
    stopifnot(
      "R must be a positive integer" =
        is.numeric(R) && length(R) == 1 && R >= 1 && R == round(R),
      "base_seed must be a single number" =
        is.numeric(base_seed) && length(base_seed) == 1
    )
  }
  if (stage == "fit") {
    stopifnot("methods must be a named list of functions" =
                is.list(methods) && !is.null(names(methods)) &&
                all(vapply(methods, is.function, logical(1))))
  }
  
  data_dir <- file.path(study_dir, "data")
  fits_dir <- file.path(study_dir, "fits")
  dir.create(if (stage == "generate") data_dir else fits_dir,
             showWarnings = FALSE, recursive = TRUE)
  
  scenario_file <- function(dir, id) file.path(dir, sprintf("scenario_%03d.rds", id))
  
  # Save to a temporary name, then rename: a file with the final name is
  # always complete, even if R crashes mid-save.
  safe_save <- function(object, file) {
    tmp <- paste0(file, ".tmp")
    saveRDS(object, tmp)
    if (!file.rename(tmp, file)) stop("could not rename ", tmp, " to ", file)
  }
  
  # ---- One scenario, generate stage -----------------------------------------
  generate_one <- function(id) {
    out_file <- scenario_file(data_dir, id)
    if (file.exists(out_file)) return("skipped")
    
    row <- design_grid[design_grid$scenario_id == id, ]
    sim <- generate_scenario(row, R = R, base_seed = base_seed)
    sim$created   <- Sys.time()
    sim$r_version <- R.version.string
    safe_save(sim, out_file)
    "done"
  }
  
  # ---- One scenario, fit stage ----------------------------------------------
  fit_one <- function(id) {
    in_file  <- scenario_file(data_dir, id)
    out_file <- scenario_file(fits_dir, id)
    if (file.exists(out_file)) return("skipped")
    if (!file.exists(in_file)) stop("no data for scenario ", id, "; run the generate stage first")
    
    sim  <- readRDS(in_file)
    fits <- fit_scenario(sim, methods)
    attr(fits, "created")   <- Sys.time()
    attr(fits, "r_version") <- R.version.string
    safe_save(fits, out_file)
    "done"
  }
  
  # ---- Wrapper: timing, error capture, progress message ---------------------
  # An error in one scenario is recorded in the status table instead of
  # stopping the remaining scenarios.
  do_stage <- if (stage == "generate") generate_one else fit_one
  
  run_one <- function(id) {
    t0 <- proc.time()[["elapsed"]]
    status <- tryCatch(do_stage(id), error = function(e) paste("error:", conditionMessage(e)))
    secs <- round(proc.time()[["elapsed"]] - t0, 1)
    message("scenario ", id, ": ", status, " (", secs, " s)")
    data.frame(scenario_id = id, status = status, seconds = secs)
  }
  
  # ---- Run every requested scenario -----------------------------------------
  if (cores == 1) {
    results <- lapply(scenarios, run_one)
  } else {
    # Workers are fresh R sessions (always on Windows), so each loads the
    # project's functions once, into its global environment, before its first
    # scenario. A flag is used rather than exists("some_function") because
    # future copies some of your functions to the workers itself, which would
    # fool that check while others (e.g. those only called inside `methods`)
    # are still missing.
    r_files <- normalizePath(sort(list.files(here::here("R"), pattern = "\\.R$",
                                             full.names = TRUE)))
    old_plan <- future::plan(future::multisession, workers = cores)
    on.exit(future::plan(old_plan), add = TRUE)
    
    results <- future.apply::future_lapply(scenarios, function(id) {
      if (!exists(".project_loaded", envir = globalenv(), inherits = FALSE)) {
        for (f in r_files) source(f, local = globalenv())
        assign(".project_loaded", TRUE, envir = globalenv())
      }
      run_one(id)
    }, future.seed = TRUE)
  }
  
  # ---- Summary ----------------------------------------------------------------
  summary <- do.call(rbind, results)
  n_err <- sum(startsWith(summary$status, "error"))
  if (n_err > 0) warning(n_err, " scenario(s) failed; see the returned status table")
  invisible(summary)
}
