# ==============================================================================
# 05_simulation_engine.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Orchestrates the full Monte Carlo simulation study by iterating over the
#   complete (scenario x n x target_dropout_rate) design grid, running B=1000
#   replicates per cell, and collecting raw replicate-level results for
#   downstream metric computation in 04_metrics.R. This file is the single
#   point of entry for the computationally intensive part of the study.
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#
#   1. CENTRALIZED, REPRODUCIBLE PER-REPLICATE SEEDING
#      The simulation uses L'Ecuyer-CMRG as the RNG kind, which is designed
#      specifically for parallel/independent streams and is the recommended
#      choice for Monte Carlo studies that may run in parallel (R Core Team;
#      L'Ecuyer 1999). Each replicate is assigned a unique, deterministic seed
#      derived from a single master seed combined with cell- and replicate-
#      level indices. This guarantees:
#        a) Full reproducibility: the same results are obtained on any
#           machine, in any cell order, regardless of whether cells are
#           processed sequentially or in parallel.
#        b) Independence between cells: the RNG state for replicate b in cell
#           c does not depend on the execution of any other cell.
#        c) Independence between replicates within a cell: each replicate
#           draws from a distinct, non-overlapping RNG stream.
#      Using set.seed(b) naively inside a loop -- a common but incorrect
#      shortcut -- would couple seeds across cells and make parallel
#      execution non-reproducible.
#
#   2. DROPOUT CALIBRATION IS PERFORMED ONCE, BEFORE THE SIMULATION LOOP
#      As established in 02_apply_dropout.R, the calibrated hazard intercepts
#      are properties of the data-generating process, not of any particular
#      Monte Carlo draw. Computing them once (at high pseudo-population size)
#      and reusing them across all N-levels and all 1000 replicates avoids
#      redundant calibration work and guarantees that every cell in the grid
#      targets the same population-level dropout probability.
#
#   3. SIMULATION GRID USES CJ() (CROSS JOIN) FROM data.table
#      The full factorial design is expressed as a single data.table cross-
#      join of all factor levels. This makes the study design self-documenting
#      (the grid object IS the protocol table) and avoids nested for-loops
#      whose depth would obscure the cell structure.
#
#   4. OPTIONAL PARALLEL EXECUTION VIA parallel::mclapply
#      Cells are independent and embarrassingly parallel. When n_cores > 1,
#      cells are distributed across cores using parallel::mclapply with
#      mc.set.seed = FALSE (because seeding is already handled per-replicate
#      via the centralized scheme described in point 1 above). On Windows,
#      mclapply falls back to lapply with a warning (Windows does not support
#      forking). This is the simplest correct parallelization strategy that
#      does not require additional dependencies (e.g. foreach, doParallel).
#
#   5. PROGRESS LOGGING TO CONSOLE AND FILE
#      Each cell's start/end time and achieved dropout rate (mean across
#      replicates) are logged to both the console and an optional log file.
#      This is essential for a simulation study of this size: if a run is
#      interrupted, the log identifies exactly how far the grid was completed,
#      and which cells need to be re-run. Elapsed time per cell is also
#      reported to allow realistic estimation of total compute time.
#
#   6. DEFENSIVE ERROR ISOLATION AT THE REPLICATE LEVEL
#      Each replicate is wrapped in tryCatch so that a single failed replicate
#      (e.g. due to a numerical edge case not caught by 03_fit_models.R's own
#      error handling) does not abort the entire cell. Failed replicates are
#      replaced with a structured NA placeholder and counted as part of the
#      failure rate reported in the log. This is important at N=100 where
#      edge-case model failures are most likely.
#
#   7. RESULTS ARE SAVED INCREMENTALLY (CELL BY CELL)
#      After each cell completes, its raw results are written to an RDS file
#      in a user-specified output directory. This provides natural checkpointing:
#      if the simulation is interrupted mid-grid, completed cells do not need
#      to be re-run. The final step assembles the full metrics table from the
#      saved RDS files and writes a single consolidated RDS for 06_tables_figures.R.
#
# DEPENDENCIES
#   data.table, parallel (base R)
#
# UPSTREAM DEPENDENCIES (assumed already sourced by run_pipeline.R)
#   01_generate_data.R   -- generate_replicate_data()
#   02_apply_dropout.R   -- build_dropout_calibration_table(),
#                           apply_monotone_dropout(),
#                           get_dropout_hazard_slopes()
#   03_fit_models.R      -- fit_all_models()
#   04_metrics.R         -- compute_cell_metrics(), compute_all_metrics()
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# ------------------------------------------------------------------------------
# 1. STUDY DESIGN CONSTANTS
# ------------------------------------------------------------------------------

#' Full simulation grid: all (scenario x n x target_dropout_rate) cells.
#' This object IS the protocol table -- every combination run by the engine.
SIM_GRID <- CJ(
  scenario            = c("null", "alternative"),
  n                   = c(100L, 200L, 500L),
  target_dropout_rate = c(0.10, 0.20, 0.30, 0.40),
  sorted              = FALSE
)

#' Number of Monte Carlo replicates per cell.
N_REPLICATES <- 1000L

#' Master RNG seed. Changing this value and re-running the full pipeline
#' will produce a statistically independent but equally valid set of results.
MASTER_SEED <- 20240101L

# ------------------------------------------------------------------------------
# 2. PER-REPLICATE SEEDING UTILITIES
# ------------------------------------------------------------------------------

#' Compute a unique, deterministic integer seed for one (cell, replicate) pair.
#'
#' @param cell_index integer; row index of the cell in SIM_GRID.
#' @param rep_index integer; replicate index (1 to N_REPLICATES).
#' @param master_seed integer; the study-level master seed.
#' @return integer scalar; unique seed for this (cell, replicate) combination.
#' @details Uses a simple Cantor-pairing-inspired formula to map the 2D
#'   (cell, rep) index to a unique integer in the range required by set.seed().
#'   The formula is exact (no collisions) for cell_index <= 1e4 and
#'   rep_index <= 1e4, both of which are well above the study's requirements.
.make_replicate_seed <- function(cell_index, rep_index, master_seed = MASTER_SEED) {
  # Cantor pairing: N(a,b) = (a+b)*(a+b+1)/2 + b, then offset by master_seed
  a <- as.integer(cell_index - 1L)
  b <- as.integer(rep_index  - 1L)
  pair <- ((a + b) * (a + b + 1L)) %/% 2L + b
  (master_seed + pair) %% .Machine$integer.max
}

# ------------------------------------------------------------------------------
# 3. LOGGING UTILITIES
# ------------------------------------------------------------------------------

#' Write a timestamped message to the console and optionally to a log file.
#'
#' @param msg character; message text.
#' @param log_path character or NULL; if non-NULL, message is appended to file.
.log <- function(msg, log_path = NULL) {
  stamped <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(stamped)
  if (!is.null(log_path)) {
    cat(stamped, "\n", file = log_path, append = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 4. SINGLE-REPLICATE RUNNER (used inside the cell loop)
# ------------------------------------------------------------------------------

#' Run one Monte Carlo replicate for a given cell.
#'
#' @param cell_index integer; row index in SIM_GRID (for seed derivation).
#' @param rep_index integer; replicate number within the cell.
#' @param n integer; number of subjects.
#' @param scenario character; "null" or "alternative".
#' @param target_dropout_rate numeric; target dropout rate.
#' @param calibration_table data.table; from build_dropout_calibration_table().
#' @param slopes list; from get_dropout_hazard_slopes().
#' @return the output of fit_all_models() on this replicate's data, or a
#'   structured NA placeholder if an unrecoverable error occurred.
.run_one_replicate <- function(cell_index,
                               rep_index,
                               n,
                               scenario,
                               target_dropout_rate,
                               calibration_table,
                               slopes) {
  seed <- .make_replicate_seed(cell_index, rep_index)
  
  tryCatch({
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    
    rep_data   <- generate_replicate_data(n = n, scenario = scenario)
    do_data    <- apply_monotone_dropout(rep_data, target_dropout_rate,
                                         calibration_table, slopes)
    fit_result <- fit_all_models(do_data)
    
    # Attach the achieved dropout rate and replicate index for diagnostics
    fit_result$achieved_dropout_rate <- do_data$achieved_dropout_rate
    fit_result$rep_index             <- rep_index
    fit_result
    
  }, error = function(e) {
    # Structured placeholder: fit_all_models()-compatible stub with failed=TRUE
    methods <- c("LME", "GEE_continuous", "GLMM", "GEE_binary")
    failed_fits <- lapply(methods, function(m) {
      list(method = m, estimate = NA_real_, se_robust = NA_real_,
           se_naive = NA_real_, pvalue = NA_real_, converged = FALSE,
           singular = NA, failed = TRUE,
           failure_reason = conditionMessage(e))
    })
    names(failed_fits) <- methods
    
    list(
      fits                 = failed_fits,
      truth                = list(beta3 = NA_real_,
                                  alpha3_conditional = NA_real_,
                                  alpha3_marginal = NA_real_),
      params               = list(scenario = scenario),
      achieved_dropout_rate = NA_real_,
      rep_index            = rep_index
    )
  })
}

# ------------------------------------------------------------------------------
# 5. SINGLE-CELL RUNNER
# ------------------------------------------------------------------------------

#' Run all N_REPLICATES replicates for one simulation cell.
#'
#' @param cell_index integer; row index in SIM_GRID.
#' @param n integer; number of subjects.
#' @param scenario character; "null" or "alternative".
#' @param target_dropout_rate numeric; target dropout rate.
#' @param calibration_table data.table; from build_dropout_calibration_table().
#' @param slopes list; from get_dropout_hazard_slopes().
#' @param n_replicates integer; number of Monte Carlo replicates.
#' @param output_dir character; directory to save per-cell RDS file.
#' @param log_path character or NULL; path for logging.
#' @return a list with:
#'   replicate_results  - list of length n_replicates (raw fit outputs),
#'   n, scenario, target_dropout_rate (cell metadata, echoed for 04_metrics.R),
#'   mean_achieved_dropout_rate - mean realized dropout rate across replicates,
#'   n_failed - count of replicates that hit the tryCatch error handler.
.run_one_cell <- function(cell_index,
                          n,
                          scenario,
                          target_dropout_rate,
                          calibration_table,
                          slopes,
                          n_replicates = N_REPLICATES,
                          output_dir   = NULL,
                          log_path     = NULL) {
  cell_label <- sprintf("scenario=%s | n=%d | dropout=%.0f%%",
                        scenario, n, target_dropout_rate * 100)
  .log(sprintf("START cell %02d/%02d: %s",
               cell_index, nrow(SIM_GRID), cell_label), log_path)
  t_start <- proc.time()["elapsed"]
  
  reps <- lapply(seq_len(n_replicates), function(b) {
    .run_one_replicate(
      cell_index          = cell_index,
      rep_index           = b,
      n                   = n,
      scenario            = scenario,
      target_dropout_rate = target_dropout_rate,
      calibration_table   = calibration_table,
      slopes              = slopes
    )
  })
  
  # Diagnostics
  n_failed <- sum(vapply(reps, function(r) {
    any(vapply(r$fits, function(f) isTRUE(f$failed), logical(1)))
  }, logical(1)))
  
  achieved_rates <- vapply(reps, function(r) {
    if (is.null(r$achieved_dropout_rate) || is.na(r$achieved_dropout_rate))
      NA_real_
    else
      r$achieved_dropout_rate
  }, numeric(1))
  mean_achieved <- mean(achieved_rates, na.rm = TRUE)
  
  elapsed <- proc.time()["elapsed"] - t_start
  .log(
    sprintf(
      "END   cell %02d/%02d: %s | elapsed=%.1fs | mean_achieved_dropout=%.3f | n_failed_reps=%d",
      cell_index, nrow(SIM_GRID), cell_label, elapsed, mean_achieved, n_failed
    ),
    log_path
  )
  
  cell_result <- list(
    replicate_results          = reps,
    n                          = n,
    scenario                   = scenario,
    target_dropout_rate        = target_dropout_rate,
    mean_achieved_dropout_rate = mean_achieved,
    n_failed                   = n_failed
  )
  
  # Incremental save: one RDS per cell for checkpointing
  if (!is.null(output_dir)) {
    rds_path <- file.path(
      output_dir,
      sprintf("cell_%02d_scenario%s_n%d_dropout%02d.rds",
              cell_index, scenario, n,
              as.integer(target_dropout_rate * 100))
    )
    saveRDS(cell_result, file = rds_path)
  }
  
  cell_result
}

# ------------------------------------------------------------------------------
# 6. FULL SIMULATION GRID RUNNER (TOP-LEVEL ENTRY POINT)
# ------------------------------------------------------------------------------

#' Run the complete simulation study across all cells in SIM_GRID.
#'
#' @param n_replicates integer; Monte Carlo replicates per cell (default 1000).
#' @param n_calib integer; pseudo-population size for dropout calibration.
#' @param output_dir character; directory for incremental RDS saves and the
#'   final consolidated metrics table. Created if it does not exist.
#' @param log_path character or NULL; path for a plain-text log file.
#' @param n_cores integer; number of parallel cores (cells are distributed
#'   across cores via parallel::mclapply). Use 1 for sequential execution.
#'   On Windows, values > 1 fall back to 1 with a warning.
#' @return a data.table of performance metrics (one row per method per cell),
#'   as produced by compute_all_metrics(). Also saved as "metrics_table.rds"
#'   in output_dir.
#' @details Run this function from run_pipeline.R. The typical call is:
#'   results <- run_simulation_grid(output_dir = "results/")
run_simulation_grid <- function(n_replicates = N_REPLICATES,
                                n_calib      = 20000L,
                                output_dir   = "results",
                                log_path     = file.path(output_dir, "simulation.log"),
                                n_cores      = 1L) {
  # ---- Setup ----------------------------------------------------------------
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (!is.null(log_path)) {
    cat("", file = log_path)   # initialize / clear log file
  }
  
  .log("============================================================", log_path)
  .log("SIMULATION STUDY: LME vs GEE vs GLMM under Monotone Dropout", log_path)
  .log(sprintf("Grid: %d cells | %d replicates | master seed: %d",
               nrow(SIM_GRID), n_replicates, MASTER_SEED), log_path)
  .log("============================================================", log_path)
  
  # Validate n_cores
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    warning("parallel::mclapply is not supported on Windows. Falling back to n_cores=1.",
            call. = FALSE)
    n_cores <- 1L
  }
  
  # ---- Dropout calibration (performed ONCE before the simulation loop) -----
  .log("Starting dropout calibration ...", log_path)
  slopes <- get_dropout_hazard_slopes()
  calibration_table <- build_dropout_calibration_table(
    n_calib = n_calib,
    seed    = MASTER_SEED
  )
  .log("Dropout calibration complete. Calibrated intercepts:", log_path)
  .log(capture.output(print(calibration_table)) |> paste(collapse = "\n"), log_path)
  
  # ---- Cell runner wrapper (index-based, for mclapply compatibility) --------
  run_cell_by_index <- function(cell_index) {
    row <- SIM_GRID[cell_index]
    .run_one_cell(
      cell_index          = cell_index,
      n                   = row$n,
      scenario            = row$scenario,
      target_dropout_rate = row$target_dropout_rate,
      calibration_table   = calibration_table,
      slopes              = slopes,
      n_replicates        = n_replicates,
      output_dir          = output_dir,
      log_path            = log_path
    )
  }
  
  # ---- Execute across grid --------------------------------------------------
  t_total_start <- proc.time()["elapsed"]
  
  all_cell_results <- if (n_cores > 1L) {
    parallel::mclapply(
      seq_len(nrow(SIM_GRID)),
      run_cell_by_index,
      mc.cores    = n_cores,
      mc.set.seed = FALSE   # seeding is already handled per-replicate
    )
  } else {
    lapply(seq_len(nrow(SIM_GRID)), run_cell_by_index)
  }
  
  total_elapsed <- proc.time()["elapsed"] - t_total_start
  .log(sprintf("All cells complete. Total elapsed: %.1f minutes.",
               total_elapsed / 60), log_path)
  
  # ---- Compute metrics from raw replicate results ---------------------------
  .log("Computing performance metrics ...", log_path)
  metrics_table <- compute_all_metrics(all_cell_results)
  
  # ---- Save consolidated outputs -------------------------------------------
  metrics_path <- file.path(output_dir, "metrics_table.rds")
  saveRDS(metrics_table, file = metrics_path)
  .log(sprintf("Metrics table saved to: %s", metrics_path), log_path)
  
  raw_path <- file.path(output_dir, "all_cell_results.rds")
  saveRDS(all_cell_results, file = raw_path)
  .log(sprintf("Raw replicate results saved to: %s", raw_path), log_path)
  
  .log("Simulation complete.", log_path)
  metrics_table[]
}

# ------------------------------------------------------------------------------
# 7. CHECKPOINT RECOVERY UTILITY
# ------------------------------------------------------------------------------

#' Reassemble all_cell_results from incrementally saved per-cell RDS files
#' in the event the simulation was interrupted before completion.
#'
#' @param output_dir character; directory containing the per-cell RDS files
#'   written by .run_one_cell().
#' @return list of cell results in SIM_GRID row order (cells that have not
#'   yet been saved are represented as NULL and will trigger a warning).
#' @details After recovery, call compute_all_metrics(recovered) to re-derive
#'   the metrics table, then re-run only the missing cells if needed.
recover_from_checkpoint <- function(output_dir = "results") {
  n_cells <- nrow(SIM_GRID)
  recovered <- vector("list", n_cells)
  
  for (i in seq_len(n_cells)) {
    row <- SIM_GRID[i]
    rds_path <- file.path(
      output_dir,
      sprintf("cell_%02d_scenario%s_n%d_dropout%02d.rds",
              i, row$scenario, row$n,
              as.integer(row$target_dropout_rate * 100))
    )
    if (file.exists(rds_path)) {
      recovered[[i]] <- readRDS(rds_path)
    } else {
      warning(sprintf("Cell %d (%s | n=%d | dropout=%.0f%%) not found -- may still be running.",
                      i, row$scenario, row$n, row$target_dropout_rate * 100),
              call. = FALSE)
    }
  }
  
  recovered
}

# ------------------------------------------------------------------------------
# 8. SMOKE TEST (manual, interactive use only -- not sourced by the pipeline)
# ------------------------------------------------------------------------------
# Run a minimal dry run (2 replicates, first cell only) to verify end-to-end
# flow before committing to a full 1000-replicate run.
#
# source("01_generate_data.R")
# source("02_apply_dropout.R")
# source("03_fit_models.R")
# source("04_metrics.R")
# source("05_simulation_engine.R")
#
# dry_run <- run_simulation_grid(
#   n_replicates = 2,
#   n_calib      = 2000,
#   output_dir   = "results_dryrun",
#   n_cores      = 1L
# )
# print(dry_run[, .(scenario, n, target_dropout_rate, method, bias, coverage)])