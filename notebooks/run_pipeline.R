# ==============================================================================
# run_pipeline.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Single entry point for the entire simulation study. Sources all module
#   files in dependency order, runs the calibration and Monte Carlo grid,
#   and generates all manuscript tables and figures. This is the ONLY file
#   that should be executed directly by the user.
#
# USAGE
#   Rscript run_pipeline.R                  # sequential, default settings
#   Rscript run_pipeline.R --cores 4        # parallel, 4 cores (non-Windows)
#
#   Or interactively in R:
#   source("run_pipeline.R")
#
# PIPELINE EXECUTION ORDER
#   run_pipeline.R
#     -> 01_generate_data.R    data-generating mechanism + marginal truth
#     -> 02_apply_dropout.R    calibrated monotone MAR dropout
#     -> 03_fit_models.R       LME, GEE (Gaussian), GLMM, GEE (Logistic)
#     -> 04_metrics.R          performance measures + MCSEs
#     -> 05_simulation_engine.R  Monte Carlo orchestration + checkpointing
#     -> 06_tables_figures.R   publication tables (LaTeX/CSV) + figures (PDF)
#
# REPRODUCIBILITY
#   The master seed (MASTER_SEED = 20240101) in 05_simulation_engine.R
#   uniquely determines all Monte Carlo results. Changing any package
#   version may alter numerical results; the session information is therefore
#   saved alongside the results as sessioninfo.txt.
#
# ESTIMATED RUNTIME
#   N = 1000 replicates x 24 cells x 4 models (sequential, modern laptop):
#   approximately 2-4 hours. With n_cores = 4: approximately 40-60 minutes.
#   A dry run (n_replicates = 2) completes in under 2 minutes.
#
# DEPENDENCIES (must be installed before running)
#   lme4, lmerTest, geepack, data.table, ggplot2, xtable
#   Install with:
#   install.packages(c("lme4","lmerTest","geepack","data.table","ggplot2","xtable"))
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY CHECK
# ------------------------------------------------------------------------------
required_packages <- c("lme4", "lmerTest", "geepack", "data.table", "ggplot2", "xtable")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are not installed:\n  ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with:\n  install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. COMMAND-LINE ARGUMENT PARSING
# ------------------------------------------------------------------------------
# Supports two optional arguments when run via Rscript:
#   --cores N       number of parallel cores (default: auto-detected)
#   --reps  N       number of Monte Carlo replicates (default: 1000)
#   --out   PATH    output directory (default: "results")
#   --dryrun        shorthand for --reps 2 (quick end-to-end verification)

.parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Defaults
  cfg <- list(
    n_cores     = NULL,    # NULL -> auto-detect below
    n_reps      = 1000L,
    output_dir  = "results",
    dry_run     = FALSE
  )
  
  i <- 1L
  while (i <= length(args)) {
    switch(args[i],
           "--cores"   = { cfg$n_cores    <- as.integer(args[i + 1L]); i <- i + 2L },
           "--reps"    = { cfg$n_reps     <- as.integer(args[i + 1L]); i <- i + 2L },
           "--out"     = { cfg$output_dir <- args[i + 1L];             i <- i + 2L },
           "--dryrun"  = { cfg$dry_run    <- TRUE;                     i <- i + 1L },
           { warning(sprintf("Unknown argument ignored: %s", args[i]), call. = FALSE)
             i <- i + 1L }
    )
  }
  
  if (cfg$dry_run) {
    cfg$n_reps <- 2L
    cfg$output_dir <- paste0(cfg$output_dir, "_dryrun")
    message("DRY RUN mode: n_replicates = 2. Output -> ", cfg$output_dir)
  }
  
  cfg
}

cfg <- .parse_args()

# ------------------------------------------------------------------------------
# 3. CORE DETECTION AND PLATFORM-AWARE PARALLELISM
# ------------------------------------------------------------------------------
# parallel::mclapply (used in 05_simulation_engine.R) relies on OS-level
# process forking and is NOT supported on Windows. On Windows, n_cores is
# forced to 1 regardless of what was requested, with an explicit warning.

.detect_cores <- function(requested) {
  on_windows <- .Platform$OS.type == "windows"
  
  if (on_windows) {
    if (!is.null(requested) && requested > 1L) {
      warning(
        "Parallel execution via mclapply is not supported on Windows.\n",
        "Falling back to sequential execution (n_cores = 1).\n",
        "For Windows parallelism, consider running this pipeline under WSL2.",
        call. = FALSE
      )
    }
    return(1L)
  }
  
  # Non-Windows: use requested value or leave one core free for the OS
  if (!is.null(requested)) {
    return(max(1L, as.integer(requested)))
  }
  
  available <- parallel::detectCores(logical = FALSE)
  # Default: leave one physical core free to keep the system responsive
  max(1L, available - 1L)
}

n_cores <- .detect_cores(cfg$n_cores)

# ------------------------------------------------------------------------------
# 4. SOURCE ALL MODULES IN DEPENDENCY ORDER
# ------------------------------------------------------------------------------
# Each file is sourced exactly once. The order below is the only valid order:
# each file may call functions defined in earlier files.

module_files <- c(
  "01_generate_data.R",
  "02_apply_dropout.R",
  "03_fit_models.R",
  "04_metrics.R",
  "05_simulation_engine.R",
  "06_tables_figures.R"
)

missing_modules <- module_files[!file.exists(module_files)]
if (length(missing_modules) > 0) {
  stop(
    "The following module files are missing from the working directory:\n  ",
    paste(missing_modules, collapse = "\n  "),
    "\nEnsure run_pipeline.R is executed from the project root directory.",
    call. = FALSE
  )
}

message("=================================================================")
message("Simulation study: LME vs GEE vs GLMM under Monotone Dropout")
message("=================================================================")
message(sprintf("Platform        : %s (%s)", R.version$platform, .Platform$OS.type))
message(sprintf("R version       : %s", R.version.string))
message(sprintf("Cores requested : %d", n_cores))
message(sprintf("Replicates/cell : %d", cfg$n_reps))
message(sprintf("Output directory: %s", cfg$output_dir))
message(sprintf("Design cells    : %d (2 scenarios x 3 N x 4 dropout rates)", 2*3*4))
message("=================================================================")

for (f in module_files) {
  message(sprintf("Sourcing %s ...", f))
  source(f, local = FALSE)
}
message("All modules loaded.")

# ------------------------------------------------------------------------------
# 5. RUN SIMULATION GRID
# ------------------------------------------------------------------------------

message("\nStarting simulation grid ...")

metrics_table <- run_simulation_grid(
  n_replicates = cfg$n_reps,
  n_calib      = if (cfg$dry_run) 2000L else 20000L,
  output_dir   = cfg$output_dir,
  log_path     = NULL,
  n_cores      = n_cores
)

message(sprintf(
  "\nSimulation complete. Metrics table: %d rows x %d columns.",
  nrow(metrics_table), ncol(metrics_table)
))

# ------------------------------------------------------------------------------
# 6. GENERATE TABLES AND FIGURES
# ------------------------------------------------------------------------------

tables_figures_dir <- file.path(cfg$output_dir, "tables_figures")

message(sprintf("\nGenerating tables and figures -> %s", tables_figures_dir))

generate_all_outputs(
  metrics_path = file.path(cfg$output_dir, "metrics_table.rds"),
  output_dir   = tables_figures_dir
)

# ------------------------------------------------------------------------------
# 7. SAVE SESSION INFORMATION
# ------------------------------------------------------------------------------
# Session info is essential for a reproducible, publicly released codebase.
# It records exact package versions so that results can be reproduced on
# another machine or at a future date.

session_path <- file.path(cfg$output_dir, "sessioninfo.txt")
writeLines(capture.output(sessionInfo()), con = session_path)
message(sprintf("Session info saved to: %s", session_path))

# ------------------------------------------------------------------------------
# 8. FINAL SUMMARY
# ------------------------------------------------------------------------------

message("\n=================================================================")
message("PIPELINE COMPLETE")
message("=================================================================")
message(sprintf("Raw results   : %s/all_cell_results.rds", cfg$output_dir))
message(sprintf("Metrics table : %s/metrics_table.rds",   cfg$output_dir))
message(sprintf("Tables (LaTeX): %s/*.tex",                tables_figures_dir))
message(sprintf("Figures (PDF) : %s/*.pdf",                tables_figures_dir))
message(sprintf("Log file      : %s/simulation.log",       cfg$output_dir))
message(sprintf("Session info  : %s/sessioninfo.txt",      cfg$output_dir))
message("=================================================================")