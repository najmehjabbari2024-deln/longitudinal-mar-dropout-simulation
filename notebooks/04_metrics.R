# ==============================================================================
# 04_metrics.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Computes all pre-specified performance measures for one simulation design
#   cell (a fixed combination of n, target dropout rate, and scenario) from
#   a list of B Monte Carlo replicate results produced by 03_fit_models.R.
#   Returns a tidy data.table of metrics, one row per method, ready for
#   aggregation and display in 06_tables_figures.R.
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#
#   1. ESTIMAND-AWARE TRUTH ASSIGNMENT
#      Each model is evaluated against its OWN correct estimand:
#        LME          -> truth$beta3              (conditional = marginal for
#        GEE_cont.    -> truth$beta3               identity link)
#        GLMM         -> truth$alpha3_conditional  (subject-specific log-OR)
#        GEE_binary   -> truth$alpha3_marginal     (population-averaged log-OR,
#                                                   pre-computed via numerical
#                                                   integration in 01_generate_data.R)
#      Evaluating GEE_binary against alpha3_conditional would be a
#      methodological error: GEE would appear biased simply because it
#      consistently estimates a different, but equally valid, quantity.
#      This mapping is implemented in .truth_for_method() below, which is the
#      single authoritative location for the estimand assignment logic.
#
#   2. PERFORMANCE MEASURES FOLLOW MORRIS, WHITE & CROWTHER (2019)
#      All metric formulas, definitions, and Monte Carlo Standard Error (MCSE)
#      expressions follow the reference simulation-study tutorial:
#        Morris TP, White IR, Crowther MJ (2019). Using simulation studies to
#        evaluate statistical methods. Statistics in Medicine, 38(11), 2074-2102.
#      This is the current methodological standard for reporting Monte Carlo
#      simulation results in statistics journals and is the appropriate
#      reference for referee enquiries about metric choice.
#
#   3. COVERAGE USES A WALD 95% CI WITH z = 1.96 FOR ALL METHODS
#      For LME, the technically correct interval uses the Satterthwaite t
#      distribution, but for a methods-comparison simulation study the
#      marginal difference from the z-approximation at N >= 100 is
#      negligible, and using a uniform z = 1.96 across all four methods
#      ensures that coverage differences reflect genuine inferential
#      differences rather than differences in the CI formula. This choice
#      is documented explicitly and is standard practice in simulation
#      papers comparing LME and GEE (e.g. Fitzmaurice, Laird & Ware 2011).
#
#   4. SE RATIO DEFINITION
#      SE Ratio = mean(reported SE) / empirical SE
#      A value of 1 means the average reported SE perfectly matches the true
#      sampling variability of the estimator (as measured by the Monte Carlo
#      empirical SD). Values < 1 indicate anti-conservative SE estimation;
#      values > 1 indicate conservative SE estimation.
#      For GEE, the PRIMARY reported SE is the ROBUST (sandwich) SE -- this
#      is the inferentially relevant quantity. The NAIVE SE is stored
#      separately and reported as "se_ratio_naive" solely to quantify how
#      far the model-based estimator is from the robust one, which is
#      itself a secondary diagnostic metric.
#      For LME and GLMM, there is only one SE (model-based), so
#      se_ratio_naive = se_ratio.
#
#   5. RELATIVE BIAS IS UNDEFINED UNDER THE NULL
#      When the true treatment effect is zero (null scenario), dividing bias
#      by |truth| produces an undefined (Inf or NaN) quantity. Relative bias
#      is therefore computed only when |truth| > 0 (alternative scenario) and
#      set explicitly to NA for the null scenario, preventing silent Inf values
#      from propagating into downstream tables and plots.
#
#   6. POWER vs. TYPE I ERROR OCCUPIES THE SAME COLUMN
#      Under the alternative scenario, mean(pvalue < alpha) estimates POWER.
#      Under the null scenario, mean(pvalue < alpha) estimates TYPE I ERROR.
#      Both quantities are stored in a single column named "rejection_rate"
#      with a separate "rejection_rate_label" column ("Power" / "Type I Error")
#      so that 06_tables_figures.R can display them appropriately without
#      duplicating any downstream logic.
#
#   7. MONTE CARLO STANDARD ERRORS (MCSE)
#      MCSEs quantify the simulation uncertainty in each estimated metric
#      and MUST be reported in a methodological paper to demonstrate that the
#      number of replicates (B = 1000) is sufficient. The formulas used are:
#        MCSE(bias)          = empirical_se / sqrt(B)
#        MCSE(empirical_se)  = empirical_se / sqrt(2*(B-1))
#        MCSE(rmse)          = sqrt(Var((theta_hat - truth)^2)) / (2*rmse*sqrt(B))
#        MCSE(coverage)      = sqrt(coverage*(1-coverage) / B)
#        MCSE(rejection_rate)= sqrt(rejection_rate*(1-rejection_rate) / B)
#      All from Morris et al. (2019), Table I.
#
#   8. FAILED/NON-CONVERGED REPLICATES
#      Replicates flagged as failed (failed=TRUE) or non-converged
#      (converged=FALSE) are EXCLUDED from all metric computations that
#      require a valid estimate or SE. They are counted separately and
#      reported in the convergence_rate and (where applicable) singularity_rate
#      metrics. This is the standard approach: including failed replicates
#      in, e.g., bias estimates would confound the model's statistical
#      performance with its computational reliability.
#
# DEPENDENCIES
#   data.table (manipulation)
#   No modelling packages required: this file only post-processes results.
#
# UPSTREAM DEPENDENCY (assumed already sourced by run_pipeline.R)
#   03_fit_models.R -- defines the structure of the fit result lists consumed
#                      here; no functions from that file are called directly.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 1. CONSTANTS
# ------------------------------------------------------------------------------

#' Nominal significance level for hypothesis test rejection decisions.
ALPHA_LEVEL <- 0.05

#' Methods whose singularity rate should be reported (random-effects models).
SINGULAR_METHODS <- c("LME", "GLMM")

# ------------------------------------------------------------------------------
# 2. ESTIMAND-TRUTH MAPPING
# ------------------------------------------------------------------------------

#' Return the correct scalar truth for a given method and truth list.
#'
#' @param method character; one of "LME", "GEE_continuous", "GLMM",
#'   "GEE_binary".
#' @param truth list; the $truth element from a fit_all_models() output,
#'   containing beta3, alpha3_conditional, and alpha3_marginal.
#' @return numeric scalar; the true parameter value the method estimates.
#' @details This is the single authoritative location for the estimand
#'   assignment logic described in the file header (Section 1). All other
#'   functions in this file call this helper; none hard-code a truth value
#'   directly.
.truth_for_method <- function(method, truth) {
  switch(method,
         LME            = truth$beta3,
         GEE_continuous = truth$beta3,
         GLMM           = truth$alpha3_conditional,
         GEE_binary     = truth$alpha3_marginal,
         stop(sprintf("Unknown method '%s' in .truth_for_method().", method))
  )
}

# ------------------------------------------------------------------------------
# 3. EXTRACT PER-REPLICATE VECTORS FOR ONE METHOD
# ------------------------------------------------------------------------------

#' Extract per-replicate scalar quantities for one method from a list of B
#' replicate fit results.
#'
#' @param replicate_results list of length B; each element is the output of
#'   fit_all_models() from 03_fit_models.R (contains $fits, $truth, $params).
#' @param method character; the method name to extract (must be a key of
#'   each replicate_results[[b]]$fits list).
#' @return a data.table with B rows (one per replicate) and columns:
#'   b (replicate index), estimate, se_robust, se_naive, pvalue, converged,
#'   singular (NA if not applicable), failed, truth_value.
#'   Rows where failed=TRUE are retained in the table but excluded from
#'   metric calculations in .compute_metrics_from_vectors().
.extract_replicate_vectors <- function(replicate_results, method) {
  B <- length(replicate_results)
  
  rows <- vector("list", B)
  for (b in seq_len(B)) {
    fit   <- replicate_results[[b]]$fits[[method]]
    truth <- replicate_results[[b]]$truth
    
    rows[[b]] <- data.table(
      b           = b,
      estimate    = if (isTRUE(fit$failed)) NA_real_ else as.numeric(fit$estimate),
      se_robust   = if (isTRUE(fit$failed)) NA_real_ else as.numeric(fit$se_robust),
      se_naive    = if (isTRUE(fit$failed)) NA_real_ else as.numeric(fit$se_naive),
      pvalue      = if (isTRUE(fit$failed)) NA_real_ else as.numeric(fit$pvalue),
      converged   = isTRUE(fit$converged),
      singular    = if (method %in% SINGULAR_METHODS) isTRUE(fit$singular) else NA,
      failed      = isTRUE(fit$failed),
      truth_value = .truth_for_method(method, truth)
    )
  }
  
  rbindlist(rows)
}

# ------------------------------------------------------------------------------
# 4. COMPUTE ALL METRICS FROM PER-REPLICATE VECTORS
# ------------------------------------------------------------------------------

#' Compute the full set of pre-specified performance measures from per-
#' replicate vectors for one method in one simulation cell.
#'
#' @param rv data.table; output of .extract_replicate_vectors() for one method.
#' @param scenario character; "null" or "alternative" (determines the
#'   rejection_rate_label and whether relative bias is defined).
#' @return a named list of scalar metrics and their MCSEs.
.compute_metrics_from_vectors <- function(rv, scenario) {
  
  B        <- nrow(rv)
  truth    <- rv$truth_value[1]   # constant within a cell
  
  # -- Subsets --------------------------------------------------------------
  # Replicates with a valid estimate/SE (not failed AND converged)
  rv_valid <- rv[!failed & converged]
  B_valid  <- nrow(rv_valid)
  
  if (B_valid == 0) {
    warning("No valid replicates found for this method/cell combination.",
            call. = FALSE)
  }
  
  # -- Core estimate quantities (valid replicates only) ---------------------
  estimates <- rv_valid$estimate
  se_robust <- rv_valid$se_robust
  se_naive  <- rv_valid$se_naive
  pvalues   <- rv_valid$pvalue
  
  mean_est     <- if (B_valid > 0) mean(estimates)         else NA_real_
  empirical_se <- if (B_valid > 1) sd(estimates)           else NA_real_
  avg_se_robust<- if (B_valid > 0) mean(se_robust)         else NA_real_
  avg_se_naive <- if (B_valid > 0) mean(se_naive)          else NA_real_
  
  # -- Bias -----------------------------------------------------------------
  bias <- if (!is.na(mean_est)) mean_est - truth else NA_real_
  
  # Relative bias: undefined when truth = 0 (null scenario)
  rel_bias <- if (!is.na(bias) && abs(truth) > .Machine$double.eps) {
    bias / abs(truth)
  } else {
    NA_real_
  }
  
  # -- RMSE -----------------------------------------------------------------
  errors_sq <- (estimates - truth)^2
  rmse      <- if (B_valid > 0) sqrt(mean(errors_sq)) else NA_real_
  
  # -- Coverage (Wald 95% CI, z = 1.96 for all methods; see file header) ---
  z_crit  <- qnorm(1 - ALPHA_LEVEL / 2)
  ci_lo   <- rv_valid$estimate - z_crit * rv_valid$se_robust
  ci_hi   <- rv_valid$estimate + z_crit * rv_valid$se_robust
  covered <- (ci_lo <= truth) & (truth <= ci_hi)
  coverage <- if (B_valid > 0) mean(covered) else NA_real_
  
  # -- SE Ratio -------------------------------------------------------------
  se_ratio       <- if (!is.na(avg_se_robust) && !is.na(empirical_se) && empirical_se > 0) {
    avg_se_robust / empirical_se
  } else {
    NA_real_
  }
  se_ratio_naive <- if (!is.na(avg_se_naive) && !is.na(empirical_se) && empirical_se > 0) {
    avg_se_naive / empirical_se
  } else {
    NA_real_
  }
  
  # -- Rejection rate (Power under alternative; Type I Error under null) ----
  rejection_rate <- if (B_valid > 0) mean(pvalues < ALPHA_LEVEL, na.rm = TRUE) else NA_real_
  rejection_label <- if (scenario == "null") "Type I Error" else "Power"
  
  # -- Convergence and singularity rates (all B replicates, not just valid) -
  convergence_rate <- mean(rv$converged)
  singularity_rate <- if (all(is.na(rv$singular))) NA_real_ else mean(rv$singular, na.rm = TRUE)
  
  # -- Monte Carlo Standard Errors (Morris et al. 2019, Table I) -----------
  mcse_bias <- if (!is.na(empirical_se)) empirical_se / sqrt(B_valid) else NA_real_
  
  mcse_empirical_se <- if (!is.na(empirical_se) && B_valid > 1) {
    empirical_se / sqrt(2 * (B_valid - 1))
  } else {
    NA_real_
  }
  
  mcse_rmse <- if (B_valid > 0 && !is.na(rmse) && rmse > 0) {
    sqrt(var(errors_sq)) / (2 * rmse * sqrt(B_valid))
  } else {
    NA_real_
  }
  
  mcse_coverage <- if (!is.na(coverage)) {
    sqrt(coverage * (1 - coverage) / B_valid)
  } else {
    NA_real_
  }
  
  mcse_rejection <- if (!is.na(rejection_rate)) {
    sqrt(rejection_rate * (1 - rejection_rate) / B_valid)
  } else {
    NA_real_
  }
  
  # -- Assemble result list -------------------------------------------------
  list(
    # Counts
    n_replicates      = B,
    n_valid           = B_valid,
    
    # Point estimate summary
    mean_estimate     = mean_est,
    truth             = truth,
    
    # Core metrics
    bias              = bias,
    relative_bias     = rel_bias,
    rmse              = rmse,
    coverage          = coverage,
    empirical_se      = empirical_se,
    avg_se_robust     = avg_se_robust,
    avg_se_naive      = avg_se_naive,
    se_ratio          = se_ratio,          # avg_se_robust / empirical_se
    se_ratio_naive    = se_ratio_naive,    # avg_se_naive  / empirical_se
    
    # Rejection rate (label distinguishes Power vs. Type I Error)
    rejection_rate      = rejection_rate,
    rejection_rate_label = rejection_label,
    
    # Computational diagnostics
    convergence_rate  = convergence_rate,
    singularity_rate  = singularity_rate,  # NA for GEE
    
    # MCSEs
    mcse_bias         = mcse_bias,
    mcse_empirical_se = mcse_empirical_se,
    mcse_rmse         = mcse_rmse,
    mcse_coverage     = mcse_coverage,
    mcse_rejection    = mcse_rejection
  )
}

# ------------------------------------------------------------------------------
# 5. TOP-LEVEL: COMPUTE METRICS FOR ALL METHODS IN ONE SIMULATION CELL
# ------------------------------------------------------------------------------

#' Compute performance measures for all four methods for a single simulation
#' design cell (fixed n, target_dropout_rate, and scenario).
#'
#' @param replicate_results list of length B; each element is the output of
#'   fit_all_models() from 03_fit_models.R.
#' @param n integer; number of subjects in this cell (for labelling).
#' @param target_dropout_rate numeric; the target dropout rate (for labelling).
#' @param scenario character; "null" or "alternative" (for labelling and
#'   relative-bias / rejection-rate logic).
#' @return a data.table with one row per method, all metric columns, and
#'   cell-identifying columns (n, target_dropout_rate, scenario).
compute_cell_metrics <- function(replicate_results,
                                 n,
                                 target_dropout_rate,
                                 scenario = c("null", "alternative")) {
  scenario <- match.arg(scenario)
  
  if (length(replicate_results) == 0) {
    stop("replicate_results must contain at least one element.")
  }
  
  methods <- c("LME", "GEE_continuous", "GLMM", "GEE_binary")
  
  method_rows <- lapply(methods, function(method) {
    rv      <- .extract_replicate_vectors(replicate_results, method)
    metrics <- .compute_metrics_from_vectors(rv, scenario)
    
    # Convert to a one-row data.table
    as.data.table(c(list(method = method), metrics))
  })
  
  result <- rbindlist(method_rows)
  
  # Prepend cell-identifying columns
  result[, `:=`(
    n                   = n,
    target_dropout_rate = target_dropout_rate,
    scenario            = scenario
  )]
  
  # Reorder: identifying columns first
  setcolorder(result, c("scenario", "n", "target_dropout_rate", "method"))
  
  result[]
}

# ------------------------------------------------------------------------------
# 6. AGGREGATE ACROSS ALL CELLS (CONVENIENCE WRAPPER)
# ------------------------------------------------------------------------------

#' Aggregate metrics across all simulation cells from the full results
#' object produced by 05_simulation_engine.R.
#'
#' @param all_cell_results named list; produced by run_simulation_grid() in
#'   05_simulation_engine.R. Each element corresponds to one (scenario, n,
#'   target_dropout_rate) cell and contains a $replicate_results list and
#'   cell-level metadata ($n, $target_dropout_rate, $scenario).
#' @return a single data.table with all cells and all methods stacked,
#'   suitable for direct use by 06_tables_figures.R.
compute_all_metrics <- function(all_cell_results) {
  if (length(all_cell_results) == 0) {
    stop("all_cell_results is empty.")
  }
  
  cell_tables <- lapply(all_cell_results, function(cell) {
    compute_cell_metrics(
      replicate_results   = cell$replicate_results,
      n                   = cell$n,
      target_dropout_rate = cell$target_dropout_rate,
      scenario            = cell$scenario
    )
  })
  
  rbindlist(cell_tables)
}

# ------------------------------------------------------------------------------
# 7. SMOKE TEST (manual, interactive use only -- not sourced by the pipeline)
# ------------------------------------------------------------------------------
# Uncomment to verify locally (requires 01-03 already sourced):
#
# source("01_generate_data.R")
# source("02_apply_dropout.R")
# source("03_fit_models.R")
# source("04_metrics.R")
#
# slopes <- get_dropout_hazard_slopes()
# calib  <- build_dropout_calibration_table(n_calib = 5000)
#
# set.seed(1)
# reps <- lapply(1:20, function(b) {
#   rep_obj <- generate_replicate_data(n = 200, scenario = "alternative")
#   do_obj  <- apply_monotone_dropout(rep_obj, 0.20, calib, slopes)
#   fit_all_models(do_obj)
# })
#
# metrics <- compute_cell_metrics(reps, n = 200, target_dropout_rate = 0.20,
#                                  scenario = "alternative")
# print(metrics[, .(method, bias, rmse, coverage, rejection_rate,
#                   convergence_rate, mcse_bias, mcse_coverage)])