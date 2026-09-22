# ==============================================================================
# 02_apply_dropout.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Applies a MONOTONE, SUBJECT-LEVEL, MISSING-AT-RANDOM (MAR) dropout
#   mechanism to the fully observed data produced by 01_generate_data.R.
#   Dropout probability at each post-baseline visit depends on the subject's
#   own PREVIOUSLY OBSERVED continuous outcome (y1), the current visit time,
#   and treatment group -- never on the (unobserved, post-dropout) outcome
#   itself, which is what makes the mechanism MAR rather than MNAR.
#
#   The hazard-model intercept is NOT chosen arbitrarily. It is calibrated
#   numerically (via root-finding) so that the realized subject-level dropout
#   proportion matches each target rate (10%, 20%, 30%, 40%) to high
#   precision, separately for each effect-size scenario (null/alternative),
#   since the y1 distribution -- and hence the dropout hazard -- differs
#   slightly between scenarios.
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#
#   1. MONOTONE, SUBJECT-LEVEL MISSINGNESS
#      Dropout is modeled as a discrete-time survival process: at each visit
#      t = 1, 2, 3, a subject who is still in the study faces a hazard of
#      dropping out. Once a subject drops out at visit t, ALL visits >= t
#      become missing for BOTH outcomes (y1 and y2), and the subject is never
#      "at risk" again. This is the standard, clinically realistic monotone
#      missing-data pattern in longitudinal RCTs and is exactly what LME/GEE/
#      GLMM under MAR are designed to handle correctly (intermittent,
#      non-monotone missingness is a materially different problem and is
#      intentionally NOT simulated here, per protocol).
#
#   2. MAR MECHANISM
#      The hazard of dropping out at visit t depends only on (a) the
#      subject's own observed y1 value at visit t-1, (b) time t, and (c)
#      treatment group -- all of which are FULLY OBSERVED for any subject
#      still in the study at visit t-1. No dependence on the unobserved
#      future outcome, and no dependence on y2 (to keep the missingness
#      mechanism for the binary outcome a function of a different, but
#      observed, variable -- a common and realistic real-world scenario,
#      e.g. a continuous biomarker driving discontinuation that also
#      happens to be correlated, through the shared visit schedule, with a
#      binary clinical endpoint). This satisfies Rubin's MAR definition:
#      P(missing | observed, unobserved) = P(missing | observed).
#
#   3. CALIBRATION VIA ROOT-FINDING WITH COMMON RANDOM NUMBERS
#      The relationship between the hazard-model intercept and the realized
#      dropout proportion is monotonic in expectation, but a naive Monte
#      Carlo evaluation of "achieved dropout rate" at two different
#      intercepts is NOISY and not exactly monotonic for a finite sample,
#      which breaks bisection-based root-finders such as uniroot(). We
#      resolve this with the COMMON RANDOM NUMBERS (CRN) variance-reduction
#      technique: a single matrix of Uniform(0,1) draws is generated ONCE
#      per calibration run and reused across every candidate intercept value
#      evaluated by uniroot(). Because the dropout probability is a strictly
#      increasing function of the intercept, and a subject drops out iff
#      U < P(intercept), CRN guarantees the achieved dropout rate is EXACTLY
#      (not just approximately) monotonically non-decreasing in the
#      intercept, making uniroot() converge reliably and reproducibly to a
#      precise calibrated value. A large pseudo-population (n_calib) is also
#      used to keep finite-sample noise small.
#
#   4. WHY CALIBRATION DOES NOT NEED TO BE REPEATED PER N
#      The subject-level dropout probability model is a property of the
#      hazard parameters and the (N-invariant) distribution of y1; it does
#      not depend on the number of subjects in a given simulation cell.
#      Calibration is therefore performed once per (scenario, target rate)
#      pair at a large pseudo-population size and the resulting intercept is
#      reused, unchanged, for N = 100, 200, and 500 and across all 1000
#      Monte Carlo replicates -- avoiding 12x (3 N) redundant calibration
#      work and guaranteeing that all N-levels target the identical
#      population dropout probability.
#
#   5. PREVIOUS-OUTCOME STANDARDIZATION
#      The previous-y1 predictor is standardized (centered/scaled by fixed,
#      data-generating-process-informed constants) before entering the
#      hazard's linear predictor. This is purely for numerical conditioning
#      of the logistic hazard model (avoiding extreme intercepts on the raw
#      ~45-55 scale of y1) and has no effect on the calibrated achieved
#      dropout rate, which is determined empirically by the root-finding
#      procedure regardless of this choice.
#
# DEPENDENCIES
#   data.table   (manipulation)
#   stats        (uniroot, rbinom, runif -- base R)
#
# UPSTREAM DEPENDENCY (assumed already sourced by run_pipeline.R)
#   01_generate_data.R   -- provides generate_replicate_data(), SIM_TIME_GRID
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 1. DROPOUT HAZARD MODEL: FIXED PARAMETERS
# ------------------------------------------------------------------------------

#' Fixed coefficients for the dropout hazard's MAR predictors.
#'
#' @return a named list of hazard-model slopes and y1-standardization
#'   constants. The intercept is deliberately excluded -- it is the sole
#'   quantity calibrated by root-finding to hit each target dropout rate.
#' @details
#'   gamma_y1_std : log-odds change in dropout hazard per 1-SD increase in
#'                  the subject's previous y1 (y1 is parameterized in
#'                  01_generate_data.R as a severity-type score where higher
#'                  values are clinically worse, so a positive coefficient
#'                  means subjects with worse previous status are more
#'                  likely to drop out -- a realistic informative-looking
#'                  but still MAR mechanism, since it depends only on the
#'                  observed previous value).
#'   gamma_time    : log-odds change in dropout hazard per unit increase in
#'                   visit time (captures general trial-fatigue / cumulative
#'                   discontinuation risk independent of outcome).
#'   gamma_group   : log-odds change in dropout hazard for active treatment
#'                   vs. control (small negative value: subjects on active
#'                   treatment are modestly less likely to discontinue,
#'                   a common pattern attributable to perceived benefit).
#'   y1_center,
#'   y1_scale      : fixed standardization constants for the previous-y1
#'                   predictor, chosen to be representative of the y1
#'                   distribution implied by 01_generate_data.R's parameters
#'                   (baseline mean 50, moderate decline, total SD on the
#'                   order of 8-10 once between- and within-subject variance
#'                   are combined). These are conditioning constants only.
get_dropout_hazard_slopes <- function() {
  list(
    gamma_y1_std = 0.40,
    gamma_time   = 0.15,
    gamma_group  = -0.10,
    y1_center    = 45,
    y1_scale     = 8
  )
}

# ------------------------------------------------------------------------------
# 2. HAZARD PROBABILITY
# ------------------------------------------------------------------------------

#' Compute the per-visit dropout hazard probability under the MAR model.
#'
#' @param prev_y1 numeric vector; previously OBSERVED y1 value(s).
#' @param time numeric vector; current visit time (the visit at risk).
#' @param group numeric/integer vector; treatment group (0/1).
#' @param intercept numeric scalar; the (calibrated) hazard intercept.
#' @param slopes list; output of get_dropout_hazard_slopes().
#' @return numeric vector of dropout probabilities in (0, 1).
dropout_hazard_probability <- function(prev_y1, time, group, intercept, slopes) {
  prev_y1_std <- (prev_y1 - slopes$y1_center) / slopes$y1_scale
  
  eta <- intercept +
    slopes$gamma_y1_std * prev_y1_std +
    slopes$gamma_time   * time +
    slopes$gamma_group  * group
  
  plogis(eta)
}

# ------------------------------------------------------------------------------
# 3. CORE MONOTONE DROPOUT SIMULATOR
# ------------------------------------------------------------------------------

#' Simulate the first dropout visit for every subject under a monotone,
#' sequential discrete-time hazard process.
#'
#' @param dt data.table; long-format data with columns id, time, group, y1
#'   (as produced by generate_replicate_data()$data). Must contain exactly
#'   the visits in SIM_TIME_GRID for every subject, sorted by id, time.
#' @param intercept numeric scalar; hazard intercept (calibrated or candidate).
#' @param slopes list; output of get_dropout_hazard_slopes().
#' @param u_matrix optional numeric matrix (n_subjects x 3) of Uniform(0,1)
#'   draws, one column per post-baseline visit (1, 2, 3). If supplied, these
#'   EXACT draws are used for the dropout comparison (Common Random Numbers
#'   -- used during calibration to guarantee monotonicity of the achieved
#'   dropout rate in `intercept`). If NULL, fresh draws are generated via
#'   runif() (used for actual Monte Carlo replicate generation).
#' @return a data.table keyed by id with columns: id, dropout_visit (1, 2, 3,
#'   or Inf if the subject completes the study without dropping out).
#' @details Visits are processed sequentially (t = 1, 2, 3) because the
#'   hazard at visit t depends on the y1 value observed at visit t-1, which
#'   is only well-defined for subjects who have not already dropped out.
#'   Each iteration is fully vectorized across subjects.
simulate_monotone_dropout <- function(dt, intercept, slopes, u_matrix = NULL) {
  ids <- sort(unique(dt$id))
  n <- length(ids)
  post_baseline_visits <- SIM_TIME_GRID[SIM_TIME_GRID > 0]  # c(1, 2, 3)
  n_visits <- length(post_baseline_visits)
  
  if (!is.null(u_matrix)) {
    if (!all(dim(u_matrix) == c(n, n_visits))) {
      stop(sprintf("u_matrix must be %d x %d (subjects x post-baseline visits).",
                   n, n_visits))
    }
  }
  
  # Wide matrix of y1 by visit (rows = subjects in `ids` order, columns =
  # SIM_TIME_GRID order), used to look up "previous observed y1" in O(1).
  y1_wide <- dcast(dt, id ~ time, value.var = "y1")
  setorder(y1_wide, id)
  if (!identical(y1_wide$id, ids)) {
    stop("Internal error: y1_wide subject ordering does not match `ids`.")
  }
  y1_mat <- as.matrix(y1_wide[, -"id"])  # columns ordered 0, 1, 2, 3 by dcast
  
  group_by_id <- unique(dt[, .(id, group)])
  setorder(group_by_id, id)
  group_vec <- group_by_id$group
  
  already_dropped <- rep(FALSE, n)
  dropout_visit <- rep(Inf, n)
  
  for (j in seq_len(n_visits)) {
    t_visit <- post_baseline_visits[j]
    prev_y1_all <- y1_mat[, j]  # column j of y1_mat corresponds to time = t_visit - 1
    
    eligible <- !already_dropped
    if (!any(eligible)) break
    
    haz_p <- dropout_hazard_probability(
      prev_y1 = prev_y1_all[eligible],
      time    = t_visit,
      group   = group_vec[eligible],
      intercept = intercept,
      slopes = slopes
    )
    
    u_draws <- if (is.null(u_matrix)) {
      runif(sum(eligible))
    } else {
      u_matrix[eligible, j]
    }
    
    newly_dropped_local <- u_draws < haz_p
    
    eligible_idx <- which(eligible)
    dropout_visit[eligible_idx[newly_dropped_local]] <- t_visit
    already_dropped[eligible_idx[newly_dropped_local]] <- TRUE
  }
  
  data.table(id = ids, dropout_visit = dropout_visit)
}

#' Compute the realized (achieved) subject-level dropout proportion.
#'
#' @param dropout_visit numeric vector as returned by
#'   simulate_monotone_dropout()$dropout_visit.
#' @param max_visit numeric scalar; the last scheduled visit (default 3,
#'   matching SIM_TIME_GRID). A subject is counted as a "dropout" if they
#'   leave the study at or before this visit (i.e. fail to complete the
#'   full protocol), which includes dropping out exactly at the final visit.
#' @return numeric scalar in [0, 1].
compute_achieved_dropout_rate <- function(dropout_visit, max_visit = max(SIM_TIME_GRID)) {
  mean(dropout_visit <= max_visit)
}

# ------------------------------------------------------------------------------
# 4. CALIBRATION: ROOT-FINDING FOR THE HAZARD INTERCEPT
# ------------------------------------------------------------------------------

#' Calibrate the dropout hazard intercept so the achieved subject-level
#' dropout rate matches a target proportion, using Common Random Numbers to
#' guarantee a numerically monotonic objective function for uniroot().
#'
#' @param target_rate numeric scalar in (0, 1); desired dropout proportion.
#' @param scenario character; "null" or "alternative" (passed to
#'   generate_replicate_data() since the y1 distribution differs slightly
#'   by scenario through beta3).
#' @param slopes list; output of get_dropout_hazard_slopes().
#' @param n_calib integer; pseudo-population size used for calibration
#'   (large enough to make residual finite-sample noise negligible at the
#'   precision required; calibration is a one-time cost, not incurred per
#'   Monte Carlo replicate).
#' @param seed integer; RNG seed for the calibration pseudo-population and
#'   the CRN matrix, ensuring fully reproducible calibrated intercepts.
#' @param intercept_interval numeric length-2 vector; initial bracket passed
#'   to uniroot() (auto-extended via extendInt if the root falls outside it).
#' @return a list: intercept (calibrated value), achieved_rate (achieved
#'   dropout rate at the calibrated intercept, for verification), target_rate.
calibrate_dropout_intercept <- function(target_rate,
                                        scenario = c("null", "alternative"),
                                        slopes = get_dropout_hazard_slopes(),
                                        n_calib = 20000,
                                        seed = 2024,
                                        intercept_interval = c(-15, 15)) {
  scenario <- match.arg(scenario)
  
  if (target_rate <= 0 || target_rate >= 1) {
    stop("target_rate must lie strictly within (0, 1).")
  }
  
  set.seed(seed)
  calib_replicate <- generate_replicate_data(n = n_calib, scenario = scenario)
  calib_dt <- calib_replicate$data
  
  n_post_visits <- length(SIM_TIME_GRID) - 1
  u_matrix <- matrix(runif(n_calib * n_post_visits), nrow = n_calib, ncol = n_post_visits)
  
  objective <- function(intercept) {
    dropout_result <- simulate_monotone_dropout(calib_dt, intercept, slopes, u_matrix = u_matrix)
    compute_achieved_dropout_rate(dropout_result$dropout_visit) - target_rate
  }
  
  root <- uniroot(objective, interval = intercept_interval, extendInt = "yes", tol = 1e-6)
  
  achieved <- objective(root$root) + target_rate
  
  list(
    intercept     = root$root,
    achieved_rate = achieved,
    target_rate   = target_rate,
    scenario      = scenario
  )
}

#' Build the full calibration lookup table across all (scenario, target
#' dropout rate) combinations required by the study protocol.
#'
#' @param scenarios character vector; defaults to c("null", "alternative").
#' @param target_rates numeric vector; defaults to c(0.10, 0.20, 0.30, 0.40).
#' @param slopes list; output of get_dropout_hazard_slopes().
#' @param n_calib integer; calibration pseudo-population size.
#' @param seed integer; base RNG seed (offset per cell for reproducible,
#'   non-degenerate CRN matrices across cells).
#' @return a data.table with columns: scenario, target_rate, intercept,
#'   achieved_rate. Intended to be computed ONCE (e.g. in run_pipeline.R or
#'   at the top of 05_simulation_engine.R) and reused for every N and every
#'   Monte Carlo replicate.
build_dropout_calibration_table <- function(scenarios = c("null", "alternative"),
                                            target_rates = c(0.10, 0.20, 0.30, 0.40),
                                            slopes = get_dropout_hazard_slopes(),
                                            n_calib = 20000,
                                            seed = 2024) {
  grid <- CJ(scenario = scenarios, target_rate = target_rates, sorted = FALSE)
  
  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    cell_seed <- seed + i  # distinct, reproducible seed per calibration cell
    calib <- calibrate_dropout_intercept(
      target_rate = grid$target_rate[i],
      scenario    = grid$scenario[i],
      slopes      = slopes,
      n_calib     = n_calib,
      seed        = cell_seed
    )
    results[[i]] <- data.table(
      scenario      = calib$scenario,
      target_rate   = calib$target_rate,
      intercept     = calib$intercept,
      achieved_rate = calib$achieved_rate
    )
  }
  
  rbindlist(results)
}

# ------------------------------------------------------------------------------
# 5. APPLY CALIBRATED DROPOUT TO A SINGLE MONTE CARLO REPLICATE
# ------------------------------------------------------------------------------

#' Apply the calibrated monotone MAR dropout mechanism to one complete-data
#' Monte Carlo replicate.
#'
#' @param replicate_obj list; output of generate_replicate_data() from
#'   01_generate_data.R (contains $data, $params, $truth).
#' @param target_rate numeric scalar; the target dropout rate for this
#'   simulation cell (must be one of the rates present in calibration_table).
#' @param calibration_table data.table; output of
#'   build_dropout_calibration_table().
#' @param slopes list; output of get_dropout_hazard_slopes() (must match
#'   the slopes used to build calibration_table).
#' @return a list:
#'   data                 - the long-format data.table with y1, y2 set to NA
#'                           from the subject's dropout visit onward, plus an
#'                           `observed` logical column,
#'   dropout_visit         - per-subject first-missing-visit table,
#'   achieved_dropout_rate - the realized dropout proportion in THIS
#'                           replicate (will fluctuate around target_rate by
#'                           Monte Carlo error at finite N -- this fluctuation
#'                           is itself a quantity worth summarizing across
#'                           replicates in 06_tables_figures.R),
#'   target_dropout_rate   - target_rate (echoed for downstream bookkeeping),
#'   truth, params         - passed through unchanged from replicate_obj.
apply_monotone_dropout <- function(replicate_obj, target_rate, calibration_table, slopes) {
  scen_val <- replicate_obj$params$scenario   
  rate_val <- target_rate                     
  
  calib_row <- calibration_table[scenario == scen_val & target_rate == rate_val]
  if (nrow(calib_row) != 1) {
    stop(sprintf(
      "Expected exactly one calibration row for scenario='%s', target_rate=%.2f; found %d.",
      scen_val, rate_val, nrow(calib_row)
    ))
  }
  intercept <- calib_row$intercept
  
  dt <- copy(replicate_obj$data)
  dropout_result <- simulate_monotone_dropout(dt, intercept, slopes, u_matrix = NULL)
  
  dt <- merge(dt, dropout_result, by = "id", sort = FALSE)
  setorder(dt, id, time)
  
  dt[, observed := time < dropout_visit]
  dt[observed == FALSE, `:=`(y1 = NA_real_, y2 = NA_integer_)]
  
  list(
    data = dt,
    dropout_visit = dropout_result,
    achieved_dropout_rate = compute_achieved_dropout_rate(dropout_result$dropout_visit),
    target_dropout_rate = target_rate,
    truth = replicate_obj$truth,
    params = replicate_obj$params
  )
}

# ------------------------------------------------------------------------------
# 6. SMOKE TEST (manual, interactive use only -- not sourced by the pipeline)
# ------------------------------------------------------------------------------
# Uncomment to verify locally (requires 01_generate_data.R already sourced):
#
# source("01_generate_data.R")
# source("02_apply_dropout.R")
#
# slopes <- get_dropout_hazard_slopes()
# calib_table <- build_dropout_calibration_table(n_calib = 20000)
# print(calib_table)   # achieved_rate should closely match target_rate
#
# set.seed(99)
# rep_obj <- generate_replicate_data(n = 500, scenario = "alternative")
# out <- apply_monotone_dropout(rep_obj, target_rate = 0.30,
#                                calibration_table = calib_table, slopes = slopes)
# cat("Achieved dropout rate (this replicate):", out$achieved_dropout_rate, "\n")
# print(out$data[id %in% 1:3])   # inspect monotone NA pattern for a few subjects