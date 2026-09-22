# ==============================================================================
# 03_fit_models.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Fits all four models specified by the study protocol to a single
#   (possibly incomplete, post-dropout) Monte Carlo replicate dataset:
#
#     Continuous outcome (y1)
#       M1: Linear Mixed Model          (LME)   -- lme4 / lmerTest
#       M2: Generalized Estimating Equations     -- geepack (gaussian)
#
#     Binary outcome (y2)
#       M3: Generalized Linear Mixed Model (GLMM) -- lme4 (binomial)
#       M4: Generalized Estimating Equations       -- geepack (binomial)
#
#   Each fitting function returns a STANDARDIZED result list so that
#   04_metrics.R can consume all four outputs identically, without any
#   model-specific branching.
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#
#   1. ML NOT REML FOR LME
#      The protocol specifies Maximum Likelihood (ML) estimation.
#      lme4 defaults to REML for variance components, which is preferred for
#      pure variance-component estimation but produces fixed-effect SEs that
#      are not directly comparable to the sandwich-based GEE SEs we are
#      evaluating alongside. Under ML (REML=FALSE), both point estimates and
#      their sampling distributions are on a common likelihood basis, making
#      cross-method comparisons in bias/RMSE/coverage well-defined.
#
#   2. EXCHANGEABLE WORKING CORRELATION FOR GEE
#      With only 4 visits per subject, an unstructured working correlation
#      matrix (6 free parameters) would be poorly identified at small N.
#      Exchangeable (compound symmetry) is the natural, parsimonious choice
#      for an RCT where no prior reason exists to assume that correlations
#      decay with visit-lag. It is also the most common choice in the
#      methodological GEE literature for this setting.
#
#   3. ROBUST (SANDWICH) SE AS PRIMARY GEE SE
#      GEE inference is valid under MAR ONLY with the robust (sandwich)
#      variance estimator -- the naive/model-based estimator is inconsistent
#      unless the working correlation structure is correctly specified, which
#      we cannot verify in any single dataset. The robust SE is therefore the
#      correct primary SE for GEE. The naive SE is stored separately and is
#      used exclusively for the SE Ratio metric in 04_metrics.R.
#
#   4. SATTERTHWAITE DEGREES OF FREEDOM FOR LME P-VALUES
#      lmerTest::summary() appends Satterthwaite df and p-values to lme4
#      models. Satterthwaite is the current consensus method for df
#      approximation in LME (preferred over Kenward-Roger for ML-fitted
#      models; Kenward-Roger was designed for REML). These p-values are used
#      for the Power / Type I Error metrics in 04_metrics.R.
#
#   5. CONVERGENCE AND SINGULARITY DETECTION
#      For lme4 models: convergence is assessed from the optimizer's own
#      diagnostic messages stored in fit@optinfo. Singularity is assessed
#      via lme4::isSingular(), which checks whether any variance component
#      or correlation parameter is at its boundary (0 for variances, +/-1
#      for correlations). Both flags are stored in the result list.
#      For GEE models: geepack stores an error code in fit$geese$error
#      (0 = converged). Singularity is not a meaningful concept for GEE
#      (no random effects) and is returned as NA.
#
#   6. TREATMENT EFFECT COEFFICIENT OF INTEREST
#      The treatment effect of interest is always the time:group interaction
#      term. The interaction term name in all model formula outputs under R's
#      default naming convention is "time:group". A single constant
#      (INTERACTION_TERM) is defined once at the top of this file and used
#      everywhere, ensuring that a future model re-parameterization (e.g.,
#      changing the formula) only needs to be updated in one place.
#
#   7. NA HANDLING
#      All four models receive the full long-format dataset with NAs already
#      applied in 02_apply_dropout.R. lme4 and geepack both use na.action =
#      na.omit by default: rows where the outcome is NA are silently excluded,
#      and each subject contributes all and only their observed visits.
#      This is the standard, correct approach for monotone-missing longitudinal
#      data under MAR.
#
# DEPENDENCIES
#   lme4, lmerTest, geepack   (must be installed)
#   data.table                (manipulation)
#
# UPSTREAM DEPENDENCY (assumed already sourced by run_pipeline.R)
#   01_generate_data.R   -- (for SIM_TIME_GRID; not called here directly)
#   02_apply_dropout.R   -- provides the data list consumed by fit_all_models()
# ==============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(geepack)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 1. CONSTANTS
# ------------------------------------------------------------------------------

#' Name of the treatment-effect interaction term as it appears in R model
#' coefficient output. Used everywhere to extract the focal coefficient.
INTERACTION_TERM <- "time:group"

# ------------------------------------------------------------------------------
# 2. SHARED HELPER: STANDARDIZED FAILED-FIT RESULT
# ------------------------------------------------------------------------------

#' Return a standardized NA-filled result list signalling a failed model fit.
#'
#' @param method character; one of "LME", "GEE_continuous", "GLMM",
#'   "GEE_binary".
#' @param reason character; brief description of why the fit failed.
#' @return named list with all result fields set to NA except converged=FALSE.
.failed_fit <- function(method, reason = "unknown") {
  warning(sprintf("[%s] Model fit failed: %s", method, reason), call. = FALSE)
  list(
    method        = method,
    estimate      = NA_real_,
    se_robust     = NA_real_,   # primary SE (robust/sandwich for GEE; model-based for LME/GLMM)
    se_naive      = NA_real_,   # naive SE (model-based; for SE Ratio metric only)
    pvalue        = NA_real_,
    converged     = FALSE,
    singular      = NA,
    failed        = TRUE,
    failure_reason = reason
  )
}

# ------------------------------------------------------------------------------
# 3. HELPER: CONVERGENCE DETECTION FOR lme4 MODELS
# ------------------------------------------------------------------------------

#' Extract convergence status from an lme4 model fit object.
#'
#' @param fit an object of class "lmerMod" or "glmerMod".
#' @return logical TRUE if the optimizer reported no convergence issues.
#' @details lme4 stores optimizer messages in fit@optinfo$conv$lme4$messages
#'   (post-fit convergence checks) and fit@optinfo$warnings (optimizer
#'   warnings). Both must be empty for the fit to be considered converged.
#'   CRITICAL FIX: fit@optinfo$warnings is version-unstable and can be
#'   non-empty even for perfectly converged models in recent lme4 releases.
#'   We now check ONLY fit@optinfo$conv$lme4$messages, which is the
#'   authoritative post-fit convergence signal from gradient/Hessian checks.
.lme4_converged <- function(fit) {
  msgs <- tryCatch(
    fit@optinfo$conv$lme4$messages,
    error = function(e) "slot access failed"
  )
  is.null(msgs) || length(msgs) == 0
}

# ------------------------------------------------------------------------------
# 3b. HELPER: CONVERGENCE DETECTION FOR geepack MODELS
# ------------------------------------------------------------------------------

#' Assess convergence of a geepack geeglm fit.
#'
#' @param fit geeglm object.
#' @return logical scalar.
#' @details geepack stores the optimizer error code in fit$geese$error
#'   (0 = converged). This helper guards against NULL or zero-length values
#'   that can arise in some geepack versions, where a bare == 0 comparison
#'   silently returns logical(0), causing isTRUE() to return FALSE and
#'   flagging every GEE fit as non-converged.
.gee_converged <- function(fit) {
  err <- tryCatch(fit$geese$error, error = function(e) NULL)
  length(err) == 1L && !is.na(err) && err == 0L
}

# ------------------------------------------------------------------------------
# 4. HELPER: COEFFICIENT EXTRACTION FROM lme4 SUMMARY
# ------------------------------------------------------------------------------

#' Extract estimate, model-based SE, and p-value for INTERACTION_TERM from
#' an lmerTest/lme4 summary table.
#'
#' @param fit lmerMod or glmerMod object.
#' @param method character; method label (used in error messages).
#' @return named numeric vector: estimate, se, pvalue.
.extract_lme4_coef <- function(fit, method) {
  coef_table <- tryCatch(
    coef(summary(fit)),
    error = function(e) NULL
  )
  if (is.null(coef_table) || !INTERACTION_TERM %in% rownames(coef_table)) {
    stop(sprintf("[%s] Coefficient '%s' not found in summary.", method, INTERACTION_TERM))
  }
  
  row <- coef_table[INTERACTION_TERM, , drop = FALSE]
  
  # Column naming differs between lmerTest (adds "df", "Pr(>|t|)") and
  # glmer (uses "Pr(>|z|)"). Robustly find the SE and p-value columns.
  se_col <- grep("^Std\\.", colnames(row), value = TRUE)
  pv_col <- grep("^Pr", colnames(row), value = TRUE)
  
  if (length(se_col) != 1 || length(pv_col) != 1) {
    stop(sprintf("[%s] Could not unambiguously identify SE/p-value columns.", method))
  }
  
  c(estimate = row[1, "Estimate"],
    se       = row[1, se_col],
    pvalue   = row[1, pv_col])
}

# ------------------------------------------------------------------------------
# 5. HELPER: COEFFICIENT EXTRACTION FROM geepack SUMMARY
# ------------------------------------------------------------------------------

#' Extract estimate, robust SE, naive SE, and p-value for INTERACTION_TERM
#' from a geeglm summary table.
#'
#' @param fit a "geeglm" object.
#' @param method character; method label (used in error messages).
#' @return named numeric vector: estimate, se_robust, se_naive, pvalue.
#' @details geepack reports the ROBUST (sandwich) SE as "Std.err" in the
#'   summary table. The naive (model-based) variance is stored separately in
#'   fit$geese$vbeta.naiv (the naive variance-covariance matrix). The p-value
#'   corresponds to the Wald chi-squared test using the robust SE.
.extract_gee_coef <- function(fit, method) {
  coef_table <- tryCatch(
    coef(summary(fit)),
    error = function(e) NULL
  )
  if (is.null(coef_table) || !INTERACTION_TERM %in% rownames(coef_table)) {
    stop(sprintf("[%s] Coefficient '%s' not found in summary.", method, INTERACTION_TERM))
  }
  
  row <- coef_table[INTERACTION_TERM, , drop = FALSE]
  
  # Robust SE (geepack labels this "Std.err")
  se_robust <- row[1, "Std.err"]
  
  # Naive SE: extracted from vbeta.naiv (the naive variance-covariance matrix)
  coef_names <- names(coef(fit))
  idx <- which(coef_names == INTERACTION_TERM)
  if (length(idx) != 1) {
    stop(sprintf("[%s] Could not locate '%s' in coefficient names for naive SE.", method, INTERACTION_TERM))
  }
  se_naive <- sqrt(fit$geese$vbeta.naiv[idx, idx])
  
  # P-value: geepack uses Wald statistic; column may be labeled "Pr(>|W|)"
  pv_col <- grep("^Pr", colnames(row), value = TRUE)
  if (length(pv_col) != 1) {
    stop(sprintf("[%s] Could not unambiguously identify p-value column.", method))
  }
  
  c(estimate  = row[1, "Estimate"],
    se_robust = se_robust,
    se_naive  = se_naive,
    pvalue    = row[1, pv_col])
}

# ------------------------------------------------------------------------------
# 6. MODEL M1: LINEAR MIXED MODEL (LME) for continuous outcome y1
# ------------------------------------------------------------------------------

#' Fit a Linear Mixed Model to the continuous outcome y1.
#'
#' @param dt data.table; long-format dataset with columns id, time, group,
#'   y1 (possibly with NAs from dropout).
#' @return standardized result list (see .failed_fit() for field definitions).
#' @details Model specification:
#'   y1 ~ time + group + time:group + (1 + time | id)
#'   Fitted by ML (REML = FALSE) using lmerTest::lmer(), which appends
#'   Satterthwaite df and p-values to the coefficient table. The time:group
#'   interaction is the treatment effect of interest.
fit_lme <- function(dt) {
  method <- "LME"
  
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        y1 ~ time + group + time:group + (1 + time | id),
        data   = dt,
        REML   = FALSE,
        control = lmerControl(optimizer = "bobyqa",
                              optCtrl  = list(maxfun = 2e5))
      ),
      warning = function(w) {
        # Let lme4 store the warning internally; do not suppress it here
        # so that .lme4_converged() can find it in fit@optinfo.
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      return(.failed_fit(method, reason = conditionMessage(e)))
    }
  )
  
  # If tryCatch returned a failed-fit list, pass it through
  if (!inherits(fit, "lmerMod")) return(fit)
  
  converged <- .lme4_converged(fit)
  singular  <- isSingular(fit)
  
  coefs <- tryCatch(
    .extract_lme4_coef(fit, method),
    error = function(e) return(.failed_fit(method, reason = conditionMessage(e)))
  )
  if (is.list(coefs)) return(coefs)  # failed-fit list from extraction error
  
  list(
    method        = method,
    estimate      = coefs["estimate"],
    se_robust     = coefs["se"],    # model-based SE; "robust" slot for uniform downstream access
    se_naive      = coefs["se"],    # same as model-based (no sandwich correction for LME)
    pvalue        = coefs["pvalue"],
    converged     = converged,
    singular      = singular,
    failed        = FALSE,
    failure_reason = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 7. MODEL M2: GEE (GAUSSIAN) for continuous outcome y1
# ------------------------------------------------------------------------------

#' Fit a GEE model with Gaussian family to the continuous outcome y1.
#'
#' @param dt data.table; long-format dataset with columns id, time, group,
#'   y1 (possibly with NAs from dropout). Must be sorted by id then time.
#' @return standardized result list.
#' @details Model specification:
#'   y1 ~ time + group + time:group
#'   Family: gaussian (identity link)
#'   Working correlation: exchangeable
#'   Variance estimator: robust (sandwich) -- default in geepack
#'
#'   For the Gaussian identity-link case, marginal and conditional treatment
#'   effects coincide, so the GEE and LME estimators target the SAME truth
#'   (beta3). Both are compared against the same truth in 04_metrics.R.
#'
#'   Data are sorted by id then time before fitting because geepack requires
#'   observations within each cluster to be contiguous and time-ordered.
fit_gee_continuous <- function(dt) {
  method <- "GEE_continuous"
  
  dt_fit <- dt[!is.na(y1)]
  setorder(dt_fit, id, time)
  
  fit <- tryCatch(
    geepack::geeglm(
      y1 ~ time + group + time:group,
      family   = gaussian(link = "identity"),
      data     = dt_fit,
      id       = id,
      corstr   = "exchangeable",
      std.err  = "san.se"   # robust sandwich SE (default; explicit for clarity)
    ),
    error = function(e) {
      return(.failed_fit(method, reason = conditionMessage(e)))
    }
  )
  
  if (!inherits(fit, "geeglm")) return(fit)
  
  # CRITICAL FIX: use robust .gee_converged() helper instead of bare == 0
  converged <- .gee_converged(fit)
  
  coefs <- tryCatch(
    .extract_gee_coef(fit, method),
    error = function(e) return(.failed_fit(method, reason = conditionMessage(e)))
  )
  if (is.list(coefs)) return(coefs)
  
  list(
    method        = method,
    estimate      = coefs["estimate"],
    se_robust     = coefs["se_robust"],
    se_naive      = coefs["se_naive"],
    pvalue        = coefs["pvalue"],
    converged     = converged,
    singular      = NA,    # not applicable for GEE
    failed        = FALSE,
    failure_reason = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 8. MODEL M3: GENERALIZED LINEAR MIXED MODEL (GLMM) for binary outcome y2
# ------------------------------------------------------------------------------

#' Fit a GLMM with binomial family and random intercept to binary outcome y2.
#'
#' @param dt data.table; long-format dataset with columns id, time, group,
#'   y2 (possibly with NAs from dropout).
#' @return standardized result list.
#' @details Model specification:
#'   y2 ~ time + group + time:group + (1 | id)
#'   Family: binomial(link = "logit")
#'   Random intercept only (no random slope, per protocol)
#'   Fitted by Laplace approximation (the lme4 default for glmer, and the
#'   accepted standard for binary GLMMs in the methodological literature;
#'   adaptive Gauss-Hermite quadrature with nAGQ > 1 is more accurate but
#'   prohibitively slow for 1000 Monte Carlo replicates).
#'
#'   The GLMM estimates the CONDITIONAL (subject-specific) log-odds ratio
#'   alpha3. This is the correct estimand for this model class and the truth
#'   against which GLMM estimates are evaluated in 04_metrics.R.
fit_glmm <- function(dt) {
  method <- "GLMM"
  
  dt_fit <- dt[!is.na(y2)]
  setorder(dt_fit, id, time)
  
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        y2 ~ time + group + time:group + (1 | id),
        family  = binomial(link = "logit"),
        data    = dt_fit,
        nAGQ    = 1,   # Laplace approximation (nAGQ=1); see details above
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl  = list(maxfun = 2e5))
      ),
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      return(.failed_fit(method, reason = conditionMessage(e)))
    }
  )
  
  if (!inherits(fit, "glmerMod")) return(fit)
  
  converged <- .lme4_converged(fit)
  singular  <- isSingular(fit)
  
  coefs <- tryCatch(
    .extract_lme4_coef(fit, method),
    error = function(e) return(.failed_fit(method, reason = conditionMessage(e)))
  )
  if (is.list(coefs)) return(coefs)
  
  list(
    method        = method,
    estimate      = coefs["estimate"],
    se_robust     = coefs["se"],
    se_naive      = coefs["se"],
    pvalue        = coefs["pvalue"],
    converged     = converged,
    singular      = singular,
    failed        = FALSE,
    failure_reason = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 9. MODEL M4: GEE (BINOMIAL) for binary outcome y2
# ------------------------------------------------------------------------------

#' Fit a GEE model with binomial family to the binary outcome y2.
#'
#' @param dt data.table; long-format dataset with columns id, time, group,
#'   y2 (possibly with NAs from dropout). Must be sorted by id then time.
#' @return standardized result list.
#' @details Model specification:
#'   y2 ~ time + group + time:group
#'   Family: binomial(link = "logit")
#'   Working correlation: exchangeable
#'   Variance estimator: robust (sandwich)
#'
#'   The GEE estimates the MARGINAL (population-averaged) log-odds ratio,
#'   which is attenuated relative to the GLMM's conditional alpha3. The
#'   correct truth for GEE bias and coverage evaluation is alpha3_marginal,
#'   pre-computed in 01_generate_data.R via numerical integration and stored
#'   in the scenario's truth list. 04_metrics.R applies the correct truth
#'   for each estimator automatically.
fit_gee_binary <- function(dt) {
  method <- "GEE_binary"
  
  dt_fit <- dt[!is.na(y2)]
  setorder(dt_fit, id, time)
  
  fit <- tryCatch(
    geepack::geeglm(
      y2 ~ time + group + time:group,
      family  = binomial(link = "logit"),
      data    = dt_fit,
      id      = id,
      corstr  = "exchangeable",
      std.err = "san.se"
    ),
    error = function(e) {
      return(.failed_fit(method, reason = conditionMessage(e)))
    }
  )
  
  if (!inherits(fit, "geeglm")) return(fit)
  
  # CRITICAL FIX: use robust .gee_converged() helper instead of bare == 0
  converged <- .gee_converged(fit)
  
  coefs <- tryCatch(
    .extract_gee_coef(fit, method),
    error = function(e) return(.failed_fit(method, reason = conditionMessage(e)))
  )
  if (is.list(coefs)) return(coefs)
  
  list(
    method        = method,
    estimate      = coefs["estimate"],
    se_robust     = coefs["se_robust"],
    se_naive      = coefs["se_naive"],
    pvalue        = coefs["pvalue"],
    converged     = converged,
    singular      = NA,
    failed        = FALSE,
    failure_reason = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 10. TOP-LEVEL: FIT ALL FOUR MODELS TO ONE REPLICATE
# ------------------------------------------------------------------------------

#' Fit all four models (LME, GEE continuous, GLMM, GEE binary) to a single
#' post-dropout Monte Carlo replicate dataset.
#'
#' @param dropout_obj list; output of apply_monotone_dropout() from
#'   02_apply_dropout.R. Must contain $data (long-format data.table with NAs
#'   applied), $truth (list with beta3, alpha3_conditional, alpha3_marginal),
#'   and $params.
#' @return a list with:
#'   fits  - named list of four standardized result lists (one per model),
#'   truth - the truth list passed through from dropout_obj (for
#'           convenience of 04_metrics.R, which needs both fits and truths),
#'   params - echoed for bookkeeping.
fit_all_models <- function(dropout_obj) {
  dt <- dropout_obj$data
  
  list(
    fits = list(
      LME            = fit_lme(dt),
      GEE_continuous = fit_gee_continuous(dt),
      GLMM           = fit_glmm(dt),
      GEE_binary     = fit_gee_binary(dt)
    ),
    truth  = dropout_obj$truth,
    params = dropout_obj$params
  )
}

# ------------------------------------------------------------------------------
# 11. SMOKE TEST (manual, interactive use only -- not sourced by the pipeline)
# ------------------------------------------------------------------------------
# Uncomment to verify locally:
#
# source("01_generate_data.R")
# source("02_apply_dropout.R")
# source("03_fit_models.R")
#
# calib <- build_dropout_calibration_table(n_calib = 5000)
# slopes <- get_dropout_hazard_slopes()
#
# set.seed(42)
# rep_obj <- generate_replicate_data(n = 200, scenario = "alternative")
# do_obj  <- apply_monotone_dropout(rep_obj, target_rate = 0.20,
#                                    calibration_table = calib, slopes = slopes)
# results <- fit_all_models(do_obj)
#
# lapply(results$fits, function(f) {
#   cat(sprintf("%-18s estimate=%.4f  se=%.4f  p=%.4f  converged=%s  singular=%s\n",
#     f$method, f$estimate, f$se_robust, f$pvalue,
#     f$converged, as.character(f$singular)))
# })