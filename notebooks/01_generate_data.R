# ==============================================================================
# 01_generate_data.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Defines the complete-data generating mechanism (DGM) for the simulation
#   study comparing LME, GEE, and GLMM under monotone dropout. This file
#   produces FULLY OBSERVED (pre-dropout) longitudinal data for both a
#   continuous outcome (y1) and a binary outcome (y2). Missingness is applied
#   downstream in 02_apply_dropout.R.
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#   1. Continuous outcome (y1): generated from a linear mixed model with a
#      JOINTLY sampled random intercept/slope pair, b_i = (b0i, b1i), drawn
#      from a bivariate normal distribution with Corr(b0,b1) = 0.30. Joint
#      sampling (via Cholesky factorization of the 2x2 covariance matrix) is
#      essential -- generating b0i and b1i independently and only inducing
#      correlation afterward is a common and statistically incorrect shortcut.
#      Because the link is identity, the marginal mean E[Y1 | time, group] is
#      IDENTICAL to the conditional mean -- so beta3 is a valid target of
#      inference for BOTH LME (conditional) and GEE (marginal) estimators.
#      No truth-adjustment is required for y1.
#
#   2. Binary outcome (y2): generated from a logistic mixed model with its
#      OWN random intercept u_i ~ N(0, sigma_u^2), independent of the
#      continuous-outcome random effects (b0i, b1i). This is a deliberate
#      design choice, not an oversight: the study evaluates LME/GEE on y1
#      and GLMM/GEE on y2 as two SEPARATE estimation problems, each with its
#      own treatment-effect target (beta3 for y1, alpha3 for y2). Critically,
#      the conditional (subject-specific) log-odds ratio alpha3 is NOT the
#      same quantity that a marginal (population-averaged) GEE model
#      estimates, because E_u[expit(eta)] != expit(E_u[eta]) for any
#      nonlinear link (Neuhaus, Kalbfleisch & Hauck 1991; Heagerty & Zeger
#      2000). Reporting GEE bias against the conditional alpha3 would be a
#      methodological error -- GEE would appear "biased" purely as an
#      artifact of comparing it to the wrong estimand. We therefore compute,
#      in closed numerical form, the TRUE marginal interaction coefficient
#      that a correctly specified large-sample GEE fit converges to. This is
#      done by:
#        a) integrating the conditional probability P(Y=1 | time, group, u)
#           over the random-intercept distribution at each design point
#           (time x group cell) via numerical quadrature (stats::integrate),
#           producing the exact marginal probability surface;
#        b) fitting a saturated logistic regression to that exact marginal
#           probability surface (treated as an infinite-sample population),
#           and extracting the time:group coefficient. This is the
#           "least-false parameter" the GEE estimator targets asymptotically
#           under the working logistic mean model (Zeger, Liang & Albert,
#           1988, Biometrics).
#      This marginal truth is a FIXED property of the DGM (computed once per
#      scenario, not per Monte Carlo draw) and is attached to the scenario
#      object for use in 04_metrics.R.
#
#   3. All random-effects and residual draws use VECTORIZED base-R RNG calls
#      operating on the full subject (or subject-visit) index. No explicit
#      for-loops over subjects/visits are used for outcome generation -- this
#      keeps the per-replicate RNG stream short, reproducible, and fast
#      enough for 1000 Monte Carlo replicates per design cell.
#
#   4. Random effects are generated via base-R chol() rather than MASS::mvrnorm().
#      This avoids an unnecessary dependency (MASS) for a one-line operation,
#      and gives identical correlated draws given the same underlying
#      standard-normal stream and seed.
#
# DEPENDENCIES: data.table (storage/manipulation only; no randomization here)
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 1. FIXED STUDY DESIGN CONSTANTS
# ------------------------------------------------------------------------------

#' Visit/time grid for the trial (Visit 0 = baseline, ..., Visit 3)
SIM_TIME_GRID <- 0:3

#' Group codes (0 = control, 1 = active treatment)
SIM_GROUP_LEVELS <- c(0L, 1L)

# ------------------------------------------------------------------------------
# 2. DATA-GENERATING PARAMETERS
# ------------------------------------------------------------------------------
# beta2 / alpha2 (the group MAIN effect) are fixed at 0 because randomization
# guarantees no baseline group difference in expectation; the entire
# treatment signal of interest is carried by the time x group interaction.

#' Build the full parameter list for one effect-size scenario.
#'
#' @param scenario character; one of "null" or "alternative".
#' @return a named list of all DGM parameters, including the scenario label,
#'   the conditional treatment effects (beta3, alpha3), and variance
#'   components for both outcomes.
#' @details Under "null", beta3 = alpha3 = 0 (used to estimate Type I error).
#'   Under "alternative", beta3 and alpha3 take fixed non-null values (used
#'   to estimate power, bias, RMSE, and coverage at a realistic effect size).
get_scenario_parameters <- function(scenario = c("null", "alternative")) {
  scenario <- match.arg(scenario)
  
  beta3_value  <- if (scenario == "null") 0 else 1.0   # continuous: slope diff
  alpha3_value <- if (scenario == "null") 0 else 0.5   # binary: log-OR diff
  
  list(
    scenario = scenario,
    
    # ---- Continuous outcome (y1): linear mixed model -----------------------
    beta0       = 50.0,        # baseline mean (e.g., a symptom severity score)
    beta1       = -1.0,        # time (control-group) slope: mean decline/visit
    beta2       = 0.0,         # group main effect at baseline (~0 by randomization)
    beta3       = beta3_value, # TREATMENT EFFECT OF INTEREST: time x group
    sigma_b0    = 5.0,         # SD of random intercept
    sigma_b1    = 1.5,         # SD of random slope
    rho_b       = 0.30,        # Corr(b0i, b1i) -- intentional, per protocol
    sigma_resid = 3.0,         # residual (within-subject) SD
    
    # ---- Binary outcome (y2): logistic mixed model --------------------------
    alpha0  = qlogis(0.30),    # baseline log-odds, ~30% baseline event rate
    alpha1  = -0.20,           # time (control-group) log-odds slope
    alpha2  = 0.0,             # group main effect at baseline (~0 by randomization)
    alpha3  = alpha3_value,    # TREATMENT EFFECT OF INTEREST: time x group
    sigma_u = 1.0              # SD of random intercept (logistic scale)
  )
}

# ------------------------------------------------------------------------------
# 3. DESIGN SKELETON (subject x visit long-format frame)
# ------------------------------------------------------------------------------

#' Construct the long-format subject-by-visit design skeleton.
#'
#' @param n integer; number of subjects (assumed even for exact 1:1 allocation).
#' @return a data.table with one row per subject-visit: id, time, group.
#' @details Treatment assignment is permuted-block randomized 1:1 at the
#'   subject level (not re-randomized per visit). Using a balanced allocation
#'   (rather than independent Bernoulli draws per subject) avoids unintended
#'   sample-size imbalance at small N, which would otherwise inflate Monte
#'   Carlo variance for reasons unrelated to the methods being compared.
build_design_skeleton <- function(n) {
  if (n %% 2 != 0) {
    warning("n is odd; 1:1 allocation will be off by one subject.")
  }
  
  n_treat   <- floor(n / 2)
  n_control <- n - n_treat
  group_assignment <- sample(c(rep(1L, n_treat), rep(0L, n_control)), size = n)
  
  n_visits <- length(SIM_TIME_GRID)
  
  data.table(
    id    = rep(seq_len(n), each = n_visits),
    time  = rep(SIM_TIME_GRID, times = n),
    group = rep(group_assignment, each = n_visits)
  )
}

# ------------------------------------------------------------------------------
# 4. CONTINUOUS OUTCOME (y1): joint random intercept/slope + LME
# ------------------------------------------------------------------------------

#' Jointly sample correlated random intercepts and slopes.
#'
#' @param n integer; number of subjects.
#' @param sigma_b0,sigma_b1 numeric; random-effect standard deviations.
#' @param rho_b numeric; correlation between intercept and slope.
#' @return an n x 2 matrix with columns b0, b1.
#' @details Random effects are drawn via Cholesky factorization of the 2x2
#'   covariance matrix applied to independent standard normal draws,
#'   guaranteeing the specified correlation exactly rather than
#'   approximately. This is the correct way to induce dependence between
#'   random effects; sampling them independently and "correlating after the
#'   fact" would not reproduce a valid bivariate normal joint distribution.
sample_joint_random_effects <- function(n, sigma_b0, sigma_b1, rho_b) {
  if (sigma_b0 <= 0 || sigma_b1 <= 0) {
    stop("sigma_b0 and sigma_b1 must be strictly positive.")
  }
  if (abs(rho_b) >= 1) {
    stop("rho_b must lie strictly within (-1, 1).")
  }
  
  Sigma_b <- matrix(
    c(sigma_b0^2,                  rho_b * sigma_b0 * sigma_b1,
      rho_b * sigma_b0 * sigma_b1, sigma_b1^2),
    nrow = 2, ncol = 2
  )
  
  L <- chol(Sigma_b)          # upper-triangular: t(L) %*% L = Sigma_b
  z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  b <- z %*% L                # rows ~ MVN(0, Sigma_b)
  
  colnames(b) <- c("b0", "b1")
  b
}

#' Simulate the continuous outcome y1 onto an existing design skeleton.
#'
#' @param dt data.table produced by build_design_skeleton().
#' @param params list returned by get_scenario_parameters().
#' @return the input data.table augmented with b0, b1 (random effects, kept
#'   for diagnostic/debugging purposes only -- not used as predictors by the
#'   fitted models) and y1 (the observed continuous outcome).
simulate_continuous_outcome <- function(dt, params) {
  dt <- copy(dt)
  ids <- unique(dt$id)
  n <- length(ids)
  
  re <- sample_joint_random_effects(n, params$sigma_b0, params$sigma_b1, params$rho_b)
  re_dt <- data.table(id = ids, b0 = re[, "b0"], b1 = re[, "b1"])
  
  dt <- merge(dt, re_dt, by = "id", sort = FALSE)
  setorder(dt, id, time)
  
  eps <- rnorm(nrow(dt), mean = 0, sd = params$sigma_resid)
  
  dt[, y1 := params$beta0 +
       params$beta1 * time +
       params$beta2 * group +
       params$beta3 * time * group +
       b0 + b1 * time + eps]
  
  dt[]
}

# ------------------------------------------------------------------------------
# 5. BINARY OUTCOME (y2): random intercept + GLMM
# ------------------------------------------------------------------------------

#' Simulate the binary outcome y2 onto an existing (continuous-augmented)
#' design skeleton.
#'
#' @param dt data.table that already contains id, time, group columns.
#' @param params list returned by get_scenario_parameters().
#' @return the input data.table augmented with u (random intercept, for
#'   diagnostics only), the true conditional probability p2, and the
#'   realized binary outcome y2.
simulate_binary_outcome <- function(dt, params) {
  dt <- copy(dt)
  ids <- unique(dt$id)
  n <- length(ids)
  
  u <- rnorm(n, mean = 0, sd = params$sigma_u)
  u_dt <- data.table(id = ids, u = u)
  
  dt <- merge(dt, u_dt, by = "id", sort = FALSE)
  setorder(dt, id, time)
  
  eta <- params$alpha0 +
    params$alpha1 * dt$time +
    params$alpha2 * dt$group +
    params$alpha3 * dt$time * dt$group +
    dt$u
  
  dt[, p2 := plogis(eta)]
  dt[, y2 := rbinom(.N, size = 1, prob = p2)]
  
  dt[]
}

# ------------------------------------------------------------------------------
# 6. MARGINAL (POPULATION-AVERAGED) TRUTH FOR alpha3 -- GEE TARGET ESTIMAND
# ------------------------------------------------------------------------------

#' Compute the exact marginal probability P(Y2 = 1 | time, group) implied by
#' the conditional logistic-normal model, by numerically integrating over the
#' random-intercept distribution.
#'
#' @param eta_fixed numeric scalar; the fixed-effects linear predictor
#'   alpha0 + alpha1*time + alpha2*group + alpha3*time*group at one design point.
#' @param sigma_u numeric; SD of the random intercept.
#' @return numeric scalar; the marginal probability at that design point.
.marginal_probability <- function(eta_fixed, sigma_u) {
  integrand <- function(u) plogis(eta_fixed + u) * dnorm(u, mean = 0, sd = sigma_u)
  # Integrate over a wide range; sigma_u is modest so +-8 SD is ample.
  result <- integrate(integrand, lower = -8 * sigma_u, upper = 8 * sigma_u,
                      rel.tol = 1e-10)
  result$value
}

#' Compute the true marginal (population-averaged) time x group interaction
#' coefficient -- the estimand that a correctly specified GEE logistic model
#' targets asymptotically -- for a given scenario's conditional parameters.
#'
#' @param params list returned by get_scenario_parameters().
#' @param time_grid numeric vector of visit times (defaults to SIM_TIME_GRID).
#' @param group_levels numeric vector of group codes (defaults to c(0,1)).
#' @return a list with: marginal_alpha3 (the scalar truth for GEE bias/coverage),
#'   and marginal_prob_table (the underlying marginal probability surface, for
#'   transparency / unit testing).
#' @details See file header (section 2 of the rationale) for full
#'   methodological justification. Briefly: exact marginal probabilities are
#'   computed at every (time, group) design cell via Gaussian quadrature
#'   (stats::integrate), then a saturated logistic regression is fit to
#'   those exact probabilities (weighted as an effectively infinite
#'   pseudo-sample) to recover the time:group log-odds coefficient implied
#'   by the true marginal mean structure. This is the "least-false
#'   parameter" targeted by a large-sample GEE fit under a logistic working
#'   mean model.
compute_marginal_alpha3 <- function(params,
                                    time_grid = SIM_TIME_GRID,
                                    group_levels = SIM_GROUP_LEVELS) {
  
  design <- expand.grid(time = time_grid, group = group_levels)
  setDT(design)
  
  design[, eta_fixed := params$alpha0 +
           params$alpha1 * time +
           params$alpha2 * group +
           params$alpha3 * time * group]
  
  design[, p_marginal := mapply(.marginal_probability,
                                eta_fixed = eta_fixed,
                                MoreArgs = list(sigma_u = params$sigma_u))]
  
  # Fit a saturated logistic regression to the exact marginal probabilities,
  # using a large pseudo-sample size so the fit recovers the population
  # coefficients to numerical precision (this is NOT a Monte Carlo draw --
  # it is a deterministic functional of params, computed once per scenario).
  pseudo_n <- 1e6
  design[, successes := p_marginal * pseudo_n]
  design[, failures  := (1 - p_marginal) * pseudo_n]
  
  fit <- glm(
    cbind(successes, failures) ~ time * group,
    family = binomial(link = "logit"),
    data = design
  )
  
  marginal_alpha3 <- unname(coef(fit)["time:group"])
  
  list(
    marginal_alpha3     = marginal_alpha3,
    marginal_prob_table = design[, .(time, group, p_marginal)]
  )
}

# ------------------------------------------------------------------------------
# 7. TOP-LEVEL GENERATOR: one fully observed Monte Carlo dataset
# ------------------------------------------------------------------------------

#' Generate one complete (pre-dropout) Monte Carlo replicate dataset.
#'
#' @param n integer; number of subjects.
#' @param scenario character; "null" or "alternative" (see get_scenario_parameters()).
#' @return a list with:
#'   data   - the long-format data.table (id, time, group, y1, y2, plus
#'            diagnostic columns b0, b1, u, p2 retained for downstream
#'            dropout-mechanism use and debugging),
#'   params - the scenario's conditional DGM parameters,
#'   truth  - a list of estimand truths: beta3 (conditional, valid for both
#'            LME and GEE on the continuous outcome), alpha3_conditional
#'            (the GLMM target), alpha3_marginal (the GEE target).
#' @details This function does NOT set the RNG seed; seeding is the
#'   responsibility of the calling simulation engine (05_simulation_engine.R),
#'   to guarantee a single, centrally controlled, reproducible RNG stream
#'   across the full N x dropout x scenario x replicate grid.
generate_replicate_data <- function(n, scenario = c("null", "alternative")) {
  scenario <- match.arg(scenario)
  
  if (!is.numeric(n) || n <= 0 || n != round(n)) {
    stop("n must be a positive integer.")
  }
  
  params <- get_scenario_parameters(scenario)
  
  dt <- build_design_skeleton(n)
  dt <- simulate_continuous_outcome(dt, params)
  dt <- simulate_binary_outcome(dt, params)
  
  marginal <- compute_marginal_alpha3(params)
  
  list(
    data = dt,
    params = params,
    truth = list(
      beta3              = params$beta3,
      alpha3_conditional = params$alpha3,
      alpha3_marginal    = marginal$marginal_alpha3
    )
  )
}

# ------------------------------------------------------------------------------
# 8. SMOKE TEST (manual, interactive use only -- not sourced by the pipeline)
# ------------------------------------------------------------------------------
# Uncomment to verify locally:
#
# source("01_generate_data.R")
# set.seed(1)
# out <- generate_replicate_data(n = 200, scenario = "alternative")
# print(out$truth)        # check alpha3_marginal is attenuated vs alpha3_conditional
# print(head(out$data))
# print(out$data[time == 0, .N, by = group])   # check ~1:1 balance