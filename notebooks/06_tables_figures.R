# ==============================================================================
# 06_tables_figures.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Reads the consolidated metrics table produced by 05_simulation_engine.R
#   and generates all publication-quality tables and figures for the
#   manuscript. All output is written to a user-specified output directory
#   as PDF (figures) and LaTeX / CSV (tables).
#
# WHY THIS IMPLEMENTATION IS STATISTICALLY APPROPRIATE
#
#   1. TABLES FOLLOW MORRIS ET AL. (2019) REPORTING CONVENTIONS
#      Every numeric metric is accompanied by its Monte Carlo Standard Error
#      (MCSE) in parentheses, formatted as "value (MCSE)". This is the
#      reporting standard for simulation studies in Statistics in Medicine,
#      Statistical Methods in Medical Research, and similar Q1/Q2 outlets.
#      Omitting MCSEs is a common referee complaint.
#
#   2. SEPARATE PANELS FOR NULL AND ALTERNATIVE SCENARIOS
#      Type I Error (null) and Power (alternative) are never placed in the
#      same table column since they answer qualitatively different questions.
#      Each major table is stratified by scenario, then n, then dropout rate.
#
#   3. FIGURES USE A COLOR-BLIND-SAFE PALETTE (Okabe-Ito)
#      The Okabe-Ito palette (Okabe & Ito 2008) is distinguishable both by
#      color-normal viewers and by people with the most common forms of color
#      vision deficiency. Method is also encoded by line type and point shape
#      for black-and-white print compatibility.
#
#   4. ESTIMAND LABELS ARE EXPLICIT IN ALL OUTPUTS
#      Every table and figure that includes binary-outcome results carries an
#      explicit footnote distinguishing "GLMM (conditional alpha3)" from
#      "GEE (Logistic) (marginal alpha3)".
#
#   5. METHOD NAMING STRATEGY (Option A)
#      The `method` column in the raw metrics table uses internal keys
#      ("LME", "GEE_continuous", "GLMM", "GEE_binary"). These are preserved
#      throughout all data manipulation as `method_original` and as the
#      factor LEVELS of `method`. Human-readable display labels are applied
#      ONLY at the ggplot2 scale layer via `labels = METHOD_LABELS`, and in
#      table formatting via direct lookup. This design means:
#        a) All subsetting and grouping logic (e.g., outcome assignment) uses
#           original keys, never display labels -- no silent mismatch.
#        b) OKI_PALETTE, METHOD_LINETYPES, METHOD_SHAPES are all keyed by
#           original names, keeping aesthetic dictionaries consistent.
#        c) Adding a new method in the future requires updating only
#           METHOD_LABELS and the aesthetic dictionaries, with no changes
#           to data manipulation code.
#
#   6. ggplot2 IS THE SOLE PLOTTING DEPENDENCY
#      Multi-panel layouts use facet_grid() / facet_wrap() exclusively.
#      No cowplot, patchwork, or other extensions are required.
#
#   7. LATEX TABLES USE booktabs FORMATTING VIA xtable
#      Tables are written as standalone .tex snippets for \input{} inclusion
#      in the manuscript source.
#
# DEPENDENCIES
#   data.table, ggplot2, xtable
#
# UPSTREAM DEPENDENCY
#   results/metrics_table.rds -- produced by run_simulation_grid() in
#                                 05_simulation_engine.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(xtable)
})

# ------------------------------------------------------------------------------
# 1. CONSTANTS AND STYLE SETTINGS
# ------------------------------------------------------------------------------

# All aesthetic dictionaries are keyed by ORIGINAL method names.
# Display labels are applied only at the ggplot2 scale layer.

#' Human-readable method labels (original name -> display label).
METHOD_LABELS <- c(
  LME            = "LME",
  GEE_continuous = "GEE (Gaussian)",
  GLMM           = "GLMM",
  GEE_binary     = "GEE (Logistic)"
)

#' Okabe-Ito color-blind-safe palette, keyed by original method name.
OKI_PALETTE <- c(
  LME            = "#E69F00",   # orange
  GEE_continuous = "#56B4E9",   # sky blue
  GLMM           = "#009E73",   # bluish green
  GEE_binary     = "#CC79A7"    # reddish purple
)

#' Line types keyed by original method name.
METHOD_LINETYPES <- c(
  LME            = "solid",
  GEE_continuous = "dashed",
  GLMM           = "dotdash",
  GEE_binary     = "dotted"
)

#' Point shapes keyed by original method name.
METHOD_SHAPES <- c(
  LME            = 16L,
  GEE_continuous = 17L,
  GLMM           = 15L,
  GEE_binary     = 18L
)

#' Methods whose outcome is continuous (y1).
CONTINUOUS_METHODS <- c("LME", "GEE_continuous")

#' Methods whose outcome is binary (y2).
BINARY_METHODS <- c("GLMM", "GEE_binary")

#' Nominal significance level (must match ALPHA_LEVEL in 04_metrics.R).
ALPHA_LEVEL <- 0.05

#' Figure dimensions (inches) for a two-column journal layout.
FIG_WIDTH  <- 7.0
FIG_HEIGHT <- 5.5

#' ggplot2 base theme for all figures.
.sim_theme <- function() {
  theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey92", color = NA),
      strip.text        = element_text(size = 9, face = "bold"),
      legend.position   = "bottom",
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = 9),
      axis.text         = element_text(size = 8),
      plot.title        = element_text(size = 11, face = "bold"),
      plot.subtitle     = element_text(size = 8, color = "grey40"),
      plot.caption      = element_text(size = 7, color = "grey50", hjust = 0)
    )
}

# ------------------------------------------------------------------------------
# 2. DATA PREPARATION
# ------------------------------------------------------------------------------

#' Load and prepare the metrics table for plotting and tabulation.
#'
#' @param metrics_path character; path to metrics_table.rds.
#' @return a data.table ready for all downstream functions.
#' @details Key design decisions (see file header, Section 5):
#'   - `method_original` is created first, before any transformation, and
#'     is used for ALL grouping/subsetting logic throughout this file.
#'   - `method` is converted to an ordered factor with original names as
#'     LEVELS (not labels). Display labels are applied only in ggplot scales.
#'   - `outcome` is assigned from `method_original`, never from `method`,
#'     to avoid any risk of label-vs-key mismatch.
prepare_metrics <- function(
    metrics_path = file.path("results", "metrics_table.rds")) {
  
  if (!file.exists(metrics_path)) {
    stop(sprintf(
      "Metrics file not found: %s\nRun 05_simulation_engine.R first.",
      metrics_path
    ))
  }
  
  dt <- readRDS(metrics_path)
  
  # -- Step 1: preserve the original key BEFORE any transformation ----------
  dt[, method_original := as.character(method)]
  
  # Validate: all method_original values must be known keys
  unknown <- setdiff(dt$method_original, names(METHOD_LABELS))
  if (length(unknown) > 0) {
    stop(sprintf(
      "Unknown method keys in metrics table: %s\nUpdate METHOD_LABELS.",
      paste(unknown, collapse = ", ")
    ))
  }
  
  # -- Step 2: convert `method` to ordered factor (levels = original names) -
  # Levels are original names; display labels are applied only in ggplot
  # scale functions via `labels = METHOD_LABELS`. This keeps data
  # manipulation code free of display-label strings.
  dt[, method := factor(method_original,
                        levels = names(METHOD_LABELS))]
  
  # -- Step 3: outcome grouping uses method_original (never method) ---------
  dt[, outcome := fifelse(
    method_original %in% CONTINUOUS_METHODS,
    "Continuous outcome (y1)",
    "Binary outcome (y2)"
  )]
  
  # -- Step 4: convenience display columns ----------------------------------
  dt[, dropout_pct := factor(
    sprintf("%d%%", as.integer(target_dropout_rate * 100)),
    levels = c("10%", "20%", "30%", "40%")
  )]
  
  dt[, n_label := factor(
    sprintf("N = %d", n),
    levels = sprintf("N = %d", c(100L, 200L, 500L))
  )]
  
  dt[, scenario_label := fifelse(
    scenario == "null",
    "Null scenario (Type I Error)",
    "Alternative scenario"
  )]
  
  dt[]
}

# ------------------------------------------------------------------------------
# 3. SHARED SCALE CONSTRUCTORS
# ------------------------------------------------------------------------------
# These return consistent ggplot2 scale layers for color, linetype, and shape,
# all keyed by original method names and displaying METHOD_LABELS.

#' ggplot2 color scale for method, with Okabe-Ito palette.
.scale_color_method <- function() {
  scale_color_manual(
    name   = "Method",
    values = OKI_PALETTE,
    labels = METHOD_LABELS
  )
}

#' ggplot2 linetype scale for method.
.scale_linetype_method <- function() {
  scale_linetype_manual(
    name   = "Method",
    values = METHOD_LINETYPES,
    labels = METHOD_LABELS
  )
}

#' ggplot2 shape scale for method.
.scale_shape_method <- function() {
  scale_shape_manual(
    name   = "Method",
    values = METHOD_SHAPES,
    labels = METHOD_LABELS
  )
}

# ------------------------------------------------------------------------------
# 4. FORMATTING HELPERS
# ------------------------------------------------------------------------------

#' Format a metric value and its MCSE as "value (MCSE)".
#'
#' @param x numeric vector of metric values.
#' @param mcse numeric vector of Monte Carlo SEs.
#' @param digits integer; decimal places for the value.
#' @param mcse_digits integer; decimal places for the MCSE.
#' @return character vector.
.fmt_mcse <- function(x, mcse, digits = 3, mcse_digits = 3) {
  val_str  <- formatC(x,    digits = digits,      format = "f", flag = " ")
  mcse_str <- formatC(mcse, digits = mcse_digits, format = "f", flag = " ")
  ifelse(is.na(x), "---",
         sprintf("%s (%s)", trimws(val_str), trimws(mcse_str)))
}

#' Format a proportion in [0,1] as a percentage with its MCSE.
#'
#' @param x numeric vector in [0, 1].
#' @param mcse numeric vector.
#' @return character vector, e.g. "94.8 (0.7)".
.fmt_pct_mcse <- function(x, mcse) {
  val_str  <- formatC(x * 100,    digits = 1, format = "f")
  mcse_str <- formatC(mcse * 100, digits = 1, format = "f")
  ifelse(is.na(x), "---",
         sprintf("%s (%s)", trimws(val_str), trimws(mcse_str)))
}

#' Save a data.table as both CSV and a booktabs LaTeX table.
#'
#' @param wide data.table; formatted wide table ready for output.
#' @param output_dir character; output directory.
#' @param base_name character; file name without extension.
#' @param caption character; LaTeX table caption.
#' @param label character; LaTeX \label{} string.
.save_table <- function(wide, output_dir, base_name, caption, label) {
  fwrite(wide, file = file.path(output_dir, paste0(base_name, ".csv")))
  
  xt <- xtable(wide, caption = caption, label = label)
  print(
    xt,
    file              = file.path(output_dir, paste0(base_name, ".tex")),
    include.rownames  = FALSE,
    booktabs          = TRUE,
    caption.placement = "top",
    sanitize.text.function = identity
  )
  invisible(wide)
}

# ------------------------------------------------------------------------------
# 5. TABLE 1: BIAS AND RMSE (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Build and save the Bias / RMSE table for the alternative scenario.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_table_bias_rmse <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative",
            .(method, n_label, dropout_pct,
              bias, mcse_bias, rmse, mcse_rmse, relative_bias)]
  
  sub[, bias_fmt     := .fmt_mcse(bias, mcse_bias, digits = 3)]
  sub[, rmse_fmt     := .fmt_mcse(rmse, mcse_rmse, digits = 3)]
  sub[, rel_bias_fmt := ifelse(is.na(relative_bias), "---",
                               formatC(relative_bias * 100, digits = 1, format = "f"))]
  
  # Use METHOD_LABELS for column headers in the wide table
  sub[, method_label := METHOD_LABELS[as.character(method)]]
  
  wide_bias <- dcast(sub, n_label + dropout_pct ~ method_label,
                     value.var = "bias_fmt")
  wide_rmse <- dcast(sub, n_label + dropout_pct ~ method_label,
                     value.var = "rmse_fmt")
  
  # Stack: bias rows then rmse rows with a panel identifier
  wide_bias[, Metric := "Bias (MCSE)"]
  wide_rmse[, Metric := "RMSE (MCSE)"]
  wide <- rbindlist(list(wide_bias, wide_rmse))
  setcolorder(wide, c("Metric", "n_label", "dropout_pct"))
  
  .save_table(
    wide, output_dir, "table_bias_rmse",
    caption = paste0(
      "Bias (Monte Carlo SE) and RMSE (Monte Carlo SE) of the ",
      "time$\\times$group interaction estimator under the alternative scenario ",
      "($\\beta_3 = 1.0$, continuous outcome; $\\alpha_3 = 0.5$ conditional, ",
      "binary outcome). GEE (Logistic) is evaluated against the marginal ",
      "(population-averaged) $\\alpha_3$; GLMM against the conditional ",
      "(subject-specific) $\\alpha_3$. Each cell: value (MCSE) over 1000 replicates."
    ),
    label = "tab:bias_rmse"
  )
}

# ------------------------------------------------------------------------------
# 6. TABLE 2: COVERAGE AND SE RATIO (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Build and save the Coverage / SE Ratio table for the alternative scenario.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_table_coverage_se <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative",
            .(method, n_label, dropout_pct,
              coverage, mcse_coverage, se_ratio, se_ratio_naive)]
  
  sub[, coverage_fmt       := .fmt_pct_mcse(coverage, mcse_coverage)]
  sub[, se_ratio_fmt       := ifelse(is.na(se_ratio),
                                     "---", formatC(se_ratio,       digits = 3, format = "f"))]
  sub[, se_ratio_naive_fmt := ifelse(is.na(se_ratio_naive),
                                     "---", formatC(se_ratio_naive, digits = 3, format = "f"))]
  
  sub[, method_label := METHOD_LABELS[as.character(method)]]
  
  wide_cov  <- dcast(sub, n_label + dropout_pct ~ method_label, value.var = "coverage_fmt")
  wide_ser  <- dcast(sub, n_label + dropout_pct ~ method_label, value.var = "se_ratio_fmt")
  wide_sern <- dcast(sub, n_label + dropout_pct ~ method_label, value.var = "se_ratio_naive_fmt")
  
  wide_cov[,  Metric := "Coverage 95\\% CI (MCSE)"]
  wide_ser[,  Metric := "SE Ratio (robust)"]
  wide_sern[, Metric := "SE Ratio (naive)"]
  
  wide <- rbindlist(list(wide_cov, wide_ser, wide_sern))
  setcolorder(wide, c("Metric", "n_label", "dropout_pct"))
  
  .save_table(
    wide, output_dir, "table_coverage_se",
    caption = paste0(
      "Coverage probability (\\%) of nominal 95\\% Wald confidence intervals ",
      "and SE Ratio (average reported SE / empirical SE) under the alternative ",
      "scenario. SE Ratio (robust): uses the sandwich SE for GEE methods and ",
      "model-based SE for LME/GLMM. SE Ratio (naive): model-based SE for all ",
      "methods (diagnostic for GEE). SE Ratio $= 1$ indicates perfect calibration; ",
      "$< 1$ anti-conservative. Monte Carlo SE in parentheses for coverage."
    ),
    label = "tab:coverage_se"
  )
}

# ------------------------------------------------------------------------------
# 7. TABLE 3: TYPE I ERROR AND POWER
# ------------------------------------------------------------------------------

#' Build and save the Type I Error / Power table (both scenarios).
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_table_rejection_rates <- function(dt, output_dir) {
  
  sub <- dt[, .(method, scenario, n_label, dropout_pct,
                rejection_rate, mcse_rejection)]
  sub[, rej_fmt      := .fmt_pct_mcse(rejection_rate, mcse_rejection)]
  sub[, method_label := METHOD_LABELS[as.character(method)]]
  
  wide_null <- dcast(sub[scenario == "null"],
                     n_label + dropout_pct ~ method_label, value.var = "rej_fmt")
  wide_alt  <- dcast(sub[scenario == "alternative"],
                     n_label + dropout_pct ~ method_label, value.var = "rej_fmt")
  
  wide_null[, Panel := "Type I Error -- Null ($\\beta_3 = \\alpha_3 = 0$)"]
  wide_alt[,  Panel := "Power -- Alternative ($\\beta_3 = 1.0$, $\\alpha_3 = 0.5$ conditional)"]
  
  wide <- rbindlist(list(wide_null, wide_alt))
  setcolorder(wide, c("Panel", "n_label", "dropout_pct"))
  
  .save_table(
    wide, output_dir, "table_rejection_rates",
    caption = paste0(
      "Rejection rates (\\%) at the nominal 5\\% significance level, with ",
      "Monte Carlo SE in parentheses. Upper panel: Type I Error (null scenario). ",
      "Lower panel: Power (alternative scenario). ",
      "Nominal Type I Error = 5\\%; conventional power threshold = 80\\%."
    ),
    label = "tab:rejection"
  )
}

# ------------------------------------------------------------------------------
# 8. TABLE 4: CONVERGENCE AND SINGULARITY RATES
# ------------------------------------------------------------------------------

#' Build and save the convergence / singularity rate diagnostic table.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_table_convergence <- function(dt, output_dir) {
  
  sub <- dt[, .(method, method_original, scenario, n_label, dropout_pct,
                convergence_rate, singularity_rate)]
  
  sub[, conv_fmt := formatC(convergence_rate * 100, digits = 1, format = "f")]
  sub[, sing_fmt := ifelse(
    is.na(singularity_rate), "N/A",
    formatC(singularity_rate * 100, digits = 1, format = "f")
  )]
  sub[, method_label := METHOD_LABELS[method_original]]
  
  wide_conv <- dcast(sub, scenario + n_label + dropout_pct ~ method_label,
                     value.var = "conv_fmt")
  wide_sing <- dcast(sub, scenario + n_label + dropout_pct ~ method_label,
                     value.var = "sing_fmt")
  
  wide_conv[, Metric := "Convergence rate (\\%)"]
  wide_sing[, Metric := "Singularity rate (\\%)"]
  
  wide <- rbindlist(list(wide_conv, wide_sing))
  setcolorder(wide, c("Metric", "scenario", "n_label", "dropout_pct"))
  
  .save_table(
    wide, output_dir, "table_convergence",
    caption = paste0(
      "Convergence rate (\\%) and singularity rate (\\%) across 1000 Monte Carlo ",
      "replicates. Singularity rate is reported only for random-effects models ",
      "(LME, GLMM); N/A indicates the concept is not applicable (GEE)."
    ),
    label = "tab:convergence"
  )
}

# ------------------------------------------------------------------------------
# 9. FIGURE 1: BIAS (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Plot absolute bias vs. dropout rate, faceted by outcome and n.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_figure_bias <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative"]
  
  p <- ggplot(sub, aes(x = dropout_pct, y = bias,
                       color = method, linetype = method,
                       shape = method, group = method)) +
    geom_hline(yintercept = 0, linewidth = 0.4,
               color = "grey55", linetype = "solid") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_errorbar(
      aes(ymin = bias - 1.96 * mcse_bias,
          ymax = bias + 1.96 * mcse_bias),
      width = 0.15, linewidth = 0.4
    ) +
    facet_grid(outcome ~ n_label, scales = "free_y") +
    .scale_color_method() +
    .scale_linetype_method() +
    .scale_shape_method() +
    labs(
      title    = "Bias of the time \u00d7 group interaction estimator",
      subtitle = "Alternative scenario. Error bars: \u00b1 1.96 \u00d7 MCSE.",
      x        = "Target dropout rate",
      y        = "Bias",
      caption  = paste0(
        "GLMM is evaluated against the conditional (subject-specific) treatment effect. ",
        "GEE (Logistic) is evaluated against the marginal (population-averaged) treatment effect."
      )
    ) +
    .sim_theme()
  
  ggsave(file.path(output_dir, "figure_bias.pdf"),
         plot = p, width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, device = "pdf")
  invisible(p)
}

# ------------------------------------------------------------------------------
# 10. FIGURE 2: COVERAGE (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Plot 95% CI coverage vs. dropout rate, faceted by outcome and n.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_figure_coverage <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative"]
  
  p <- ggplot(sub, aes(x = dropout_pct, y = coverage * 100,
                       color = method, linetype = method,
                       shape = method, group = method)) +
    geom_hline(yintercept = 95, linewidth = 0.4,
               color = "grey55", linetype = "dashed") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_errorbar(
      aes(ymin = (coverage - 1.96 * mcse_coverage) * 100,
          ymax = (coverage + 1.96 * mcse_coverage) * 100),
      width = 0.15, linewidth = 0.4
    ) +
    facet_grid(outcome ~ n_label) +
    scale_y_continuous(limits = c(NA, 100), breaks = seq(80, 100, by = 5)) +
    .scale_color_method() +
    .scale_linetype_method() +
    .scale_shape_method() +
    labs(
      title    = "Coverage probability of nominal 95% Wald confidence intervals",
      subtitle = "Alternative scenario. Dashed line: nominal 95% level. Error bars: \u00b1 1.96 \u00d7 MCSE.",
      x        = "Target dropout rate",
      y        = "Coverage (%)",
      caption  = paste0(
        "GLMM evaluated against conditional \u03b1\u2083; ",
        "GEE (Logistic) evaluated against marginal \u03b1\u2083."
      )
    ) +
    .sim_theme()
  
  ggsave(file.path(output_dir, "figure_coverage.pdf"),
         plot = p, width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, device = "pdf")
  invisible(p)
}

# ------------------------------------------------------------------------------
# 11. FIGURE 3: RMSE (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Plot RMSE vs. dropout rate, faceted by outcome and n.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_figure_rmse <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative"]
  
  p <- ggplot(sub, aes(x = dropout_pct, y = rmse,
                       color = method, linetype = method,
                       shape = method, group = method)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_errorbar(
      aes(ymin = rmse - 1.96 * mcse_rmse,
          ymax = rmse + 1.96 * mcse_rmse),
      width = 0.15, linewidth = 0.4
    ) +
    facet_grid(outcome ~ n_label, scales = "free_y") +
    .scale_color_method() +
    .scale_linetype_method() +
    .scale_shape_method() +
    labs(
      title    = "Root Mean Squared Error (RMSE)",
      subtitle = "Alternative scenario. Error bars: \u00b1 1.96 \u00d7 MCSE.",
      x        = "Target dropout rate",
      y        = "RMSE",
      caption  = paste0(
        "GLMM evaluated against conditional \u03b1\u2083; ",
        "GEE (Logistic) evaluated against marginal \u03b1\u2083."
      )
    ) +
    .sim_theme()
  
  ggsave(file.path(output_dir, "figure_rmse.pdf"),
         plot = p, width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, device = "pdf")
  invisible(p)
}

# ------------------------------------------------------------------------------
# 12. FIGURE 4: REJECTION RATES (TYPE I ERROR AND POWER)
# ------------------------------------------------------------------------------

#' Plot rejection rates for both scenarios as a dual-panel figure.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_figure_rejection <- function(dt, output_dir) {
  
  sub <- dt[, .(method, scenario_label, n_label, dropout_pct,
                rejection_rate, mcse_rejection)]
  
  ref_lines <- data.table(
    scenario_label = c("Null scenario (Type I Error)", "Alternative scenario"),
    yref           = c(ALPHA_LEVEL * 100, 80)
  )
  
  p <- ggplot(sub, aes(x = dropout_pct, y = rejection_rate * 100,
                       color = method, linetype = method,
                       shape = method, group = method)) +
    geom_hline(
      data = ref_lines,
      aes(yintercept = yref),
      linewidth = 0.4, color = "grey55", linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    geom_errorbar(
      aes(ymin = (rejection_rate - 1.96 * mcse_rejection) * 100,
          ymax = (rejection_rate + 1.96 * mcse_rejection) * 100),
      width = 0.15, linewidth = 0.4
    ) +
    facet_grid(scenario_label ~ n_label, scales = "free_y") +
    .scale_color_method() +
    .scale_linetype_method() +
    .scale_shape_method() +
    labs(
      title    = "Rejection rates at the 5% nominal significance level",
      subtitle = paste0(
        "Upper panels: Type I Error (null, \u03b2\u2083 = \u03b1\u2083 = 0). ",
        "Lower panels: Power (alternative). ",
        "Dashed lines: 5% and 80% reference levels."
      ),
      x        = "Target dropout rate",
      y        = "Rejection rate (%)",
      caption  = "Error bars: \u00b1 1.96 \u00d7 MCSE."
    ) +
    .sim_theme()
  
  ggsave(file.path(output_dir, "figure_rejection_rates.pdf"),
         plot = p, width = FIG_WIDTH, height = FIG_HEIGHT + 2, device = "pdf")
  invisible(p)
}

# ------------------------------------------------------------------------------
# 13. FIGURE 5: SE RATIO (ALTERNATIVE SCENARIO)
# ------------------------------------------------------------------------------

#' Plot SE Ratio (avg robust SE / empirical SE) vs. dropout rate.
#'
#' @param dt data.table; output of prepare_metrics().
#' @param output_dir character.
make_figure_se_ratio <- function(dt, output_dir) {
  
  sub <- dt[scenario == "alternative",
            .(method, outcome, n_label, dropout_pct, se_ratio)]
  
  p <- ggplot(sub, aes(x = dropout_pct, y = se_ratio,
                       color = method, linetype = method,
                       shape = method, group = method)) +
    geom_hline(yintercept = 1, linewidth = 0.4,
               color = "grey55", linetype = "dashed") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    facet_grid(outcome ~ n_label) +
    .scale_color_method() +
    .scale_linetype_method() +
    .scale_shape_method() +
    labs(
      title    = "SE Ratio: average reported SE / empirical SE",
      subtitle = paste0(
        "Alternative scenario. Dashed line: perfect calibration (ratio = 1). ",
        "Values < 1 indicate anti-conservative standard error estimation."
      ),
      x        = "Target dropout rate",
      y        = "SE Ratio",
      caption  = "For GEE methods, SE Ratio uses the robust (sandwich) SE."
    ) +
    .sim_theme()
  
  ggsave(file.path(output_dir, "figure_se_ratio.pdf"),
         plot = p, width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, device = "pdf")
  invisible(p)
}

# ------------------------------------------------------------------------------
# 14. TOP-LEVEL ENTRY POINT
# ------------------------------------------------------------------------------

#' Generate all tables and figures for the manuscript.
#'
#' @param metrics_path character; path to metrics_table.rds.
#' @param output_dir character; directory for all PDF and .tex outputs.
#' @return invisible NULL; all outputs written to disk as side effects.
generate_all_outputs <- function(
    metrics_path = file.path("results", "metrics_table.rds"),
    output_dir   = file.path("results", "tables_figures")) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  message("Loading and preparing metrics table ...")
  dt <- prepare_metrics(metrics_path)
  
  message("Building tables ...")
  make_table_bias_rmse(dt, output_dir)
  make_table_coverage_se(dt, output_dir)
  make_table_rejection_rates(dt, output_dir)
  make_table_convergence(dt, output_dir)
  
  message("Building figures ...")
  make_figure_bias(dt, output_dir)
  make_figure_coverage(dt, output_dir)
  make_figure_rmse(dt, output_dir)
  make_figure_rejection(dt, output_dir)
  make_figure_se_ratio(dt, output_dir)
  
  message(sprintf("Done. All outputs written to: %s", normalizePath(output_dir)))
  invisible(NULL)
}