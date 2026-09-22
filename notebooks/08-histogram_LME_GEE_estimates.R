# ============================================================
# Histogram of Estimates for time:group effect (True value = 1.5)
# LME vs GEE_cont
# ============================================================

library(ggplot2)
library(data.table)

# Load simulation results (assuming 500 replications per scenario)
# Adjust the path if needed
results <- readRDS("simulation_results.rds")

# Filter for continuous outcome methods and keep only converged replications
estimates <- results[method %in% c("LME", "GEE_cont") & converged == TRUE]

# Check the number of estimates
cat("Number of LME estimates:", nrow(estimates[method == "LME"]), "\n")
cat("Number of GEE_cont estimates:", nrow(estimates[method == "GEE_cont"]), "\n")

# Create a combined histogram with faceting
p <- ggplot(estimates, aes(x = estimate, fill = method)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~ method, ncol = 2, scales = "free_y") +
  labs(
    title = "Distribution of Estimated time×group Effect (True = 1.5)",
    x = "Estimate",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

# Display the plot
print(p)

# Save high-resolution version for paper
ggsave("histogram_LME_GEE_estimates.png", p, width = 10, height = 6, dpi = 300, bg = "white")

# Optional: summary statistics
summary_stats <- estimates[, .(
  Mean = mean(estimate),
  SD = sd(estimate),
  `2.5%` = quantile(estimate, 0.025),
  `97.5%` = quantile(estimate, 0.975),
  Pct_below_1.5 = mean(estimate < 1.5) * 100,
  Pct_above_1.5 = mean(estimate > 1.5) * 100
), by = method]

cat("\nSummary statistics:\n")
print(summary_stats)