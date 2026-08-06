# =============================================================================
# twfe_pretreatment_averaging.R
#
# Simulation study confirming a basic TWFE identity.
#
# Claim: With two PRE-treatment periods (t = 1, 2) and one POST-treatment
# period (t = 3), the two-way fixed-effects (TWFE) DiD coefficient uses the
# AVERAGE of the two pre-treatment periods as the baseline for its "first
# difference" -- not the last pre-period alone.
#
# Formally, with treated-group cell means a1,a2,a3 and control-group cell means
# c1,c2,c3, we prove and then verify by simulation that
#
#     beta_TWFE = [ a3 - (a1 + a2)/2 ] - [ c3 - (c1 + c2)/2 ].
#
# The bracketed treated term is the treated group's first difference: its post
# value minus the mean of its two pre-treatment values. The identity holds
# EXACTLY in every draw (not just on average), because the treatment dummy only
# varies at the group x time level, so FWL collapses the estimator onto the six
# cell means regardless of the idiosyncratic noise.
#
# Base R only (lm with factor FEs); no external data or packages required.
# Writes output/twfe_pretreatment_averaging.csv.
# =============================================================================

set.seed(20260806)

# ------------------------------------------------------------------------------
# Data-generating process
# ------------------------------------------------------------------------------
# One treated group, one control group; n units each; periods t = 1, 2, 3.
# Treatment switches on for the treated group in period 3 only.
#
#   Y_it = alpha_i + lambda_t + tau * D_it + eps_it
#
# alpha_i : unit fixed effect (persistent level differences)
# lambda_t: common time shocks (parallel trends hold by construction)
# tau     : the true treatment effect, applied to treated units at t = 3
# D_it    : 1 iff unit i is treated AND t == 3
# ------------------------------------------------------------------------------
simulate_panel <- function(n = 200, tau = 3, time_fx = c(0, 1.5, 4),
                           treat_gap = 5, sd_unit = 2, sd_noise = 1) {
  units  <- seq_len(2 * n)
  treated <- units <= n                       # first n units are treated
  alpha  <- rnorm(2 * n, mean = ifelse(treated, treat_gap, 0), sd = sd_unit)

  grid <- expand.grid(id = units, t = 1:3)
  grid$treated <- grid$id <= n
  grid$alpha   <- alpha[grid$id]
  grid$lambda  <- time_fx[grid$t]
  grid$D       <- as.integer(grid$treated & grid$t == 3)
  grid$Y       <- grid$alpha + grid$lambda + tau * grid$D +
                  rnorm(nrow(grid), sd = sd_noise)
  grid
}

# ------------------------------------------------------------------------------
# Estimators
# ------------------------------------------------------------------------------

# TWFE DiD: unit FE + time FE + treatment dummy.
twfe_beta <- function(df) {
  fit <- lm(Y ~ D + factor(id) + factor(t), data = df)
  unname(coef(fit)["D"])
}

# The "averaging" formula built directly from the six group x time cell means.
averaging_beta <- function(df) {
  cell <- tapply(df$Y, list(df$treated, df$t), mean)   # rows: FALSE/TRUE, cols: 1,2,3
  c1 <- cell["FALSE", "1"]; c2 <- cell["FALSE", "2"]; c3 <- cell["FALSE", "3"]
  a1 <- cell["TRUE",  "1"]; a2 <- cell["TRUE",  "2"]; a3 <- cell["TRUE",  "3"]
  treated_fd <- a3 - (a1 + a2) / 2      # treated first difference (avg baseline)
  control_fd <- c3 - (c1 + c2) / 2      # control first difference (avg baseline)
  list(beta = treated_fd - control_fd,
       treated_fd = treated_fd, control_fd = control_fd,
       a1 = a1, a2 = a2, a3 = a3, c1 = c1, c2 = c2, c3 = c3)
}

# Contrast: naive DiD that uses only the LAST pre-period (t = 2) as the baseline.
last_pre_beta <- function(df) {
  cell <- tapply(df$Y, list(df$treated, df$t), mean)
  (cell["TRUE", "3"] - cell["TRUE", "2"]) - (cell["FALSE", "3"] - cell["FALSE", "2"])
}

# ------------------------------------------------------------------------------
# 1. Single draw: show the estimator and the formula agree to machine precision
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(" TWFE with 2 pre-periods + 1 post-period: baseline averaging\n")
cat("============================================================\n\n")

df  <- simulate_panel()
b_twfe <- twfe_beta(df)
avg    <- averaging_beta(df)

cat("== 1. Single simulated panel (n = 200 per group) ==\n")
cat(sprintf("  treated cell means : a1=%.4f  a2=%.4f  a3=%.4f\n", avg$a1, avg$a2, avg$a3))
cat(sprintf("  control cell means : c1=%.4f  c2=%.4f  c3=%.4f\n", avg$c1, avg$c2, avg$c3))
cat(sprintf("  treated first diff  a3 - (a1+a2)/2 = %.6f\n", avg$treated_fd))
cat(sprintf("  control first diff  c3 - (c1+c2)/2 = %.6f\n", avg$control_fd))
cat("  ------------------------------------------------------------\n")
cat(sprintf("  TWFE coefficient (lm)            : %.10f\n", b_twfe))
cat(sprintf("  Averaging formula                : %.10f\n", avg$beta))
cat(sprintf("  |difference|                     : %.3e\n", abs(b_twfe - avg$beta)))
cat(sprintf("  Identity holds (< 1e-8)?         : %s\n\n", abs(b_twfe - avg$beta) < 1e-8))

# For reference: what the last-pre-period-only DiD would give instead.
b_lastpre <- last_pre_beta(df)
cat(sprintf("  For contrast, DiD using ONLY t=2 as baseline: %.6f\n", b_lastpre))
cat("  (differs from TWFE -- TWFE averages t=1 and t=2, it does not\n")
cat("   drop t=1)\n\n")

# ------------------------------------------------------------------------------
# 2. Monte Carlo: the identity is exact in EVERY replication
# ------------------------------------------------------------------------------
cat("== 2. Monte Carlo over 1,000 independent panels ==\n")
n_rep <- 1000
diffs <- numeric(n_rep)
twfe_draws <- numeric(n_rep)
for (r in seq_len(n_rep)) {
  d <- simulate_panel()
  bt <- twfe_beta(d)
  ba <- averaging_beta(d)$beta
  twfe_draws[r] <- bt
  diffs[r] <- abs(bt - ba)
}
cat(sprintf("  max |TWFE - averaging formula| over 1,000 draws : %.3e\n", max(diffs)))
cat(sprintf("  identity holds in all draws (< 1e-8)?           : %s\n", all(diffs < 1e-8)))
cat(sprintf("  mean TWFE estimate (true tau = 3)               : %.4f\n", mean(twfe_draws)))
cat(sprintf("  sd of TWFE estimate                             : %.4f\n\n", sd(twfe_draws)))

# ------------------------------------------------------------------------------
# 3. Why the averaging matters: a linear pre-trend
# ------------------------------------------------------------------------------
# Averaging the two pre-periods only recovers tau when the pre-periods share the
# treated group's period-3 baseline. Add a treated-specific linear trend g per
# period (a violation of parallel trends). The averaged baseline sits BELOW the
# correctly extrapolated t=3 counterfactual, so TWFE is biased by exactly the
# trend carried from the averaged pre-mean (t = 1.5) to the post period (t = 3),
# i.e. a distance of 1.5 periods: bias = 1.5 * g.
# ------------------------------------------------------------------------------
cat("== 3. Interpretation under a treated-group linear pre-trend ==\n")
simulate_trend <- function(n = 200, tau = 3, g = 1.0, sd_noise = 1) {
  df <- simulate_panel(n = n, tau = tau, sd_noise = sd_noise)
  df$Y <- df$Y + ifelse(df$treated, g * df$t, 0)   # add treated linear trend g*t
  df
}
for (g in c(0.0, 0.5, 1.0, 2.0)) {
  set.seed(1000 + round(g * 10))
  bt <- mean(replicate(200, twfe_beta(simulate_trend(g = g))))
  cat(sprintf("  trend g = %.1f : mean TWFE = %.4f  (true tau = 3, predicted bias 1.5*g = %.3f => %.3f)\n",
              g, bt, 1.5 * g, 3 + 1.5 * g))
}
cat("  => The bias equals 1.5 * g: the averaged pre-mean sits at t = 1.5, so\n")
cat("     the trend is extrapolated over 1.5 periods to reach t = 3. This is\n")
cat("     the concrete cost of averaging both pre-periods instead of using the\n")
cat("     last one (which would extrapolate over only 1 period).\n\n")

# ------------------------------------------------------------------------------
# Persist a tidy summary
# ------------------------------------------------------------------------------
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
summary_df <- data.frame(
  quantity = c("twfe_beta_lm", "averaging_formula", "abs_diff_single",
               "mc_max_abs_diff", "mc_mean_twfe", "mc_sd_twfe",
               "last_pre_only_diff"),
  value = c(b_twfe, avg$beta, abs(b_twfe - avg$beta),
            max(diffs), mean(twfe_draws), sd(twfe_draws), b_lastpre)
)
out_path <- file.path(out_dir, "twfe_pretreatment_averaging.csv")
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("Wrote %s\n", out_path))
