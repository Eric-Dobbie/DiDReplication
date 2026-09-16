# =============================================================================
# verify_cs_sa_equivalence.R
#
# Verifies the claim used in the paper's figures (plot_cohort_level.R,
# plot_weight_vs_beta.R) that, under
#     (a) never-treated units as the control group, and
#     (b) a balanced panel (every cohort x period cell populated),
# the Callaway & Sant'Anna (CS) and Sun & Abraham (SA) estimators reduce to the
# SAME underlying 2x2 DiD comparisons and their aggregation weights collapse to
# the same values -- so both the disaggregated "atoms" AND the aggregate ATT
# coincide numerically (not merely approximately).
#
# The script is self-contained: it SIMULATES a balanced staggered-adoption panel
# with never-treated controls and known, heterogeneous, dynamic effects, so the
# result does not depend on the (git-ignored) P&P data being present. If you set
# DIDREP_DATA to a folder containing dta.csv, an appendix block at the bottom
# repeats the headline aggregate comparison on the real fatal-encounters sample.
#
# It then BREAKS each assumption one at a time and confirms divergence:
#   (1) unbalance the panel  (drop cells)          -> aggregate weights diverge
#   (2) switch CS to control_group = "notyettreated" -> cell estimates diverge
#
# Packages: did (CS), fixest (SA).  Run:  Rscript R/verify_cs_sa_equivalence.R
# =============================================================================

suppressMessages({
  library(did)
  library(fixest)
})

set.seed(20240517)
options(width = 110)

`%||%` <- function(a, b) if (is.null(a)) b else a
OUT_DIR <- Sys.getenv("DIDREP_OUT", unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 0. Simulate a BALANCED staggered panel with never-treated controls
# -----------------------------------------------------------------------------
# cohorts: first-treatment year in {2004, 2007, 2010}; never-treated coded 0.
# Balanced: every unit is observed in every year 2000..2014 (no missing cells).
# Effects are cohort- and event-time-specific (heterogeneous + dynamic), so CS,
# SA and TWFE are NOT mechanically equal -- the CS==SA match below is a real
# coincidence of construction, not a degenerate "all effects zero" artifact.
# -----------------------------------------------------------------------------
make_panel <- function(n_per_group = 40, years = 2000:2014,
                       cohorts = c(2004, 2007, 2010)) {
  groups <- c(0, cohorts)                       # 0 = never-treated
  units  <- data.frame(
    id = seq_len(n_per_group * length(groups)),
    g  = rep(groups, each = n_per_group)
  )
  units$fe_u <- rnorm(nrow(units), 0, 1)        # unit fixed effect
  yr_fe <- setNames(rnorm(length(years), 0, 0.5), years)

  # heterogeneous dynamic effect: cohort-specific slope in event time
  coh_slope <- c("2004" = 0.06, "2007" = -0.03, "2010" = 0.10)

  panel <- do.call(rbind, lapply(seq_len(nrow(units)), function(i) {
    g <- units$g[i]
    e <- if (g == 0) rep(NA_integer_, length(years)) else years - g
    treat_post <- !is.na(e) & e >= 0
    te <- rep(0, length(years))
    if (g != 0) te[treat_post] <- coh_slope[as.character(g)] * (e[treat_post] + 1)
    data.frame(
      id   = units$id[i],
      g    = g,
      year = years,
      y    = units$fe_u[i] + yr_fe[as.character(years)] + te +
             rnorm(length(years), 0, 0.4)
    )
  }))
  panel
}

panel <- make_panel()
stopifnot(!any(is.na(panel$y)))
# confirm balance: every id has every year exactly once
tab <- table(table(panel$id))
cat("== simulated panel ==  units:", length(unique(panel$id)),
    " years:", paste(range(panel$year), collapse = "-"),
    " rows:", nrow(panel), "\n")
cat("   balance check (every unit has", max(panel$year) - min(panel$year) + 1,
    "obs):", all(table(panel$id) == (max(panel$year) - min(panel$year) + 1)),
    "\n")
cat("   cohorts (g):", paste(sort(unique(panel$g)), collapse = ", "),
    "  (0 = never-treated)\n\n")

# -----------------------------------------------------------------------------
# Helper: CS ATT(g,t) as a tidy (g, t, e=t-g, estimate) frame
# -----------------------------------------------------------------------------
cs_cells <- function(data, control_group = "nevertreated",
                     base_period = "universal") {
  mp <- did::att_gt(yname = "y", tname = "year", idname = "id", gname = "g",
                    xformla = ~1, data = data, control_group = control_group,
                    base_period = base_period, bstrap = FALSE, cband = FALSE)
  agg <- did::aggte(mp, type = "simple", na.rm = TRUE)
  cells <- data.frame(g = mp$group, t = mp$t, att = mp$att, se = mp$se)
  cells$e <- cells$t - cells$g
  list(cells = cells, overall = as.numeric(agg$overall.att),
       overall_se = as.numeric(agg$overall.se), mp = mp)
}

# -----------------------------------------------------------------------------
# Helper: SA CATT(g,e) as a tidy (g, e, estimate) frame  (never-treated -> ref)
# -----------------------------------------------------------------------------
sa_cells <- function(data, ref_c = 10000) {
  d <- data
  d$coh <- ifelse(d$g == 0, ref_c, d$g)
  res <- fixest::feols(y ~ sunab(coh, year) | id + year, data = d,
                       cluster = ~id, notes = FALSE)
  att <- fixest::coeftable(summary(res, agg = "att"))["ATT", ]
  ct  <- summary(res, agg = FALSE)$coeftable
  cn  <- rownames(ct)
  mm  <- regmatches(cn, regexec("::(-?[0-9]+):.*::([0-9]+)$", cn))
  e   <- as.integer(vapply(mm, function(z) z[2], character(1)))
  g   <- as.integer(vapply(mm, function(z) z[3], character(1)))
  cells <- data.frame(g = g, e = e, catt = ct[, 1], se = ct[, 2])
  list(cells = cells, overall = as.numeric(att[1]),
       overall_se = as.numeric(att[2]), res = res)
}

# =============================================================================
# 1. BASELINE: balanced panel + never-treated controls  ->  CS == SA
# =============================================================================
cat("=========================================================\n")
cat(" (1) BASELINE: balanced panel, never-treated controls\n")
cat("=========================================================\n")

cs <- cs_cells(panel, control_group = "nevertreated", base_period = "universal")
sa <- sa_cells(panel)

# --- cell-level join on (cohort g, event time e) --------------------------
m <- merge(cs$cells[, c("g", "e", "att")], sa$cells[, c("g", "e", "catt")],
           by = c("g", "e"))
m <- m[order(m$g, m$e), ]
m$abs_diff <- abs(m$att - m$catt)

cat("\nCell-level ATT(g,t=g+e) [CS]  vs  CATT(g,e) [SA]",
    "(joined on cohort g and event time e):\n")
print(transform(m, att = round(att, 6), catt = round(catt, 6),
                abs_diff = signif(abs_diff, 3)), row.names = FALSE)

cat(sprintf("\n  # cells compared : %d\n", nrow(m)))
cat(sprintf("  max |CS - SA|    : %.3e\n", max(m$abs_diff)))
cat(sprintf("  post-period cells (e>=0) max |CS - SA| : %.3e\n",
            max(m$abs_diff[m$e >= 0])))

cat(sprintf("\n  Aggregate ATT   CS = %.10f   SA = %.10f   |diff| = %.3e\n",
            cs$overall, sa$overall, abs(cs$overall - sa$overall)))

baseline_cell_ok <- max(m$abs_diff) < 1e-6
baseline_agg_ok  <- abs(cs$overall - sa$overall) < 1e-6
cat(sprintf("  --> cell-level match: %s   aggregate match: %s\n",
            baseline_cell_ok, baseline_agg_ok))

# =============================================================================
# 2. BREAK (a): UNBALANCE the panel  ->  aggregation weights diverge
# =============================================================================
cat("\n=========================================================\n")
cat(" (2) BREAK BALANCE: drop a chunk of cells (unbalanced panel)\n")
cat("=========================================================\n")

# Drop ~12% of post-treatment rows at random. Two things then break at once:
#   (i)  did::att_gt defaults to allow_unbalanced_panel = FALSE, so it REBALANCES
#        by dropping any unit not observed in every period (watch for its
#        "converting to balanced panel by dropping them" warning) -- SA's
#        regression instead keeps every available row, so the two no longer even
#        use the same observations; and
#   (ii) cohort x event obs counts become uneven, so SA's obs-count weights
#        (colSums(sign(mm))) no longer equal CS's fixed group-share weights (pg).
# Either channel alone is enough to break the numeric equivalence; here both do.
unb <- panel
post_rows <- which(unb$g != 0 & unb$year >= unb$g)
drop <- sample(post_rows, size = floor(0.12 * length(post_rows)))
unb  <- unb[-drop, ]
cat(sprintf("  dropped %d of %d post-treatment rows; panel now unbalanced\n",
            length(drop), length(post_rows)))

cs_u <- cs_cells(unb, control_group = "nevertreated", base_period = "universal")
sa_u <- sa_cells(unb)
cat(sprintf("  Aggregate ATT   CS = %.10f   SA = %.10f   |diff| = %.3e\n",
            cs_u$overall, sa_u$overall, abs(cs_u$overall - sa_u$overall)))
unb_diverges <- abs(cs_u$overall - sa_u$overall) > 1e-6
cat(sprintf("  --> aggregates diverge as expected: %s\n", unb_diverges))

# =============================================================================
# 3. BREAK (b): switch CS control group to not-yet-treated -> cells diverge
# =============================================================================
cat("\n=========================================================\n")
cat(" (3) BREAK CONTROL GROUP: CS not-yet-treated vs SA (still never-treated)\n")
cat("=========================================================\n")

# sunab identifies CATT off the never-treated (reference) cohort. Pointing CS at
# a DIFFERENT control pool (not-yet-treated) changes the 2x2 comparison for every
# interior cohort, so the per-cell estimates no longer coincide.
cs_ny <- cs_cells(panel, control_group = "notyettreated",
                  base_period = "universal")
m2 <- merge(cs_ny$cells[, c("g", "e", "att")], sa$cells[, c("g", "e", "catt")],
            by = c("g", "e"))
m2$abs_diff <- abs(m2$att - m2$catt)
m2_post <- m2[m2$e >= 0, ]
cat(sprintf("  post-period cells (e>=0): max |CS_nyt - SA| = %.3e  (n=%d)\n",
            max(m2_post$abs_diff), nrow(m2_post)))
cat(sprintf("  Aggregate ATT   CS_nyt = %.10f   SA = %.10f   |diff| = %.3e\n",
            cs_ny$overall, sa$overall, abs(cs_ny$overall - sa$overall)))
nyt_diverges <- max(m2_post$abs_diff) > 1e-6
cat(sprintf("  --> cells diverge as expected: %s\n", nyt_diverges))

# =============================================================================
# 4. Summary + machine-readable output
# =============================================================================
summary_df <- data.frame(
  scenario = c("baseline_balanced_nevertreated",
               "unbalanced_nevertreated",
               "balanced_notyettreated"),
  cs_overall = c(cs$overall, cs_u$overall, cs_ny$overall),
  sa_overall = c(sa$overall, sa_u$overall, sa$overall),
  agg_abs_diff = c(abs(cs$overall - sa$overall),
                   abs(cs_u$overall - sa_u$overall),
                   abs(cs_ny$overall - sa$overall)),
  cell_max_abs_diff = c(max(m$abs_diff), NA, max(m2_post$abs_diff)),
  equivalence_expected = c("YES", "no (unbalanced)", "no (diff control group)")
)
cat("\n================= SUMMARY =================\n")
print(summary_df, row.names = FALSE, digits = 6)

write.csv(m, file.path(OUT_DIR, "cs_sa_cell_equivalence.csv"), row.names = FALSE)
write.csv(summary_df, file.path(OUT_DIR, "cs_sa_equivalence_summary.csv"),
          row.names = FALSE)
cat("\nWrote:", file.path(OUT_DIR, "cs_sa_cell_equivalence.csv"), "\n")
cat("Wrote:", file.path(OUT_DIR, "cs_sa_equivalence_summary.csv"), "\n")

# assertions (so the script fails loudly if the claim ever breaks)
stopifnot(
  "baseline cells match"       = baseline_cell_ok,
  "baseline aggregate matches" = baseline_agg_ok,
  "unbalanced diverges"        = unb_diverges,
  "not-yet-treated diverges"   = nyt_diverges
)
cat("\nAll equivalence assertions passed.\n")

# =============================================================================
# 5. (optional) Repeat the headline aggregate check on the real P&P data
# =============================================================================
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "")
if (nzchar(DATA_DIR) && file.exists(file.path(DATA_DIR, "dta.csv"))) {
  cat("\n=========================================================\n")
  cat(" (5) Real P&P fatal-encounters sample (aggregate ATT)\n")
  cat("=========================================================\n")
  dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
  dta$agency.num <- as.numeric(factor(dta$agency.id))
  dta$year.changed <- ave(dta$year.changed, dta$agency.id,
    FUN = function(x) { v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_ })
  dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                      ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987,
                             dta$year.changed))
  d <- dta[!is.na(dta$any.fatalities) &
             (dta$year.changed == 0 | dta$year.changed >= 2001), ]
  names(d)[names(d) == "any.fatalities"] <- "y"
  names(d)[names(d) == "year"]           <- "year"
  names(d)[names(d) == "agency.num"]     <- "id"
  names(d)[names(d) == "year.changed"]   <- "g"
  cs_r <- cs_cells(d, control_group = "nevertreated", base_period = "universal")
  sa_r <- sa_cells(d)
  cat(sprintf("  Aggregate ATT   CS = %.10f   SA = %.10f   |diff| = %.3e\n",
              cs_r$overall, sa_r$overall, abs(cs_r$overall - sa_r$overall)))
}
