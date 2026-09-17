# =============================================================================
# stacked_loo_both.R
#
# Leave-one-cohort-out for the stacked regression under the CORRECTED interactive
# fixed effects (agency x stack + year x cohort), run on BOTH stacks so the
# sensitivity of the pooled coefficient to any single sub-experiment can be
# compared across control-attachment rules:
#
#   * shipped-11        : the authors' stacked_fatal.csv. Controls attach only
#                         where g-4>=2000, so 2000-2003 are treated-only; the
#                         plain pooled coef = -0.0989. The `weights` (entropy-
#                         balancing) column is present, so the m4/Table-3 spec is
#                         run here.
#   * reconstructed-15  : all 15 change-cohorts (>=2000) given the 741 clean
#                         controls (unweighted Table A5 col-3 construction). The
#                         ebal weights are NOT reconstructible, so the covariate
#                         spec here is m3 (covariates, unweighted), not m4.
#
# For each stack, each row drops one cohort's entire sub-experiment and refits.
# `plain` (no controls, no weights) is the decomposition target and the apples-
# to-apples comparison across the two stacks. Clustered SEs on agency.id.
# =============================================================================

suppressMessages({library(lfe)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
covs <- "log.pop + log.med.inc + pct.white + pct.white.officers.imputed"

# ---- shipped stack (has ebal weights) --------------------------------------
ship <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
ship$agency_stack <- paste(ship$agency.id, ship$cohort, sep = "__")

# ---- reconstructed 15-cohort stack (no ebal weights) -----------------------
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
build15 <- function() {
  cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= 2000]))
  nochange <- dta[dta$change.type == "No Change", ]
  parts <- list()
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g, ]; tr$cohort <- g; tr$treat <- 1L
    parts[[length(parts) + 1L]] <- tr
    ct <- nochange; ct$cohort <- g; ct$treat <- 0L; parts[[length(parts) + 1L]] <- ct
  }
  d <- do.call(rbind, parts)
  d$year.cohort  <- as.numeric(paste0(d$year, d$cohort))
  d$agency_stack <- paste(d$agency.id, d$cohort, sep = "__")
  d
}
rec15 <- build15()

# ---- fitters ---------------------------------------------------------------
grab <- function(m, v = "no.req") c(est = unname(coef(m)[v]), se = unname(m$se[v]))
fit_plain <- function(d)
  felm(any.fatalities ~ no.req | agency_stack + year.cohort | 0 | agency.id, data = d)
fit_cov <- function(d, weighted) {
  f <- as.formula(paste0("any.fatalities ~ no.req + ", covs,
                         " | agency_stack + year.cohort | 0 | agency.id"))
  if (weighted) felm(f, data = d, weights = d$weights) else felm(f, data = d)
}

loo_table <- function(d, weighted, cov_label) {
  cohorts <- sort(unique(d$cohort))
  one <- function(sub, lab, ns) {
    p <- grab(fit_plain(sub)); m <- grab(fit_cov(sub, weighted))
    data.frame(dropped = lab, cohorts = ns,
               plain_est = p["est"], plain_se = p["se"],
               cov_est = m["est"], cov_se = m["se"],
               cov_label = cov_label, row.names = NULL)
  }
  base <- one(d, "(none - full)", length(cohorts))
  loo  <- do.call(rbind, lapply(cohorts, function(g)
    one(d[d$cohort != g, ], as.character(g), length(cohorts) - 1L)))
  out <- rbind(base, loo)
  out$plain_delta <- out$plain_est - base$plain_est
  out$cov_delta   <- out$cov_est   - base$cov_est
  out
}

ship_loo <- transform(loo_table(ship,  TRUE,  "m4 (cov+weights)"), stack = "shipped-11")
rec_loo  <- transform(loo_table(rec15, FALSE, "m3 (cov, no wts)"), stack = "reconstructed-15")

options(width = 210)
show <- function(x) transform(x,
  plain_est = round(plain_est, 4), plain_se = round(plain_se, 4), plain_delta = round(plain_delta, 4),
  cov_est = round(cov_est, 4), cov_se = round(cov_se, 4), cov_delta = round(cov_delta, 4))
cat("== LOO, shipped-11 stack (corrected FE; cov spec = m4 weighted) ==\n")
print(show(ship_loo)[, c("dropped","cohorts","plain_est","plain_se","plain_delta","cov_est","cov_se","cov_delta")], row.names = FALSE)
cat("\n== LOO, reconstructed-15 stack (corrected FE; cov spec = m3 unweighted) ==\n")
print(show(rec_loo)[, c("dropped","cohorts","plain_est","plain_se","plain_delta","cov_est","cov_se","cov_delta")], row.names = FALSE)

res <- rbind(ship_loo, rec_loo)
write.csv(res, file.path(OUT_DIR, "stacked_loo_both.csv"), row.names = FALSE)

# ---- LaTeX (both panels) ---------------------------------------------------
f3 <- function(x) sprintf("%+.3f", x)
panel <- function(sub, cap, lab, covhead) {
  L <- c("\\begin{table}[t]\\centering",
         paste0("\\caption{", cap, "}"),
         paste0("\\label{", lab, "}"),
         "\\begin{tabular}{lrcc}",
         "\\toprule",
         paste0("Cohort dropped & Cohorts & Plain & ", covhead, " \\\\"),
         "\\midrule")
  for (i in seq_len(nrow(sub))) {
    lab_i <- if (i == 1) "\\textit{None (full)}" else sub$dropped[i]
    L <- c(L, sprintf("%s & %d & %s (%.3f) & %s (%.3f) \\\\",
                      lab_i, sub$cohorts[i], f3(sub$plain_est[i]), sub$plain_se[i],
                      f3(sub$cov_est[i]), sub$cov_se[i]))
    if (i == 1) L <- c(L, "\\midrule")
  }
  c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}
lines <- c(
  panel(ship_loo,
    paste0("Leave-one-cohort-out on the shipped stacked dataset, corrected ",
           "interactive fixed effects. Each row drops one cohort's sub-experiment ",
           "and refits \\texttt{no.req}; clustered SEs in parentheses. Plain has no ",
           "controls; m4 adds the paper's four covariates and the entropy-balancing ",
           "weights (Table 3)."),
    "tab:loo-shipped", "m4 (Table 3)"),
  "",
  panel(rec_loo,
    paste0("Leave-one-cohort-out on the reconstructed 15-cohort stack (all change-",
           "cohorts $\\ge 2000$ given the 741 clean controls). The entropy-balancing ",
           "weights are not reconstructible, so the covariate column is m3 ",
           "(covariates, unweighted). Same corrected fixed effects and clustering."),
    "tab:loo-recon15", "m3 (cov, no wts)"))
writeLines(lines, file.path(OUT_DIR, "stacked_loo_both.tex"))
cat("\nWrote output/stacked_loo_both.csv and .tex\n")
