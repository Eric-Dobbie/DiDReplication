# =============================================================================
# stacked_loo.R
#
# Leave-one-cohort-out for the stacked regression under the CORRECTED
# interactive fixed effects (agency_stack = agency.id x cohort, + year.cohort).
# For each treatment cohort g, drop that cohort's entire sub-experiment (all rows
# with cohort == g) and refit the pooled no.req coefficient on the remaining
# stacks. Reports two specs:
#   * plain     : any.fatalities ~ no.req                | agency_stack + year.cohort
#                 (the decomposition target; full-sample = -0.0989)
#   * m4 (Tab3) : any.fatalities ~ no.req + <4 covariates> | agency_stack + year.cohort,
#                 weights = balancing weights (full-sample = -0.0891)
# Both cluster SEs on the real agency.id. Data unchanged (P&P full panel).
# =============================================================================

suppressMessages({library(lfe); library(dplyr)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")

st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
covs <- "log.pop + log.med.inc + pct.white + pct.white.officers.imputed"

fit_plain <- function(d)
  felm(any.fatalities ~ no.req | agency_stack + year.cohort | 0 | agency.id, data = d)
fit_m4 <- function(d)
  felm(as.formula(paste0("any.fatalities ~ no.req + ", covs,
                         " | agency_stack + year.cohort | 0 | agency.id")),
       data = d, weights = d$weights)
grab <- function(m) c(est = unname(coef(m)["no.req"]), se = unname(m$se["no.req"]))

run <- function(d, drop_label, n_stacks) {
  p <- grab(fit_plain(d)); m <- grab(fit_m4(d))
  data.frame(dropped = drop_label, stacks = n_stacks,
             N_plain = fit_plain(d)$N,
             plain_est = p["est"], plain_se = p["se"],
             m4_est = m["est"], m4_se = m["se"], row.names = NULL)
}

cohorts <- sort(unique(st$cohort))
base <- run(st, "(none - full sample)", length(cohorts))
loo  <- do.call(rbind, lapply(cohorts, function(g)
  run(st[st$cohort != g, ], as.character(g), length(cohorts) - 1L)))
res <- rbind(base, loo)
res$plain_delta <- res$plain_est - base$plain_est
res$m4_delta    <- res$m4_est - base$m4_est

options(width = 200)
cat("== Leave-one-cohort-out, corrected (agency x stack) FE ==\n")
print(transform(res,
        plain_est = round(plain_est, 4), plain_se = round(plain_se, 4),
        m4_est = round(m4_est, 4), m4_se = round(m4_se, 4),
        plain_delta = round(plain_delta, 4), m4_delta = round(m4_delta, 4)),
      row.names = FALSE)

write.csv(res, file.path(OUT_DIR, "stacked_loo.csv"), row.names = FALSE)

# ---- LaTeX table -----------------------------------------------------------
f3 <- function(x) sprintf("%+.3f", x)
lines <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Leave-one-cohort-out for the stacked regression with the corrected",
  "interactive fixed effects ($\\text{agency}\\times\\text{stack} + \\text{year}\\times\\text{cohort}$).",
  "Each row drops one treatment cohort's sub-experiment and refits the pooled",
  "coefficient on \\texttt{no.req}; clustered SEs in parentheses. `Plain' has no",
  "controls (the decomposition target); `m4' adds the paper's four covariates and",
  "the entropy-balancing weights (Table 3).}",
  "\\label{tab:stacked-loo}",
  "\\begin{tabular}{lrcc}",
  "\\toprule",
  "Cohort dropped & Cohorts & Plain & m4 (Table 3) \\\\",
  "\\midrule")
for (i in seq_len(nrow(res))) {
  lab <- if (i == 1) "\\textit{None (full)}" else res$dropped[i]
  lines <- c(lines, sprintf("%s & %d & %s (%.3f) & %s (%.3f) \\\\",
                            lab, res$stacks[i], f3(res$plain_est[i]), res$plain_se[i],
                            f3(res$m4_est[i]), res$m4_se[i]))
  if (i == 1) lines <- c(lines, "\\midrule")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(OUT_DIR, "stacked_loo.tex"))
cat("\nWrote output/stacked_loo.csv and output/stacked_loo.tex\n")
