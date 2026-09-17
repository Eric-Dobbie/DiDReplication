# =============================================================================
# stacked_full15.R
#
# The stacked specification run "exactly as the paper describes it": the FULL
# 15-cohort stack in which EVERY change-cohort (>=2000) is given all 741 pure
# control cities across the full panel (no window) -- i.e. controls are added to
# the 2000-2003 stacks that the shipped stacked_fatal.csv leaves treated-only.
# This reproduces the paper's own Table A.5 col-3 sample (N = 233,520).
#
# Estimates no.req on any.fatalities under both the interacted (agency x stack +
# year x cohort) and additive (agency.id + year x cohort) FE, with and without
# the four P&P controls; SE clustered on agency.id; unweighted (the entropy-
# balancing weights are not defined for the added early-cohort controls).
#
# Writes:
#   output/stacked_full15.csv         -- the four-cell summary
#   output/stacked_full15.tex         -- spec-summary table (coef / SE / N)
#   output/stacked_full15_regtable.tex-- standard regression-output table
# =============================================================================
suppressMessages(library(lfe))
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
COVv <- c("log.pop", "log.med.inc", "pct.white", "pct.white.officers.imputed")
COV  <- paste(COVv, collapse = " + ")

# ---- build the full 15-cohort stack (controls on every cohort, full panel) --
cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= 2000]))
nochange <- dta[dta$change.type == "No Change", ]
parts <- list()
for (g in cohorts) {
  tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g, ]; tr$cohort <- g; tr$treat <- 1L
  parts[[length(parts) + 1L]] <- tr
  ct <- nochange; ct$cohort <- g; ct$treat <- 0L; parts[[length(parts) + 1L]] <- ct
}
st <- do.call(rbind, parts)
st$year.cohort  <- as.numeric(paste0(st$year, st$cohort))
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")

fit <- function(unit, cov) {
  rhs <- if (cov) paste0("no.req + ", COV) else "no.req"
  felm(as.formula(sprintf("any.fatalities ~ %s | %s + year.cohort | 0 | agency.id", rhs, unit)),
       data = st)
}
info <- function(m, cov) {
  ct <- summary(m)$coefficients["no.req", ]
  # estimation sample for mean outcome / agency count
  samp <- !is.na(st$any.fatalities)
  if (cov) samp <- samp & stats::complete.cases(st[, COVv])
  list(coef = unname(ct["Estimate"]), se = unname(ct["Cluster s.e."]),
       p = unname(ct["Pr(>|t|)"]), N = m$N,
       nagency = length(unique(st$agency.id[samp])),
       ybar = mean(st$any.fatalities[samp]))
}
specs <- list(
  list(unit="agency_stack", cov=FALSE, fe="interacted"),
  list(unit="agency_stack", cov=TRUE,  fe="interacted"),
  list(unit="agency.id",    cov=FALSE, fe="additive"),
  list(unit="agency.id",    cov=TRUE,  fe="additive"))
R <- lapply(specs, function(s) c(s, info(fit(s$unit, s$cov), s$cov)))

# ---- CSV + spec-summary table ----------------------------------------------
summ <- do.call(rbind, lapply(R, function(r) data.frame(
  FE = r$fe, covariates = ifelse(r$cov, "4 P&P", "none"),
  coef = round(r$coef,4), cluster_se = round(r$se,4), N = r$N,
  n_agencies = r$nagency, mean_outcome = round(r$ybar,4))))
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT_DIR, "stacked_full15.csv"), row.names = FALSE)

star <- function(p) if (is.na(p)) "" else if (p < 0.05) "$^{*}$" else ""
f3 <- function(x) sprintf("%.3f", x)
Lspec <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Stacked fatal-encounters regression on the \\emph{full} 15-cohort stack:",
  "every change-cohort ($\\ge 2000$) is given all 741 pure control cities over the full",
  "panel (controls added to 2000--2003; no window), reproducing the paper's Table~A.5",
  "col.\\ 3 sample ($N=233{,}520$). \\texttt{no.req} on \\texttt{any.fatalities}; SEs",
  "clustered on \\texttt{agency.id}; unweighted. The interacted-FE + controls cell is",
  "the specification ``as described'' ($-0.092$); the published Table~5 col.\\ 3",
  "($-0.103$) instead uses the additive \\texttt{agency.id} FE on the shipped stack",
  "($N=171{,}444$, controls only on 2004--2020).}",
  "\\label{tab:stacked-full15}",
  "\\begin{tabular}{llrrr}",
  "\\toprule",
  "Fixed effects & Covariates & Coef. & Cluster SE & $N$ \\\\",
  "\\midrule")
for (r in R) {
  fe_l <- if (r$fe == "interacted") "agency$\\times$stack + yr$\\times$coh"
          else "agency.id + yr$\\times$coh"
  Lspec <- c(Lspec, sprintf("%s & %s & %s%s & %s & %s \\\\",
             fe_l, ifelse(r$cov, "4 P\\&P", "none"),
             f3(r$coef), star(r$p), f3(r$se),
             formatC(r$N, format="d", big.mark=",")))
}
Lspec <- c(Lspec, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(Lspec, file.path(OUT_DIR, "stacked_full15.tex"))

# ---- standard regression-output table (stargazer-style, 4 model columns) ----
cc <- function(i) R[[i]]
row_fe <- function(name, vals) paste0(name, " & ", paste(vals, collapse = " & "), " \\\\")
Lreg <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Residency Requirements and Probability of Fatal Encounter --- full",
  "15-cohort stacked approach (controls added to all cohorts, no window). Dependent",
  "variable: \\texttt{any.fatalities}. Robust standard errors clustered by city in",
  "parentheses. $^{*}$p$<$0.05.}",
  "\\label{tab:fatal-full15-reg}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & (1) & (2) & (3) & (4) \\\\",
  "\\midrule",
  sprintf("Requirement Dropped & %s%s & %s%s & %s%s & %s%s \\\\",
          f3(cc(1)$coef), star(cc(1)$p), f3(cc(2)$coef), star(cc(2)$p),
          f3(cc(3)$coef), star(cc(3)$p), f3(cc(4)$coef), star(cc(4)$p)),
  sprintf(" & (%s) & (%s) & (%s) & (%s) \\\\",
          f3(cc(1)$se), f3(cc(2)$se), f3(cc(3)$se), f3(cc(4)$se)),
  "\\midrule",
  row_fe("Agency $\\times$ Stack FEs", c("Yes","Yes","","")),
  row_fe("Agency FEs",                 c("","","Yes","Yes")),
  row_fe("Year $\\times$ Cohort FEs",  c("Yes","Yes","Yes","Yes")),
  row_fe("Controls",                   c("","Yes","","Yes")),
  row_fe("Mean Outcome", sprintf("%.3f", sapply(1:4, function(i) cc(i)$ybar))),
  row_fe("Num. Agencies", sapply(1:4, function(i) cc(i)$nagency)),
  row_fe("Observations", formatC(sapply(1:4, function(i) cc(i)$N), format="d", big.mark=",")),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}")
writeLines(Lreg, file.path(OUT_DIR, "stacked_full15_regtable.tex"))
cat("\nWrote output/stacked_full15.{csv,tex} and stacked_full15_regtable.tex\n")
