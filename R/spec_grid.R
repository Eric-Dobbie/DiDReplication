# =============================================================================
# spec_grid.R
#
# Specification grid reproducing P&P's fatal-encounters stacked table (their
# "Table 3" in police_residency_main.R; columns 3-4 = the Stacked Approach with
# controls). Estimates no.req on any.fatalities in stacked_fatal.csv across all
# combinations of:
#   FE         : additive agency.id + year.cohort   vs   agency.id x cohort
#                (= agency_stack) + year.cohort
#   covariates : none   vs   the 4 P&P controls
#                (log.pop + log.med.inc + pct.white + pct.white.officers.imputed)
#   weights    : unweighted   vs   native entropy-balancing `weights` column
# SE clustered on agency.id. NB felm's weighted `$se` is the NON-clustered
# analytic SE; the clustered SE is `$cse` -- we report $cse throughout.
#
# Paper targets (from the PDF, p.25): col 3 = -0.103 (0.039), col 4 = -0.091
# (0.037). Both reproduce in the ADDITIVE-FE + 4-controls row (col 3 unweighted,
# col 4 weighted), confirming the paper uses the shared agency.id FE, not the
# leakage-free agency x stack FE (see §9). [The col-3 SE is 0.039, not 0.037.]
# =============================================================================
suppressMessages(library(lfe))
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
COV <- "log.pop + log.med.inc + pct.white + pct.white.officers.imputed"

cell <- function(fe, cov, wt) {
  unit <- if (fe == "additive") "agency.id" else "agency_stack"
  rhs  <- if (cov) paste0("no.req + ", COV) else "no.req"
  f <- as.formula(sprintf("any.fatalities ~ %s | %s + year.cohort | 0 | agency.id", rhs, unit))
  m <- if (wt) felm(f, data = st, weights = st$weights) else felm(f, data = st)
  b <- unname(coef(m)["no.req"]); s <- unname(m$cse["no.req"])          # clustered SE
  hit <- ""
  if (round(b,3) == -0.103 && round(s,3) == 0.039) hit <- "col 3"
  if (round(b,3) == -0.091 && round(s,3) == 0.037) hit <- "col 4"
  data.frame(fe = fe, covariates = ifelse(cov, "4 P&P", "none"),
             weights = ifelse(wt, "weighted", "unweighted"),
             coef = b, cse = s, N = m$N, match = hit, stringsAsFactors = FALSE)
}
grid <- expand.grid(fe = c("additive", "interacted"), cov = c(FALSE, TRUE),
                    wt = c(FALSE, TRUE), stringsAsFactors = FALSE)
res <- do.call(rbind, Map(cell, grid$fe, grid$cov, grid$wt))
res <- res[order(res$fe, res$covariates, res$weights), ]

options(width = 200)
print(transform(res, coef = round(coef, 4), cse = round(cse, 4)), row.names = FALSE)
write.csv(transform(res, coef = round(coef, 6), cse = round(cse, 6)),
          file.path(OUT_DIR, "spec_grid.csv"), row.names = FALSE)

# ---- LaTeX -----------------------------------------------------------------
fe_lab <- function(x) ifelse(x == "additive",
  "additive $\\text{agency}+\\text{yr}\\times\\text{coh}$",
  "$\\text{agency}\\times\\text{coh}+\\text{yr}\\times\\text{coh}$")
L <- c(
  "\\begin{table}[t]\\centering",
  "\\caption{Specification grid for the stacked fatal-encounters regression on",
  "\\texttt{stacked\\_fatal.csv}: \\texttt{no.req} on \\texttt{any.fatalities} across",
  "fixed effects (additive \\texttt{agency.id} vs \\texttt{agency.id}$\\times$\\texttt{cohort},",
  "both with \\texttt{year.cohort}), covariates (none vs P\\&P's four controls), and",
  "weights (unweighted vs the native entropy-balancing column). SEs clustered on",
  "\\texttt{agency.id}. P\\&P Table columns 3 and 4 (both with controls) are",
  "$-0.103\\,(0.039)$ and $-0.091\\,(0.037)$; both reproduce in the additive-FE +",
  "controls row (col.\\ 3 unweighted, col.\\ 4 weighted), confirming the paper's",
  "shared \\texttt{agency.id} FE rather than the leakage-free $\\text{agency}\\times",
  "\\text{stack}$ FE.}",
  "\\label{tab:spec-grid}",
  "\\begin{tabular}{lllrrrl}",
  "\\toprule",
  "Fixed effects & Covariates & Weights & Coef. & Cluster SE & $N$ & Matches \\\\",
  "\\midrule")
for (i in seq_len(nrow(res))) {
  mtag <- if (nzchar(res$match[i])) sprintf("\\textbf{%s}", res$match[i]) else ""
  L <- c(L, sprintf("%s & %s & %s & %s%.4f%s & %.4f & %s & %s \\\\",
                    fe_lab(res$fe[i]), res$covariates[i], res$weights[i],
                    if (nzchar(res$match[i])) "\\textbf{" else "", res$coef[i],
                    if (nzchar(res$match[i])) "}" else "",
                    res$cse[i], formatC(res$N[i], format = "d", big.mark = ","), mtag))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(L, file.path(OUT_DIR, "spec_grid.tex"))
cat("\nWrote output/spec_grid.{csv,tex}\n")
