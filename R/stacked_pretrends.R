# =============================================================================
# stacked_pretrends.R
#
# Per-sub-experiment (per-cohort) pre-trend test for the stacked design, under
# the CORRECTED interactive fixed effects (agency x stack + year x cohort). For
# each treatment cohort g we fit an event study WITHIN that sub-experiment,
#
#     any.fatalities ~ i(event_time, treat, ref = -1) | agency.id + year
#
# (within one cohort agency_stack == agency.id and year.cohort == year), on the
# +/-K balancing window (K = 4, the paper's window), and jointly Wald-test the
# pre-treatment LEAD coefficients {e < 0, e != -1} = 0. A low p-value flags a
# pre-existing differential trend between that cohort's treated units and its
# clean controls; a high p-value is consistent with parallel pre-trends.
#
# Run on BOTH stacks:
#   * shipped     : the authors' stacked_fatal.csv (controls only for g-4>=2000)
#   * reconstructed-15 : all 15 change-cohorts (>=2000) given the 741 clean
#                        controls (the unweighted Table A5 col-3 construction)
#
# A cohort is testable only if it has (a) clean controls and (b) at least one
# observed treated pre-period lead; otherwise the leads are collinear/absent and
# no test is possible (reported as NA with a reason). Clustered on agency.id.
# =============================================================================

suppressMessages({library(fixest)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
K        <- 4L

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))

# ---- build the two stacks --------------------------------------------------
build_stack <- function(gate15) {
  cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed) & dta$year.changed >= 2000]))
  nochange <- dta[dta$change.type == "No Change", ]
  parts <- list()
  for (g in cohorts) {
    tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g, ]; tr$cohort <- g; tr$treat <- 1L
    parts[[length(parts) + 1L]] <- tr
    if (gate15 || (g - 4 >= 2000)) { ct <- nochange; ct$cohort <- g; ct$treat <- 0L
                                     parts[[length(parts) + 1L]] <- ct }
  }
  d <- do.call(rbind, parts)
  d$etime <- d$year - d$cohort
  d
}

# ---- per-cohort pre-trend test ---------------------------------------------
pretrend_one <- function(d, g, K) {
  s <- d[d$cohort == g & !is.na(d$any.fatalities) & abs(d$etime) <= K, , drop = FALSE]
  nt   <- length(unique(s$agency.id[s$treat == 1]))
  nc   <- length(unique(s$agency.id[s$treat == 0]))
  pre_treat <- sort(unique(s$etime[s$treat == 1 & s$etime < 0]))   # observed treated leads
  testable_leads <- setdiff(pre_treat, -1L)                        # -1 is reference
  base <- data.frame(cohort = g, n_treat = nt, n_ctrl = nc,
                     n_leads = length(testable_leads),
                     wald_stat = NA_real_, pretrend_p = NA_real_,
                     pretrend_p_hetero = NA_real_,
                     reason = NA_character_, stringsAsFactors = FALSE)
  if (nc == 0)                    { base$reason <- "no clean controls (leads collinear)"; return(base) }
  if (length(testable_leads) == 0){ base$reason <- "no observed treated pre-period";     return(base) }
  # cluster-robust (agency.id) is primary; hetero-robust is a robustness column
  # because these cohorts have few treated clusters (often 1-3), where
  # cluster-robust inference degenerates and can overstate the Wald statistic.
  fit    <- tryCatch(
    fixest::feols(any.fatalities ~ i(etime, treat, ref = -1) | agency.id + year,
                  data = s, cluster = ~agency.id, warn = FALSE, notes = FALSE),
    error = function(e) NULL)
  fit_hc <- tryCatch(
    fixest::feols(any.fatalities ~ i(etime, treat, ref = -1) | agency.id + year,
                  data = s, vcov = "hetero", warn = FALSE, notes = FALSE),
    error = function(e) NULL)
  if (is.null(fit)) { base$reason <- "event study failed to fit"; return(base) }
  w  <- tryCatch(fixest::wald(fit,    "etime::-", print = FALSE), error = function(e) NULL)
  wh <- tryCatch(fixest::wald(fit_hc, "etime::-", print = FALSE), error = function(e) NULL)
  if (is.null(w) || is.na(w$p)) { base$reason <- "leads not estimable in fit"; return(base) }
  base$wald_stat <- w$stat; base$pretrend_p <- w$p
  base$pretrend_p_hetero <- if (!is.null(wh)) wh$p else NA_real_
  base$reason <- "tested"
  base
}

run_stack <- function(gate15, tag) {
  d <- build_stack(gate15)
  cohorts <- sort(unique(d$cohort))
  res <- do.call(rbind, lapply(cohorts, function(g) pretrend_one(d, g, K)))
  res$stack <- tag
  res
}

ship  <- run_stack(FALSE, "shipped-11")
rec15 <- run_stack(TRUE,  "reconstructed-15")
res <- rbind(ship, rec15)
res <- res[, c("stack", "cohort", "n_treat", "n_ctrl", "n_leads",
               "wald_stat", "pretrend_p", "pretrend_p_hetero", "reason")]

options(width = 200)
cat("== Per-sub-experiment pre-trend tests (corrected agency x stack FE, +/-", K, "window) ==\n", sep = "")
print(transform(res, wald_stat = round(wald_stat, 2),
                pretrend_p = round(pretrend_p, 4),
                pretrend_p_hetero = round(pretrend_p_hetero, 4)),
      row.names = FALSE)

write.csv(res, file.path(OUT_DIR, "stacked_pretrends.csv"), row.names = FALSE)

# ---- LaTeX (shipped stack; the testable cohorts) ---------------------------
f <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
p <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
mklines <- function(sub, cap, lab) {
  L <- c("\\begin{table}[t]\\centering",
         paste0("\\caption{", cap, "}"),
         paste0("\\label{", lab, "}"),
         "\\begin{tabular}{lrrrrll}",
         "\\toprule",
         "Cohort & Treated & Controls & Leads & Wald $\\chi^2$ & $p$ (cluster) & $p$ (hetero) \\\\",
         "\\midrule")
  for (i in seq_len(nrow(sub))) {
    if (sub$reason[i] == "tested") {
      L <- c(L, sprintf("%d & %d & %d & %d & %s & %s & %s \\\\",
                        sub$cohort[i], sub$n_treat[i], sub$n_ctrl[i], sub$n_leads[i],
                        f(sub$wald_stat[i]), p(sub$pretrend_p[i]), p(sub$pretrend_p_hetero[i])))
    } else {
      L <- c(L, sprintf("%d & %d & %d & %d & \\multicolumn{3}{l}{\\textit{%s}} \\\\",
                        sub$cohort[i], sub$n_treat[i], sub$n_ctrl[i], sub$n_leads[i],
                        sub$reason[i]))
    }
  }
  c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}
lines <- c(
  mklines(ship,
    paste0("Per-sub-experiment pre-trend tests on the shipped stacked dataset, ",
           "corrected interactive fixed effects ($\\text{agency}\\times\\text{stack}+",
           "\\text{year}\\times\\text{cohort}$). Each row fits an event study within one ",
           "cohort's sub-experiment on the $\\pm", K, "$ window and jointly Wald-tests the ",
           "pre-treatment leads; clustered on \\texttt{agency.id}. Control-free early ",
           "cohorts have no estimable leads."),
    "tab:pretrend-shipped"),
  "",
  mklines(rec15,
    paste0("Per-sub-experiment pre-trend tests on the reconstructed 15-cohort stack ",
           "(all change-cohorts $\\ge 2000$ given the 741 clean controls; unweighted ",
           "Table~A5 col.\\ 3 construction), same specification. Attaching controls to ",
           "2001--2003 makes their leads testable; 2000 (no observed pre-period) and 2016 ",
           "(absorbed reversal) remain untestable."),
    "tab:pretrend-recon15"))
writeLines(lines, file.path(OUT_DIR, "stacked_pretrends.tex"))
cat("\nWrote output/stacked_pretrends.csv and .tex\n")
