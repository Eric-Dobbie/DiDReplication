# =============================================================================
# harmonize_atoms.R
#
# Put all three estimators' disaggregated atoms on a COMMON event-time grid so
# they can be compared and so per-sub-experiment pre-trends can be read off a
# single file. Writes output/atoms_harmonized.csv with one row per
#   (estimator, cohort, event_time, estimate, se, is_pre, identified).
#
#   CS      : ATT(g,t) from atoms_long.csv, event_time = t - g   (t<g = placebo)
#   SA      : CATT(g,e) from atoms_long.csv, event_time = e       (already rel.)
#   stacked : per-sub-experiment (per-cohort) event study under the CORRECTED
#             interactive FE (agency x stack + year x cohort), reference e = -1,
#             restricted to the +/-K balancing window (K = 4, the paper's window).
#             cohort = the sub-experiment; control-free / single-treated cohorts
#             whose leads/lags are collinear are marked identified = FALSE.
#
# The harmonized CS/SA rows carry MARGINAL SEs (the full covariance needed for a
# joint pre-trend Wald test is not storable in a long CSV); stacked_pretrends.R
# recomputes the joint tests by refitting. Point estimates and marginal
# significance are readable directly here.
# =============================================================================

suppressMessages({library(fixest); library(dplyr)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
K        <- 4L                                    # balancing / harmonization window

al <- read.csv(file.path(OUT_DIR, "atoms_long.csv"))

# ---- CS: event_time = t - g ------------------------------------------------
cs <- al[al$estimator == "CS", ]
cs_h <- data.frame(estimator = "CS", cohort = cs$group,
                   event_time = cs$time - cs$group,
                   estimate = cs$estimate, se = cs$se,
                   identified = TRUE, stringsAsFactors = FALSE)

# ---- SA: event_time = e (already relative) ---------------------------------
sa <- al[al$estimator == "SA", ]
sa_h <- data.frame(estimator = "SA", cohort = sa$group,
                   event_time = sa$time,
                   estimate = sa$estimate, se = sa$se,
                   identified = TRUE, stringsAsFactors = FALSE)

# ---- stacked: per-cohort event study on the shipped stack ------------------
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
st$etime <- st$year - st$cohort

stack_event_study <- function(d, g, K) {
  s <- d[d$cohort == g & !is.na(d$any.fatalities) & abs(d$etime) <= K, , drop = FALSE]
  # within one cohort, agency_stack = agency.id and year.cohort = year
  ntre <- length(unique(s$agency.id[s$treat == 1]))
  nctl <- length(unique(s$agency.id[s$treat == 0]))
  # need control variation (nctl > 0) or the treat x etime terms are collinear
  fit <- tryCatch(
    fixest::feols(any.fatalities ~ i(etime, treat, ref = -1) | agency.id + year,
                  data = s, warn = FALSE, notes = FALSE),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- tryCatch(fixest::coeftable(fit), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  rn <- rownames(ct)
  ee <- as.integer(sub("etime::(-?[0-9]+):treat", "\\1", rn))
  keep <- !is.na(ee)
  data.frame(estimator = "stacked", cohort = g, event_time = ee[keep],
             estimate = ct[keep, 1], se = ct[keep, 2],
             identified = nctl > 0 & ntre >= 1, stringsAsFactors = FALSE)
}

cohorts <- sort(unique(st$cohort))
st_h <- do.call(rbind, lapply(cohorts, function(g) stack_event_study(st, g, K)))

# ---- combine + write -------------------------------------------------------
harm <- rbind(cs_h, sa_h, st_h)
harm$is_pre <- harm$event_time < 0
harm <- harm[order(harm$estimator, harm$cohort, harm$event_time), ]
write.csv(harm, file.path(OUT_DIR, "atoms_harmonized.csv"), row.names = FALSE)

cat("Wrote", file.path(OUT_DIR, "atoms_harmonized.csv"),
    "(", nrow(harm), "rows )\n")
cat("  rows by estimator:\n"); print(table(harm$estimator))
cat("  stacked cohorts with an identified event study:",
    paste(sort(unique(st_h$cohort[st_h$identified])), collapse = ", "), "\n")
