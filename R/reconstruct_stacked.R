# =============================================================================
# reconstruct_stacked.R
#
# Rebuild stacked_fatal.csv from dta.csv "from first principles" and check
# whether it matches the provided file. We do NOT have the authors' stack-
# construction code, so this tests how well the construction rule is recovered.
#
# Recovered rule (inferred from the provided file; see RECOMBINATION_NOTES.md):
#   * Stacks are indexed by treatment cohort g = a policy CHANGE year.
#   * Cohorts included: change years with g >= 2000 (the outcome, any.fatalities,
#     is only observed from 2000, so an earlier change has no observable
#     treatment period).
#   * Treated rows of stack g = the dta rows with year.changed == g. Because
#     year.changed is stored row-wise, the two REVERSIBLE agencies (memphis TN,
#     portsmouth NH, which both adopted and dropped) are automatically split at
#     their change point: the pre-change segment lands in the adopt-cohort, the
#     post-change segment in the drop-cohort.
#   * Clean controls = the 741 never-changing ("No Change") agencies, full
#     panels, attached ONLY when a full 4-period pre-window is observable, i.e.
#     g - 4 >= 2000  (=> cohorts 2004..2020). This is the "implicit window":
#     +/-4 (4 pre, 4 post). Note it gates CONTROL ELIGIBILITY and defines
#     scaled.year / the balancing sample -- it does NOT trim rows: every unit
#     keeps its full 1987-2020 panel.
#   * Derived columns: treat (1 treated / 0 control), scaled.year = year - g on
#     [-4, 4] else NA, year.cohort = paste0(year, g).
#
# The one column this cannot reproduce is `weights` (entropy-balancing weights
# for the controls); that requires the authors' ebal() specification.
# =============================================================================

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
dta  <- read.csv(file.path(DATA_DIR, "dta.csv"))
prov <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))

WINDOW    <- 4                         # pre/post periods of the implicit window
FIRST_OUT <- 2000                      # first year the outcome is observed

cohorts  <- sort(unique(dta$year.changed[!is.na(dta$year.changed) &
                                         dta$year.changed >= FIRST_OUT]))
nochange <- dta[dta$change.type == "No Change", ]     # 741 never-changers

parts <- list()
for (g in cohorts) {
  tr <- dta[!is.na(dta$year.changed) & dta$year.changed == g, ]   # treated segment
  tr$cohort <- g; tr$treat <- 1L
  parts[[length(parts) + 1L]] <- tr
  if (g - WINDOW >= FIRST_OUT) {                        # control eligibility
    ct <- nochange
    ct$cohort <- g; ct$treat <- 0L
    parts[[length(parts) + 1L]] <- ct
  }
}
recon <- do.call(rbind, parts)
recon$scaled.year <- ifelse(abs(recon$year - recon$cohort) <= WINDOW,
                            recon$year - recon$cohort, NA_real_)
recon$year.cohort <- as.numeric(paste0(recon$year, recon$cohort))

# ---- compare to the provided file ------------------------------------------
key <- c("cohort", "agency.id", "year")
rk <- do.call(paste, c(recon[key], sep = "|"))
pk <- do.call(paste, c(prov[key], sep = "|"))

cat("reconstructed rows:", nrow(recon), " | provided rows:", nrow(prov), "\n")
cat("cohorts identical :", identical(sort(unique(recon$cohort)), sort(unique(prov$cohort))), "\n")
cat("row-set identical :", setequal(rk, pk),
    " (recon-only:", length(setdiff(rk, pk)), ", provided-only:", length(setdiff(pk, rk)), ")\n")

# align on the shared key and compare the constructed / carried columns
mg <- merge(recon[c(key, "treat", "no.req", "any.fatalities", "scaled.year", "year.cohort")],
            prov[c(key, "treat", "no.req", "any.fatalities", "scaled.year", "year.cohort")],
            by = key, suffixes = c(".re", ".pr"))
eqv <- function(a, b) sum(!((a == b) | (is.na(a) & is.na(b))))
cat("\nmismatches on shared rows (should all be 0):\n")
for (v in c("treat", "no.req", "any.fatalities", "scaled.year", "year.cohort"))
  cat(sprintf("  %-15s %d\n", v, eqv(mg[[paste0(v, ".re")]], mg[[paste0(v, ".pr")]])))

cat("\nNOTE: the `weights` (entropy-balancing) column is NOT reconstructed here",
    "\n      -- it requires the authors' ebal() specification. Treated weights are 1;",
    "\n      control weights sum to ~70 per cohort (their normalization).\n")
