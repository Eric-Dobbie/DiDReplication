# =============================================================================
# fe_diagnostic.R
#
# Verifies the fixed-effects specification of the stacked regression on the P&P
# full panel (stacked_fatal.csv == reconstruct_stacked.R output).
#
# The paper's m3/m4 use ADDITIVE FEs `agency.id + year.cohort`. `year.cohort` is
# already year-crossed-with-cohort (a clean stack x time indicator), but
# `agency.id` is the raw agency, pooled across every stack it appears in. Since
# each never-treated control appears in 11 stacks, a shared agency FE leaks
# levels across stacks. The leakage-free spec interacts the unit FE with the
# stack: `agency_stack (= agency.id x cohort) + year.cohort`.
#
# Finding: the POOLED ATT barely moves (~0.005 plain/m3, ~0.002 m4), well within
# the SE -- a non-issue for the paper's headline. But the PER-COHORT effects are
# artificially homogenized to ~-0.11 under the shared FE and reveal their true
# heterogeneity (+0.02 to -0.37) under the interacted FE. run_extraction.R uses
# the interacted (correct) FE for the decomposition; see RECOMBINATION_NOTES.md.
# =============================================================================

suppressMessages({library(lfe); library(dplyr)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))

cat("== Step 1: encodings ==\n")
cat("year.cohort == concat(year, cohort):",
    all(as.numeric(paste0(st$year, st$cohort)) == st$year.cohort), "\n")
cat("distinct year.cohort per (year,cohort) pair (1 => it IS year x cohort):",
    max((st %>% count(year, cohort, year.cohort) %>% count(year, cohort))$n), "\n")
cat("max cohorts a single agency.id appears in:",
    max((st %>% count(agency.id, cohort) %>% count(agency.id))$n), "(=> agency.id is NOT stack-specific)\n")

cat("\n== Step 2: agency reuse across stacks ==\n")
reuse <- st %>% group_by(agency.id) %>%
  summarise(n_stacks = n_distinct(cohort), ever_treated = any(treat == 1), .groups = "drop")
cat("agencies:", nrow(reuse), "| in >1 stack:", sum(reuse$n_stacks > 1),
    sprintf("(%.1f%%)\n", 100 * mean(reuse$n_stacks > 1)))
print(reuse %>% group_by(ever_treated) %>%
        summarise(n = n(), median_stacks = median(n_stacks), max_stacks = max(n_stacks), .groups = "drop"))

cat("\n== Step 3: pooled no.req, original vs interacted unit FE ==\n")
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
covs <- "log.pop + log.med.inc + pct.white + pct.white.officers.imputed"
fml  <- function(fe, w = FALSE) {
  f <- as.formula(sprintf("any.fatalities ~ no.req + %s | %s + year.cohort | 0 | agency.id", covs, fe))
  if (w) felm(f, data = st, weights = st$weights) else felm(f, data = st)
}
fmlp <- function(fe) felm(as.formula(sprintf("any.fatalities ~ no.req | %s + year.cohort | 0 | agency.id", fe)), data = st)
grab <- function(m) c(est = unname(coef(m)["no.req"]), se = unname(m$se["no.req"]))
tab <- rbind(
  `m3   original  agency.id`      = grab(fml("agency.id")),
  `m3   corrected agency x stack` = grab(fml("agency_stack")),
  `plain original  agency.id`     = grab(fmlp("agency.id")),
  `plain corrected agency x stack`= grab(fmlp("agency_stack")),
  `m4   original  agency.id`      = grab(fml("agency.id", w = TRUE)),
  `m4   corrected agency x stack` = grab(fml("agency_stack", w = TRUE)))
print(round(tab, 6))
