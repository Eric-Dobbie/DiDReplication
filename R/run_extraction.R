# =============================================================================
# run_extraction.R  --  driver for the atom/weight extraction
#
# Loads the Payson & Parinandi (2024) replication data, supplies the ACTUAL
# column names (hardcoded here, not in the functions), calls the three
# data-agnostic extractors in extract_atoms.R, runs the recombination check,
# and writes two CSVs:
#     output/atoms_long.csv          -- all atoms + weights, all three estimators
#     output/recombination_check.csv -- reconstructed vs package pooled estimate
#
# Outcome of interest: any.fatalities (probability of a fatal civilian
# encounter). This is the fatal-encounters analysis (Table 3), NOT the racial
# diversity analysis.
#
# Data location: set DATA_DIR to the folder holding dta.csv and
# stacked_fatal.csv (defaults to ./data, override with the DIDREP_DATA env var).
# =============================================================================

# locate extract_atoms.R next to this script (works under Rscript and source())
.this_file <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else sys.frame(1)$ofile
}, error = function(e) NULL)
.here <- if (!is.null(.this_file)) dirname(.this_file) else "R"
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- P&P column names (the ONLY place these are hardcoded) ------------------
Y      <- "any.fatalities"                # outcome
TIME   <- "year"                          # calendar time
UNIT   <- "agency.id"                     # unit id (string)
UNITN  <- "agency.num"                    # numeric unit id (needed by did)
COHORT <- "year.changed"                  # first-treatment (drop) year; NA = never
# stacked file:
S_TREAT <- "no.req"                       # treatment indicator (requirement absent)
S_TIME  <- "year.cohort"                  # stack-specific time FE (year x cohort)
S_STACK <- "cohort"                       # sub-experiment / stack id

# =============================================================================
# Load + diagnostics
# =============================================================================
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta[[UNITN]] <- as.numeric(factor(dta[[UNIT]]))

# Enforce irreversible treatment: two agencies (memphis TN, portsmouth NH) both
# ADOPTED and later DROPPED a requirement, so year.changed is not unit-constant.
# The paper flags exactly these two (callaway_replication.R drops them), and
# did >= 2.2 now errors on reversible gname. Collapse each unit to its FIRST
# treatment year (standard staggered-adoption timing) -- this reproduces the
# behavior of did 2.1.2 for CS and applies the same definition to SA, so the two
# estimators use a consistent cohort assignment.
dta[[COHORT]] <- ave(dta[[COHORT]], dta[[UNIT]],
                     FUN = function(x) { v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_ })

# P&P coding of the cohort variable for CS/SA (from callaway_replication.R):
#   never-treated (always has a requirement)  -> 0   (clean controls)
#   always-no-requirement (treated pre-sample)-> 1987 (dropped by att_gt: no pre)
dta[[COHORT]] <- ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]] == 0, 0,
                 ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]] == 1, 1987,
                        dta[[COHORT]]))

cat("== dta.csv ==  rows:", nrow(dta),
    " units:", length(unique(dta[[UNIT]])),
    " years:", paste(range(dta[[TIME]]), collapse = "-"), "\n")
cat("   outcome non-missing years:",
    paste(range(dta[[TIME]][!is.na(dta[[Y]])]), collapse = "-"), "\n")

# =============================================================================
# Step 2 -- CS atoms   (att_gt auto-drops already-treated + missing-outcome rows)
# =============================================================================
cs <- extract_cs_atoms(dta, yname = Y, tname = TIME, idname = UNITN,
                       gname = COHORT, xformla = ~1,
                       control_group = "nevertreated", clustervars = UNITN)

# =============================================================================
# Step 3 -- SA atoms
# Comparable sample: never-treated controls + cohorts estimable within the
# outcome window (first-treated >= 2001); drop units treated before the sample
# (1987..2000), exactly the set att_gt drops as "already treated in first period".
# =============================================================================
d_sa <- dta[!is.na(dta[[Y]]) &
              (dta[[COHORT]] == 0 | dta[[COHORT]] >= 2001), , drop = FALSE]
sa <- extract_sa_atoms(d_sa, yname = Y, idname = UNITN, tname = TIME,
                       gname = COHORT, never_value = 0, ref_c = 10000)

# =============================================================================
# Step 4 -- Stacked atoms  (their pre-built clean-control stack)
# =============================================================================
stacked <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))   # P&P full panel = reconstruct_stacked.R output
cat("== stacked_fatal.csv ==  rows:", nrow(stacked),
    " stacks:", length(unique(stacked[[S_STACK]])), "\n")

# LEAKAGE-FREE (correct) stacked FE: interact the unit FE with the stack, so a
# control agency (which appears in 11 stacks) is demeaned separately within each
# stack -- the fully-interacted Cengiz FE (agency x stack + year x cohort). The
# paper's m3/m4 instead use a shared agency.id FE, which leaks levels across
# stacks and artificially homogenizes the per-cohort effects (see fe_diagnostic.R
# and RECOMBINATION_NOTES.md). Clustering stays on the real agency.id.
STACKUNIT <- "agency_stack"
stacked[[STACKUNIT]] <- paste(stacked[[UNIT]], stacked[[S_STACK]], sep = "__")

st <- extract_stacked_atoms(stacked, yname = Y, idname = STACKUNIT,
                            stacktime = S_TIME, stackid = S_STACK,
                            treatname = S_TREAT, cluster = UNIT,
                            resid_on = "estimation")

# also compute the naive "residualize on all rows" version for the diagnostic
st_naive <- extract_stacked_atoms(stacked, yname = Y, idname = STACKUNIT,
                                  stacktime = S_TIME, stackid = S_STACK,
                                  treatname = S_TREAT, cluster = UNIT,
                                  resid_on = "all")

# =============================================================================
# Step 5 -- Recombination check
# =============================================================================
recomb_row <- function(df, label = attr(df, "estimator_label")) {
  recon  <- sum(df$weight * df$estimate)
  pooled <- attr(df, "pooled_estimate")
  data.frame(
    estimator            = df$estimator[1],
    pooled_label         = attr(df, "pooled_label"),
    reconstructed_pooled = recon,
    package_pooled       = pooled,
    abs_diff             = abs(recon - pooled),
    rel_diff             = abs(recon - pooled) / abs(pooled),
    stringsAsFactors = FALSE
  )
}

check <- rbind(recomb_row(cs), recomb_row(sa), recomb_row(st))

# extra diagnostic line: naive stacked residualization (all rows)
naive_recon <- sum(st_naive$weight * st_naive$estimate)
check_naive <- data.frame(
  estimator            = "stacked_naive_residALLrows",
  pooled_label         = attr(st, "pooled_label"),
  reconstructed_pooled = naive_recon,
  package_pooled       = attr(st, "pooled_estimate"),
  abs_diff             = abs(naive_recon - attr(st, "pooled_estimate")),
  rel_diff             = abs(naive_recon - attr(st, "pooled_estimate")) /
                           abs(attr(st, "pooled_estimate")),
  stringsAsFactors = FALSE)
check_full <- rbind(check, check_naive)

cat("\n================ RECOMBINATION CHECK ================\n")
print(check_full, row.names = FALSE, digits = 8)

# =============================================================================
# Step 6 -- Output
# =============================================================================
atoms_long <- rbind(cs, sa, st)
write.csv(atoms_long, file.path(OUT_DIR, "atoms_long.csv"), row.names = FALSE)
write.csv(check_full, file.path(OUT_DIR, "recombination_check.csv"), row.names = FALSE)
cat("\nWrote:", file.path(OUT_DIR, "atoms_long.csv"), "(", nrow(atoms_long), "atoms )\n")
cat("Wrote:", file.path(OUT_DIR, "recombination_check.csv"), "\n")
