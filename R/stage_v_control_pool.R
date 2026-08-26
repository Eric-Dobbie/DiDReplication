# =============================================================================
# stage_v_control_pool.R  --  Stage V: Control-pool verification (GATE)
#
# Verifies the hypothesis that the stacked design's 741 change.type=="No Change"
# control agencies are mostly units permanently in the treated state
# (no.req==1 in every period), while CS/SA use only the 41 with no.req==0
# throughout.
#
# V1  Tabulate the 741 (and all 787) by treatment history; cross-tab against
#     the callaway_replication.R gname recode; confirm the partitions agree.
# V2  Identifying variation: within-unit variation in no.req among controls on
#     the stacked estimation sample; what the always-no.req==1 units contribute.
# V3  Trend divergence 2000-2020: any.fatalities time path for the always-1 vs
#     always-0 non-changers; slopes, equal-slope test, fitted 2000->2020 change.
# V4  Internal inconsistency: agencies treated as already-treated by
#     callaway_replication.R AND as clean controls by the stacked construction.
#
# Base R only (no did/fixest/lfe needed for V1-V4). Deterministic; seed set for
# form's sake. Reads data from DIDREP_DATA, writes to DIDREP_OUT.
# =============================================================================

set.seed(0)
options(stringsAsFactors = FALSE)

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

logf <- file.path(OUT_DIR, "run.log")
lg <- function(...) { cat(..., "\n"); cat(..., "\n", file = logf, append = TRUE) }
cat("", file = logf)  # truncate
lg("================ Stage V run.log ================")
lg("timestamp:", format(Sys.time(), tz = "UTC", usetz = TRUE))
lg("seed: set.seed(0)")
lg("DIDREP_DATA:", DATA_DIR)
lg("DIDREP_OUT :", OUT_DIR)
lg("R version  :", R.version.string)

# rows accumulated for control_pool_verification.csv
V <- list()
add <- function(check, metric, value, note = "") {
  V[[length(V) + 1L]] <<- data.frame(check = check, metric = metric,
                                     value = as.character(value), note = note)
}

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
lg("\n-- loading dta.csv --")
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
lg("dta rows:", nrow(dta), " unique agency.id:", length(unique(dta$agency.id)),
   " years:", paste(range(dta$year), collapse = "-"))
stopifnot(nrow(dta) == 26758L, length(unique(dta$agency.id)) == 787L)

# -----------------------------------------------------------------------------
# Treatment-history classification (per agency, over ALL rows/years)
#   no.req is (for non-changers) constant across the agency's panel.
# -----------------------------------------------------------------------------
hist_class <- function(x) {
  u <- unique(x[!is.na(x)])
  if (length(u) == 0) "allNA"
  else if (length(u) == 1 && u == 0) "always0"
  else if (length(u) == 1 && u == 1) "always1"
  else "mixed"
}
ag_hist <- tapply(dta$no.req, dta$agency.id, hist_class)

# change.type per agency: the two reversibles (memphis, portsmouth) carry both
# "Adopted" and "Dropped" rows; collapse to a set-membership flag for No Change.
ag_nochange <- tapply(dta$change.type, dta$agency.id,
                      function(ct) all(ct == "No Change"))
# agencies whose EVERY row is "No Change" == the 741 non-changers
nochange_ids <- names(ag_nochange)[ag_nochange]

lg("\n== change.type (agency-level) ==")
ct_ag <- tapply(dta$change.type, dta$agency.id,
                function(ct) paste(sort(unique(ct)), collapse = "+"))
lg(paste(capture.output(print(table(ct_ag))), collapse = "\n"))

# =============================================================================
# V1. Tabulate 741 (and 787) by treatment history + cross-tab vs gname recode
# =============================================================================
lg("\n================ V1 ================")

# --- classification of the 741 No Change agencies ---
h741 <- ag_hist[nochange_ids]
t741 <- table(factor(h741, levels = c("always0", "always1", "mixed", "allNA")))
lg("741 No-Change agencies by no.req history:")
lg(paste(capture.output(print(t741)), collapse = "\n"))
add("V1", "n_nochange_agencies", length(nochange_ids), "change.type=='No Change' in every row")
add("V1", "nochange_always0", t741["always0"], "No Change & no.req==0 always")
add("V1", "nochange_always1", t741["always1"], "No Change & no.req==1 always")
add("V1", "nochange_mixed",   t741["mixed"],   "should be 0")

# --- classification of all 787 agencies ---
t787 <- table(factor(ag_hist, levels = c("always0", "always1", "mixed", "allNA")))
lg("\nAll 787 agencies by no.req history:")
lg(paste(capture.output(print(t787)), collapse = "\n"))
add("V1", "n_all_agencies", length(ag_hist), "")
add("V1", "all_always0", t787["always0"], "")
add("V1", "all_always1", t787["always1"], "")
add("V1", "all_mixed",   t787["mixed"],   "the changers (Dropped/Adopted)")

# --- gname recode from callaway_replication.R (line 16-17) ---
# first-treatment collapse (run_extraction.R): reversibles -> first treat year
d <- dta
d$year.changed <- ave(d$year.changed, d$agency.id,
  FUN = function(x){ v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_ })
# recode: NA & no.req==0 -> 0 ; NA & no.req==1 -> 1987 ; else keep (row-wise)
d$gname <- ifelse(is.na(d$year.changed) & d$no.req == 0, 0,
           ifelse(is.na(d$year.changed) & d$no.req == 1, 1987, d$year.changed))
# agency-level gname (constant within agency for non-changers; for changers it
# is the first-treatment year -- constant after the collapse)
ag_gname <- tapply(d$gname, d$agency.id, function(x) {
  u <- unique(x); if (length(u) == 1) u else paste(sort(u), collapse = "/") })

# cross-tab: no.req history  x  gname recode, over the 741 No Change agencies
part_hist <- ifelse(h741 == "always0", "hist:always0",
              ifelse(h741 == "always1", "hist:always1", "hist:other"))
part_gname <- ag_gname[nochange_ids]
part_gname_lab <- ifelse(part_gname == "0", "gname:0",
                  ifelse(part_gname == "1987", "gname:1987", "gname:other"))
xt <- table(part_hist, part_gname_lab)
lg("\nCross-tab (741 No Change): no.req history  x  gname recode")
lg(paste(capture.output(print(xt)), collapse = "\n"))

# confirm partitions agree exactly: always0<->0 , always1<->1987
agree <- (h741 == "always0" & part_gname == "0") |
         (h741 == "always1" & part_gname == "1987")
n_disagree <- sum(!agree)
lg("\nPartitions agree (always0<->gname0, always1<->gname1987)?  disagreements:", n_disagree)
add("V1", "xt_always0_gname0",    xt["hist:always0","gname:0"],    "cross-tab cell")
add("V1", "xt_always1_gname1987", xt["hist:always1","gname:1987"], "cross-tab cell")
add("V1", "partition_disagreements", n_disagree, "agencies not fitting either partition")
if (n_disagree > 0) {
  bad <- nochange_ids[!agree]
  lg("  disagreeing agencies:", paste(bad, collapse = ", "))
  add("V1", "disagreeing_agencies", paste(bad, collapse = ";"), "")
}

# also cross-tab over all 787 for completeness
ag_gname_all <- ag_gname
part_hist_all <- ifelse(ag_hist == "always0", "hist:always0",
                 ifelse(ag_hist == "always1", "hist:always1", "hist:mixed"))
part_gname_all <- ifelse(ag_gname_all == "0", "gname:0",
                  ifelse(ag_gname_all == "1987", "gname:1987", "gname:cohort>=2000"))
xt_all <- table(part_hist_all, part_gname_all)
lg("\nCross-tab (all 787): no.req history  x  gname recode")
lg(paste(capture.output(print(xt_all)), collapse = "\n"))

# split verdict
split_ok <- (as.integer(t741["always0"]) == 41L) && (as.integer(t741["always1"]) == 700L)
add("V1", "split_is_41_700", split_ok, "TRUE => hypothesis premise holds")
lg("\nV1 verdict: split is 41/700 ? ", split_ok)

# =============================================================================
# V2. Identifying variation within the stacked file
# =============================================================================
lg("\n================ V2 ================")
lg("-- loading stacked_fatal.csv --")
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
lg("stacked rows:", nrow(st), " cols:", ncol(st),
   " stacks:", length(unique(st$cohort)))
stopifnot(nrow(st) == 278324L)

st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")

# estimation sample = non-missing outcome
est <- st[!is.na(st$any.fatalities), , drop = FALSE]
lg("estimation-sample rows (non-missing any.fatalities):", nrow(est))
add("V2", "stacked_rows", nrow(st), "")
add("V2", "estimation_rows", nrow(est), "non-missing any.fatalities")

# control rows are treat==0 ; treated are treat==1
lg("estimation rows by treat:")
lg(paste(capture.output(print(table(est$treat))), collapse = "\n"))
add("V2", "est_control_rows", sum(est$treat == 0), "treat==0")
add("V2", "est_treated_rows", sum(est$treat == 1), "treat==1")

# within-unit variation in no.req among CONTROL units.
# Two granularities: (a) per agency.id (pooled across stacks),
#                    (b) per agency_stack (unit x stack, the FE unit).
ctrl <- est[est$treat == 0, , drop = FALSE]
var_by_agency <- tapply(ctrl$no.req, ctrl$agency.id,
                        function(x) length(unique(x[!is.na(x)])) > 1)
var_by_unitstack <- tapply(ctrl$no.req, ctrl$agency_stack,
                        function(x) length(unique(x[!is.na(x)])) > 1)
n_ctrl_agencies <- length(var_by_agency)
n_ctrl_agencies_varying <- sum(var_by_agency)
n_ctrl_unitstacks <- length(var_by_unitstack)
n_ctrl_unitstacks_varying <- sum(var_by_unitstack)
lg("control agencies (unique agency.id):", n_ctrl_agencies,
   " with within-unit no.req variation:", n_ctrl_agencies_varying)
lg("control agency_stacks (unit x stack):", n_ctrl_unitstacks,
   " with within-unit no.req variation:", n_ctrl_unitstacks_varying)
add("V2", "control_agencies", n_ctrl_agencies, "unique agency.id, treat==0, est sample")
add("V2", "control_agencies_with_noreq_variation", n_ctrl_agencies_varying, "")
add("V2", "control_unitstacks", n_ctrl_unitstacks, "agency x stack, treat==0")
add("V2", "control_unitstacks_with_noreq_variation", n_ctrl_unitstacks_varying, "")

# confirm the always-no.req==1 controls contribute NO no.req variation
ctrl_agencies <- names(var_by_agency)
ctrl_hist <- ag_hist[ctrl_agencies]
always1_ctrl <- ctrl_agencies[ctrl_hist == "always1"]
always0_ctrl <- ctrl_agencies[ctrl_hist == "always0"]
n_always1_varying <- sum(var_by_agency[always1_ctrl])
lg("always-no.req==1 control agencies:", length(always1_ctrl),
   " of which vary no.req within unit:", n_always1_varying, "(expect 0)")
add("V2", "always1_control_agencies", length(always1_ctrl), "")
add("V2", "always1_controls_with_noreq_variation", n_always1_varying, "expect 0")
add("V2", "always0_control_agencies", length(always0_ctrl), "")

# What the always-1 units DO contribute:
#  - agency_stack unit fixed effects (they anchor unit means within each stack)
#  - year.cohort (stack x time) fixed effects (they populate every cohort-year)
#  - their SHARE of estimation rows used to fit year.cohort
rows_always1 <- sum(est$agency.id %in% always1_ctrl)
rows_always0 <- sum(est$agency.id %in% always0_ctrl)
share_always1 <- rows_always1 / nrow(est)
n_yearcohort <- length(unique(est$year.cohort))
# how many distinct year.cohort levels have an always-1 control present
yc_touched_by_always1 <- length(unique(est$year.cohort[est$agency.id %in% always1_ctrl]))
lg("estimation rows from always-1 controls:", rows_always1,
   sprintf(" (%.1f%% of estimation sample)", 100*share_always1))
lg("distinct year.cohort levels:", n_yearcohort,
   " touched by always-1 controls:", yc_touched_by_always1)
add("V2", "est_rows_from_always1_controls", rows_always1, "")
add("V2", "share_est_rows_always1_controls", round(share_always1, 4),
    "share of rows fitting year.cohort FE")
add("V2", "n_yearcohort_levels", n_yearcohort, "")
add("V2", "yearcohort_levels_touched_by_always1", yc_touched_by_always1, "")

# =============================================================================
# V3. Trend divergence 2000-2020
# =============================================================================
lg("\n================ V3 ================")
# use dta (one row per agency-year); outcome observed 2000-2020
d3 <- dta[dta$year >= 2000 & dta$year <= 2020 & !is.na(dta$any.fatalities), ]
d3 <- d3[d3$agency.id %in% nochange_ids, ]          # non-changers only
d3$grp <- ifelse(d3$agency.id %in% always1_ctrl, "always1",
          ifelse(d3$agency.id %in% always0_ctrl, "always0", NA))
d3 <- d3[!is.na(d3$grp), ]
lg("V3 sample rows:", nrow(d3),
   " agencies:", length(unique(d3$agency.id)),
   " always1:", length(unique(d3$agency.id[d3$grp=="always1"])),
   " always0:", length(unique(d3$agency.id[d3$grp=="always0"])))
add("V3", "n_obs", nrow(d3), "2000-2020, non-changers, outcome observed")
add("V3", "n_units_always1", length(unique(d3$agency.id[d3$grp=="always1"])), "")
add("V3", "n_units_always0", length(unique(d3$agency.id[d3$grp=="always0"])), "")

# annual means by group
am <- aggregate(any.fatalities ~ year + grp, data = d3, FUN = mean)
am_w <- reshape(am, idvar = "year", timevar = "grp", direction = "wide")
names(am_w) <- sub("any.fatalities.", "mean_", names(am_w))
am_w <- am_w[order(am_w$year), ]
lg("\nAnnual means (any.fatalities) by group:")
lg(paste(capture.output(print(am_w, row.names = FALSE)), collapse = "\n"))
write.csv(am_w, file.path(OUT_DIR, "v3_annual_means.csv"), row.names = FALSE)

# linear trend per group (cluster-free OLS; SE from lm)
fit_a <- lm(any.fatalities ~ year, data = d3[d3$grp=="always1", ])
fit_b <- lm(any.fatalities ~ year, data = d3[d3$grp=="always0", ])
sa <- summary(fit_a)$coefficients["year", ]
sb <- summary(fit_b)$coefficients["year", ]
lg(sprintf("\nalways1 slope: %.6f (SE %.6f), p=%.4g", sa[1], sa[2], sa[4]))
lg(sprintf("always0 slope: %.6f (SE %.6f), p=%.4g", sb[1], sb[2], sb[4]))
add("V3", "slope_always1", round(sa[1],6), sprintf("SE=%.6f", sa[2]))
add("V3", "slope_always0", round(sb[1],6), sprintf("SE=%.6f", sb[2]))

# test of equal slopes: interaction grp:year
d3$grp <- factor(d3$grp, levels = c("always0","always1"))
fit_i <- lm(any.fatalities ~ year * grp, data = d3)
si <- summary(fit_i)$coefficients
int_row <- si["year:grpalways1", ]
lg(sprintf("equal-slope test (year:grpalways1): diff=%.6f SE=%.6f t=%.3f p=%.4g",
           int_row[1], int_row[2], int_row[3], int_row[4]))
add("V3", "slope_diff_always1_minus_always0", round(int_row[1],6),
    sprintf("SE=%.6f", int_row[2]))
add("V3", "equal_slope_test_t", round(int_row[3],4), "")
add("V3", "equal_slope_test_p", signif(int_row[4],4), "H0: equal slopes")

# fitted 2000->2020 change per group
fitted_change <- function(fit) as.numeric(coef(fit)["year"]) * (2020 - 2000)
fc_a <- fitted_change(fit_a); fc_b <- fitted_change(fit_b)
lg(sprintf("fitted 2000->2020 change: always1=%.5f  always0=%.5f  diff=%.5f",
           fc_a, fc_b, fc_a - fc_b))
add("V3", "fitted_change_2000_2020_always1", round(fc_a,5), "")
add("V3", "fitted_change_2000_2020_always0", round(fc_b,5), "")
add("V3", "fitted_change_diff", round(fc_a - fc_b,5), "always1 minus always0")

# =============================================================================
# V4. Internal inconsistency check
# =============================================================================
lg("\n================ V4 ================")
# (a) already-treated in CS: agencies recoded gname==1987 by callaway_replication.R
cs_alreadytreated <- names(ag_gname)[ag_gname == "1987"]
lg("agencies recoded gname==1987 (CS already-treated, dropped from control pool):",
   length(cs_alreadytreated))
# (b) clean controls in the stacked file: treat==0 agencies (change.type==No Change)
stacked_controls <- unique(st$agency.id[st$treat == 0])
lg("distinct agencies appearing as stacked clean controls (treat==0):",
   length(stacked_controls))
# intersection: receive BOTH treatments
both <- intersect(cs_alreadytreated, stacked_controls)
lg("agencies that are BOTH CS-already-treated AND stacked clean controls:",
   length(both))
add("V4", "cs_already_treated_gname1987", length(cs_alreadytreated),
    "callaway_replication.R:17 recode; att_gt drops as already-treated")
add("V4", "stacked_clean_controls", length(stacked_controls),
    "treat==0 in stacked_fatal.csv")
add("V4", "agencies_receiving_both_treatments", length(both),
    "already-treated by CS AND clean control by stack")
# sanity: these should be exactly the 700 always-no.req==1 non-changers
same_as_always1 <- setequal(both, always1_ctrl)
lg("both-set == the 700 always-no.req==1 non-changers ?", same_as_always1)
add("V4", "both_equals_always1_set", same_as_always1, "")
add("V4", "cite_callaway_line", "callaway_replication.R:17 (recode) + :52-58 (att_gt)",
    "NA & no.req==1 -> 1987; att_gt drops gname==1987")
add("V4", "cite_stacked_construction",
    "reconstruct_stacked.R: nochange<-dta[change.type=='No Change',]; treat<-0",
    "authors' create_stacked_data.R not in zip; reconstruct_stacked.R replicates & matches")

# =============================================================================
# Write outputs
# =============================================================================
vdf <- do.call(rbind, V)
write.csv(vdf, file.path(OUT_DIR, "control_pool_verification.csv"), row.names = FALSE)
lg("\nWrote:", file.path(OUT_DIR, "control_pool_verification.csv"), "(", nrow(vdf), "rows )")

lg("\n================ sessionInfo() ================")
lg(paste(capture.output(print(sessionInfo())), collapse = "\n"))

cat("\nDONE Stage V\n")
