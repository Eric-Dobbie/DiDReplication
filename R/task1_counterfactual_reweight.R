# =============================================================================
# task1_counterfactual_reweight.R
#
# Counterfactual reweighting: is the CS-vs-stacked disagreement driven by the
# ATOMS (per-cohort effects) or by the aggregation WEIGHTS? Pure recombination of
# already-estimated, committed atoms -- NO raw data or re-fit required:
#   * CS cohort atoms A_CS(g)      = within-cohort simple average of ATT(g,t)
#                                    over post cells, from output/atoms_long.csv
#   * CS simple cohort weights     = sum of the simple per-cell weights over g
#   * stacked atoms beta_s + w_s   = from output/fwl_decomp_unweighted.csv
#     (corrected agency x stack FE, no covariates)
# Restrict to the common cohort set = the cohorts the stacked spec identifies
# (w_s > 0): 2002, 2004, 2008, 2009, 2012, 2013, 2014, 2017, 2018, 2019, 2020.
#
# Writes: output/task1_cohort_inputs.csv (hand-checkable inputs),
#         output/task1_reweight_2x2.csv  (2x2 table + channel decomposition).
# =============================================================================
OUT <- Sys.getenv("DIDREP_OUT", unset = "output")
al <- read.csv(file.path(OUT, "atoms_long.csv"), stringsAsFactors = FALSE)
fw <- read.csv(file.path(OUT, "fwl_decomp_unweighted.csv"), stringsAsFactors = FALSE)

# ---- CS cohort atoms + simple cohort weights (post cells carry weight>0) -----
cs <- al[al$estimator == "CS" & al$weight > 0, ]
agg <- aggregate(cbind(wx = cs$weight * cs$estimate, w = cs$weight,
                       n = 1) ~ group, data = cs, FUN = sum)
A_CS   <- setNames(agg$wx / agg$w, agg$group)     # within-cohort simple average
wCS_raw<- setNames(agg$w, agg$group)              # cohort total simple weight
npost  <- setNames(agg$n, agg$group)

# ---- stacked cohort atoms + FWL weights (identified cohorts only) ------------
fi <- fw[fw$w_s > 0 & !is.na(fw$beta_s), ]
A_ST   <- setNames(fi$beta_s, fi$stack)
wST_raw<- setNames(fi$w_s,   fi$stack)

common <- sort(intersect(names(A_ST), names(A_CS)))
norm   <- function(v) v[common] / sum(v[common])
WCS <- norm(wCS_raw); WST <- norm(wST_raw)
combo <- function(A, W) sum(W[common] * A[common])

csR  <- combo(A_CS, WCS); stR  <- combo(A_ST, WST)
csST <- combo(A_CS, WST); stCS <- combo(A_ST, WCS)
csFULL <- sum(wCS_raw * A_CS[names(wCS_raw)]) / sum(wCS_raw)   # all CS cohorts
stFULL <- sum(wST_raw * A_ST[names(wST_raw)]) / sum(wST_raw)   # = restricted

# ---- gap decomposition: stacked_restr - CS_restr ----------------------------
tot   <- stR - csR
wch   <- sum((WST[common] - WCS[common]) * A_CS[common])          # weight channel
ach   <- sum(WCS[common] * (A_ST[common] - A_CS[common]))         # atom channel
inter <- sum((WST[common] - WCS[common]) * (A_ST[common] - A_CS[common]))

cat(sprintf("common cohorts (n=%d): %s\n", length(common), paste(common, collapse=", ")))
cat("\n== 2x2 (restricted to common set) ==\n")
cat(sprintf("  CS atoms , CS simple  : %+.5f   (unrestricted CS = %+.5f)\n", csR, csFULL))
cat(sprintf("  CS atoms , stacked FWL: %+.5f\n", csST))
cat(sprintf("  ST atoms , CS simple  : %+.5f\n", stCS))
cat(sprintf("  ST atoms , stacked FWL: %+.5f   (unrestricted stacked = %+.5f)\n", stR, stFULL))
cat(sprintf("\n== gap decomposition (total %+.5f) ==\n", tot))
cat(sprintf("  weight channel : %+.5f\n  atom channel   : %+.5f\n  interaction    : %+.5f\n  sum check      : %+.5f\n",
            wch, ach, inter, wch + ach + inter))

# ---- write outputs ----------------------------------------------------------
ci <- data.frame(cohort = common,
                 A_CS = A_CS[common], A_stacked = A_ST[common],
                 wCS_raw = wCS_raw[common], wStacked_raw = wST_raw[common],
                 WCS_restr = WCS[common], WStacked_restr = WST[common],
                 n_post_CS = npost[common], row.names = NULL)
write.csv(ci, file.path(OUT, "task1_cohort_inputs.csv"), row.names = FALSE)

res <- data.frame(
  quantity = c("CS_atoms__CS_simple","CS_atoms__stacked_FWL","stacked_atoms__CS_simple",
               "stacked_atoms__stacked_FWL","CS_simple_unrestricted","stacked_unrestricted",
               "gap_total","channel_weight","channel_atom","channel_interaction"),
  value = round(c(csR, csST, stCS, stR, csFULL, stFULL, tot, wch, ach, inter), 6))
write.csv(res, file.path(OUT, "task1_reweight_2x2.csv"), row.names = FALSE)
cat("\nWrote output/task1_cohort_inputs.csv and task1_reweight_2x2.csv\n")
