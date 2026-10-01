# =============================================================================
# task2_task4_from_atoms.R
#
# Pure-atom recombinations that need no raw data or re-fit (reproduce the CSVs
# generated while the container's R/data were wiped):
#   Task 2 (partial): CS `simple` and `group` aggregates + implied cohort weights
#                     from output/atoms_long.csv  -> task2_cs_aggregation_partial.csv
#                     (dynamic/calendar aggregates and all SEs need an aggte re-fit)
#   Task 4.1 (CS LOO): leave-one-cohort-out CS simple aggregate -- CS's never-
#                     treated-control ATT(g,t) are independent across cohorts, so
#                     dropping a cohort = renormalize the simple weights over the
#                     rest. Point estimates exact; SEs need a re-fit.
#                     -> task4_cs_loo.csv
# =============================================================================
OUT <- Sys.getenv("DIDREP_OUT", unset = "output")
al <- read.csv(file.path(OUT, "atoms_long.csv"), stringsAsFactors = FALSE)
cs <- al[al$estimator == "CS" & al$weight > 0, ]
agg <- aggregate(cbind(wx = cs$weight * cs$estimate, w = cs$weight, n = 1) ~ group,
                 data = cs, FUN = sum)
A    <- setNames(agg$wx / agg$w, agg$group)     # within-cohort simple avg = ATT_CS(g)
wsmp <- setNames(agg$w, agg$group)              # simple cohort weight = pg * n_post
npost<- setNames(agg$n, agg$group)
pg   <- wsmp / npost                             # group cohort weight = cohort size

g <- sort(as.integer(names(A))); k <- as.character(g)
Ws <- wsmp[k] / sum(wsmp); Wg <- pg[k] / sum(pg)
simple <- sum(Ws * A[k]); group <- sum(Wg * A[k])

t2 <- data.frame(cohort = c(g, NA, NA), n_post = c(npost[k], NA, NA),
                 A_CS = round(c(A[k], simple, group), 6),
                 w_simple = round(c(Ws, NA, NA), 6),
                 w_group  = round(c(Wg, NA, NA), 6))
t2$cohort <- c(as.character(g), "OVERALL_simple", "OVERALL_group")
write.csv(t2, file.path(OUT, "task2_cs_aggregation_partial.csv"), row.names = FALSE)
cat(sprintf("CS simple = %+.5f | CS group = %+.5f\n", simple, group))
cat(sprintf("2002 weight simple %.4f group %.4f | 2009 weight simple %.4f group %.4f\n",
            Ws["2002"], Wg["2002"], Ws["2009"], Wg["2009"]))

# ---- Task 4.1 CS leave-one-out ----------------------------------------------
full <- sum(wsmp[k] * A[k]) / sum(wsmp[k])
loo <- data.frame(dropped_cohort = "(none-full)", cs_overall = round(full,6),
                  delta = 0, flips_sign = "", stringsAsFactors = FALSE)
for (g0 in k) {
  keep <- setdiff(k, g0); v <- sum(wsmp[keep] * A[keep]) / sum(wsmp[keep])
  loo <- rbind(loo, data.frame(dropped_cohort = g0, cs_overall = round(v,6),
    delta = round(v - full,6), flips_sign = ifelse((v<0)!=(full<0),"YES",""),
    stringsAsFactors = FALSE))
}
write.csv(loo, file.path(OUT, "task4_cs_loo.csv"), row.names = FALSE)
cat(sprintf("\nCS full = %+.5f ; drop 2002 -> %+.5f (sign flip)\n",
            full, loo$cs_overall[loo$dropped_cohort=="2002"]))
cat("Wrote task2_cs_aggregation_partial.csv and task4_cs_loo.csv\n")
