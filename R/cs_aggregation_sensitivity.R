# =============================================================================
# cs_aggregation_sensitivity.R  (Task 2)
#
# How much of the CS result depends on the aggregation scheme rather than the
# estimator? One CS fit (never-treated controls, no covariates, universal base
# period = R5 conditions), aggregated four ways via did::aggte: simple, group,
# dynamic, calendar. Reports point estimate, SE, and the implied per-cohort
# weight under each scheme (the share of the overall estimate attributable to
# each treated cohort). Cell weights W(g,t) are reconstructed from did's formulas
# and verified against aggte's own overall.att before being summed to cohorts.
#
# Writes output/task2_aggregation_sensitivity.csv (one row per scheme: overall,
# se, and the 2002 / 2009 cohort weights) and
# output/task2_cohort_weights.csv (cohort x scheme weight matrix).
# =============================================================================
suppressMessages({library(did)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")

# ---- same recode as run_extraction.R ---------------------------------------
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987, dta$year.changed))

set.seed(0)
mp <- did::att_gt(yname = "any.fatalities", tname = "year", idname = "agency.num",
                  gname = "year.changed", xformla = ~1, control_group = "nevertreated",
                  clustervars = "agency.num", data = dta)

# pg = P(G=g) from the fitted object's one-row-per-unit data
dp <- mp$DIDparams; idata <- dp$data
keepu <- !duplicated(idata[[dp$idname]]); ug <- idata[[dp$gname]][keepu]
att <- mp$att; G <- mp$group; Tt <- mp$t
ok <- !is.na(att)
att <- att[ok]; G <- G[ok]; Tt <- Tt[ok]
groups <- sort(unique(G[G > 0]))
pg <- sapply(groups, function(g) mean(ug == g)); names(pg) <- groups

# ---- reconstruct cell weights W(g,t) for each scheme ------------------------
cellw <- function(scheme) {
  w <- setNames(rep(0, length(att)), seq_along(att))
  if (scheme == "simple") {
    k <- Tt >= G
    w[k] <- pg[as.character(G[k])]; w <- w / sum(w)
  } else if (scheme == "group") {
    for (g in groups) {
      k <- which(G == g & Tt >= g)
      w[k] <- pg[as.character(g)] / length(k)     # uniform within cohort, pg across
    }
    w <- w / sum(w)
  } else if (scheme == "dynamic") {
    es <- Tt - G
    evs <- sort(unique(es[es >= 0]))
    for (e in evs) {
      k <- which(es == e)
      den <- sum(pg[as.character(G[k])])
      w[k] <- (pg[as.character(G[k])] / den) / length(evs)   # average over event times
    }
  } else if (scheme == "calendar") {
    ts <- sort(unique(Tt[Tt >= G]))               # calendar periods with >=1 post cell
    valid <- ts
    for (tt in valid) {
      k <- which(Tt == tt & G <= tt & Tt >= G)
      if (length(k) == 0) next
      den <- sum(pg[as.character(G[k])])
      w[k] <- (pg[as.character(G[k])] / den) / length(valid)
    }
  }
  w
}

schemes <- c("simple", "group", "dynamic", "calendar")
rows <- list(); cohW <- matrix(0, length(groups), length(schemes),
                               dimnames = list(groups, schemes))
for (s in schemes) {
  ag <- did::aggte(mp, type = s, na.rm = TRUE)
  w  <- cellw(s)
  recon <- sum(w * att)
  cw <- tapply(w, G, sum)[as.character(groups)]; cw[is.na(cw)] <- 0
  cohW[, s] <- cw
  rows[[s]] <- data.frame(scheme = s, overall = as.numeric(ag$overall.att),
                          se = as.numeric(ag$overall.se), recon_check = recon,
                          abs_diff = abs(recon - as.numeric(ag$overall.att)),
                          w_2002 = cw["2002"], w_2009 = cw["2009"])
}
res <- do.call(rbind, rows); row.names(res) <- NULL
options(width = 170)
cat("== Task 2: CS aggregation-scheme sensitivity ==\n")
print(transform(res, overall = round(overall,5), se = round(se,5),
                recon_check = round(recon_check,5), abs_diff = signif(abs_diff,2),
                w_2002 = round(w_2002,4), w_2009 = round(w_2009,4)), row.names = FALSE)

write.csv(res, file.path(OUT_DIR, "task2_aggregation_sensitivity.csv"), row.names = FALSE)
cwdf <- data.frame(cohort = groups, round(cohW, 6))
write.csv(cwdf, file.path(OUT_DIR, "task2_cohort_weights.csv"), row.names = FALSE)
cat("\nWrote output/task2_aggregation_sensitivity.csv and task2_cohort_weights.csv\n")
cat("\ncohort-weight matrix:\n"); print(round(cohW, 4))
