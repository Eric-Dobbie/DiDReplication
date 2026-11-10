# =============================================================================
# cs_loo.R  (Task 4.1, with SEs)
#
# CS leave-one-cohort-out WITH standard errors: for each treated cohort, drop
# that cohort's treated units (keep the never-treated controls), refit att_gt,
# and re-aggregate (simple). Confirms the known result that dropping 2002 flips
# CS's sign, now with multiplier-bootstrap SEs. Point estimates match the atom-
# based LOO (output/task4_cs_loo.csv) because CS's never-treated-control
# ATT(g,t) are independent across treated cohorts. Writes output/cs_loo_se.csv.
# =============================================================================
suppressMessages({library(did)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987, dta$year.changed))

fit_simple <- function(d) {
  set.seed(0)
  mp <- did::att_gt(yname = "any.fatalities", tname = "year", idname = "agency.num",
                    gname = "year.changed", xformla = ~1, control_group = "nevertreated",
                    clustervars = "agency.num", data = d)
  ag <- did::aggte(mp, type = "simple", na.rm = TRUE)
  c(est = as.numeric(ag$overall.att), se = as.numeric(ag$overall.se))
}

cohorts <- sort(unique(dta$year.changed[dta$year.changed >= 2001]))
# count treated units per cohort (for the few-cluster flag)
ntr <- sapply(cohorts, function(g) length(unique(dta$agency.id[dta$year.changed == g])))

base <- fit_simple(dta)
rows <- list(data.frame(dropped = "(none-full)", n_treated = NA,
                        est = base["est"], se = base["se"], delta = 0, flips = ""))
for (i in seq_along(cohorts)) {
  g <- cohorts[i]
  d <- dta[dta$year.changed != g, ]          # drop this cohort's treated units
  v <- fit_simple(d)
  rows[[length(rows)+1]] <- data.frame(
    dropped = as.character(g), n_treated = ntr[i], est = v["est"], se = v["se"],
    delta = v["est"] - base["est"],
    flips = ifelse((v["est"] < 0) != (base["est"] < 0), "YES", ""))
}
res <- do.call(rbind, rows); row.names(res) <- NULL
options(width = 160)
cat("== CS leave-one-cohort-out (simple aggregation, with SE) ==\n")
print(transform(res, est = round(est,4), se = round(se,4), delta = round(delta,4)),
      row.names = FALSE)
write.csv(res, file.path(OUT_DIR, "cs_loo_se.csv"), row.names = FALSE)
cat("\nWrote output/cs_loo_se.csv\n")
