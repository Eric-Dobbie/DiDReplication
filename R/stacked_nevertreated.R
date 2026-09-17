# =============================================================================
# stacked_nevertreated.R
#
# Stacked regression restricting the control pool to the STRICTLY never-treated
# agencies only -- the 41 that always have a residency requirement (no.req == 0
# in every row) -- dropping the 700 "always-treated" non-changers (no.req == 1
# always) that the shipped stack also carries as controls (see stack_composition.R
# and assumptions §6/§12). Corrected interactive FE (agency x stack + year x
# cohort), clustered on agency.id. Reports the never-treated-only estimate next to
# the full No-Change-pool baseline, plain and with the paper's covariates.
#
# The entropy-balancing weights (m4/Table 3) were constructed for the full 741-
# control pool, so they do NOT validly rebalance the 41-control pool; only the
# unweighted plain and covariate (m3) specs are reported here.
# =============================================================================
suppressMessages(library(lfe))
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")
covs <- "log.pop + log.med.inc + pct.white + pct.white.officers.imputed"

# never-treated controls = control rows (treat==0) with no.req==0 in ALL rows
ctl    <- st[st$treat == 0, ]
nt_ids <- names(which(tapply(ctl$no.req, ctl$agency.id, function(v) all(v == 0))))
nt     <- st[st$treat == 1 | st$agency.id %in% nt_ids, ]   # treated + never-treated ctrl
cat("never-treated control agencies:", length(nt_ids), "\n")

fit <- function(d, cov) {
  f <- if (cov) as.formula(paste0("any.fatalities ~ no.req + ", covs,
                                  " | agency_stack + year.cohort | 0 | agency.id"))
       else      any.fatalities ~ no.req | agency_stack + year.cohort | 0 | agency.id
  m <- felm(f, data = d)
  c(est = unname(coef(m)["no.req"]), se = unname(m$se["no.req"]), N = m$N)
}

res <- rbind(
  data.frame(control_pool = "never-treated only (41)", spec = "plain", t(fit(nt, FALSE))),
  data.frame(control_pool = "never-treated only (41)", spec = "m3",    t(fit(nt, TRUE))),
  data.frame(control_pool = "full No-Change (741)",    spec = "plain", t(fit(st, FALSE))),
  data.frame(control_pool = "full No-Change (741)",    spec = "m3",    t(fit(st, TRUE))))
row.names(res) <- NULL
res$est <- round(res$est, 6); res$se <- round(res$se, 6)
print(res, row.names = FALSE)
write.csv(res, file.path(OUT_DIR, "stacked_nevertreated.csv"), row.names = FALSE)
cat("\nWrote output/stacked_nevertreated.csv\n")
