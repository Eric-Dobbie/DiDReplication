# =============================================================================
# fwl_decomp.R  --  FWL variance-share decomposition of the R2 stacked spec
#
# R2 spec: unit-by-stack (agency.id x cohort) + time-by-stack (year.cohort) FE.
# Verifies the exact TWFE identity  pooled_beta = sum_s w_s * beta_s , where
#   v_hat = residual of no.req partialled on the two stack-interacted FE sets
#           (explicit regression -- the stacks are unbalanced, so the two-way
#            demeaning closed form does not apply),
#   V_s   = sum_i v_hat_i^2 within stack s,   w_s = V_s / sum_s V_s,
#   beta_s= per-stack regression of the outcome on no.req with unit + time FE.
#
# The identity is EXACT here because both FE dimensions are stack-interacted, so
# the pooled normal equations are block-diagonal across stacks (up to the single
# shared no.req coefficient). Steps 7-8 stress it: entropy-balancing weights
# (identity preserved) and additive (non-interacted) unit FE (identity broken,
# because the 743/776 agencies shared across stacks tie the blocks together).
#
# Sample: the pooled beta is fit listwise on rows with an observed outcome, so
# steps 1-7 use that estimation sample. An aside reports the all-rows gap.
#
# This script REPORTS faithfully; it does not reconcile or "fix" discrepancies.
# Data via DIDREP_DATA (default ./data); CSV summaries to DIDREP_OUT (./output).
# =============================================================================
suppressMessages({library(lfe)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))

Y <- "any.fatalities"; TR <- "no.req"; STK <- "cohort"
st$agency_stack <- paste(st$agency.id, st$cohort, sep = "__")   # unit x stack
# time x stack FE = year.cohort (already = paste0(year, cohort), verified)

# ---- sample: the pooled regression is listwise on observed outcome ----------
full <- st
est  <- st[!is.na(st[[Y]]), , drop = FALSE]
cat(sprintf("Rows: full=%d  estimation(observed outcome)=%d\n", nrow(full), nrow(est)))
cat("All steps 1-7 use the ESTIMATION sample (the sample the pooled beta is fit on).\n\n")

fpr <- function(x) formatC(x, format = "e", digits = 15)

# =============================================================================
# steps 1-6, given a data frame and optional weights
# =============================================================================
decomp <- function(d, wname = NULL) {
  w <- if (is.null(wname)) rep(1, nrow(d)) else d[[wname]]

  # -- step 1: residualize no.req on unit-by-stack + time-by-stack (explicit reg)
  if (is.null(wname)) {
    r1 <- felm(no.req ~ 1 | agency_stack + year.cohort, data = d)
  } else {
    r1 <- felm(no.req ~ 1 | agency_stack + year.cohort, data = d, weights = w)
  }
  d$vhat <- as.numeric(r1$residuals[, 1])

  # -- step 2: within-stack mean of vhat (weighted if applicable)
  stacks <- sort(unique(d[[STK]]))
  wmean <- sapply(stacks, function(s) {
    ix <- d[[STK]] == s; sum(w[ix] * d$vhat[ix]) / sum(w[ix]) })
  maxabs_mean <- max(abs(wmean))

  # -- step 3: V_s, w_s, shares, N
  Vs <- sapply(stacks, function(s) { ix <- d[[STK]] == s; sum(w[ix] * d$vhat[ix]^2) })
  V  <- sum(Vs); ws <- Vs / V
  trsh <- sapply(stacks, function(s) { ix <- d[[STK]] == s; sum(w[ix]*d$treat[ix])/sum(w[ix]) })
  posh <- sapply(stacks, function(s) { ix <- d[[STK]] == s; sum(w[ix]*(d$year[ix] >= s))/sum(w[ix]) })
  Ns   <- sapply(stacks, function(s) sum(d[[STK]] == s))
  tab <- data.frame(stack = stacks, V_s = Vs, w_s = ws,
                    treated_share = trsh, post_share = posh, N = Ns, row.names = NULL)

  # -- step 4: beta_s (per stack: y ~ no.req | unit + time), only that stack
  betas <- sapply(stacks, function(s) {
    ds <- d[d[[STK]] == s, , drop = FALSE]
    m <- tryCatch({
      if (is.null(wname)) felm(any.fatalities ~ no.req | agency.id + year, data = ds)
      else felm(any.fatalities ~ no.req | agency.id + year, data = ds, weights = ds[[wname]])
    }, error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    b <- coef(m)["no.req"]; if (is.na(b)) NA_real_ else unname(b)
  })
  tab$beta_s <- betas

  # -- step 5: pooled beta (unit-by-stack + time-by-stack)
  if (is.null(wname)) {
    mp <- felm(any.fatalities ~ no.req | agency_stack + year.cohort, data = d)
  } else {
    mp <- felm(any.fatalities ~ no.req | agency_stack + year.cohort, data = d, weights = w)
  }
  pooled <- unname(coef(mp)["no.req"])

  # -- step 6: identity sum_s w_s beta_s (V_s=0 stacks contribute 0)
  contrib <- ifelse(tab$V_s == 0 | is.na(tab$beta_s), 0, tab$w_s * tab$beta_s)
  recon <- sum(contrib)

  list(tab = tab, maxabs_mean = maxabs_mean, pooled = pooled, recon = recon,
       zero_stacks = stacks[Vs == 0], vhat = d$vhat)
}

report <- function(res, header) {
  cat("========================================================\n", header, "\n",
      "========================================================\n", sep = "")
  cat(sprintf("[step 2] max|within-stack mean of vhat| = %s\n\n", fpr(res$maxabs_mean)))
  cat("[step 3-4] per-stack table:\n")
  pt <- res$tab
  pt$V_s <- signif(pt$V_s, 8); pt$w_s <- round(pt$w_s, 6)
  pt$treated_share <- round(pt$treated_share, 5); pt$post_share <- round(pt$post_share, 5)
  pt$beta_s <- round(pt$beta_s, 6)
  print(pt, row.names = FALSE)
  cat(sprintf("\n[step 3] stacks with V_s = 0 (no within-stack residual variation): %s\n",
              if (length(res$zero_stacks)) paste(res$zero_stacks, collapse = ", ") else "(none)"))
  cat(sprintf("\n[step 5] pooled beta (unit x stack + time x stack) = %s\n", fpr(res$pooled)))
  cat(sprintf("[step 6] sum_s w_s * beta_s              = %s\n", fpr(res$recon)))
  cat(sprintf("[step 6] difference (recon - pooled)     = %s\n\n", fpr(res$recon - res$pooled)))
}

# =============================================================================
# UNWEIGHTED (steps 1-6)
# =============================================================================
u <- decomp(est, wname = NULL)
report(u, "UNWEIGHTED  (steps 1-6)")

# aside: residualize over ALL rows (incl. pre-2000 missing-outcome) -> gap
r1_all <- felm(no.req ~ 1 | agency_stack + year.cohort, data = full)
full$vhat_all <- as.numeric(r1_all$residuals[, 1])
fa <- full[!is.na(full[[Y]]), ]
stacks <- sort(unique(est$cohort))
Vs_all <- sapply(stacks, function(s) sum(fa$vhat_all[fa$cohort == s]^2))
ws_all <- Vs_all / sum(Vs_all)
recon_all <- sum(ifelse(Vs_all == 0 | is.na(u$tab$beta_s), 0, ws_all * u$tab$beta_s))
cat(sprintf("[aside] if vhat is residualized over ALL rows (incl. missing-outcome) then\n"))
cat(sprintf("        summed over the observed-outcome rows: sum w_s' beta_s = %s  (diff from pooled = %s)\n\n",
            fpr(recon_all), fpr(recon_all - u$pooled)))

# =============================================================================
# step 7: WEIGHTED (entropy balancing weights) throughout
# =============================================================================
wd <- decomp(est, wname = "weights")
report(wd, "WEIGHTED by entropy-balancing weights  (step 7)")

# =============================================================================
# step 8: additive unit FE (agency.id) + time-by-stack FE
# =============================================================================
cat("========================================================\n step 8: additive unit FE + time-by-stack FE\n========================================================\n")
m8 <- felm(any.fatalities ~ no.req | agency.id + year.cohort, data = est)
pooled8 <- unname(coef(m8)["no.req"])
nmulti <- sum(tapply(est$cohort, est$agency.id, function(x) length(unique(x))) > 1)
cat(sprintf("[step 8] pooled beta (agency.id additive + year.cohort) = %s\n", fpr(pooled8)))
cat(sprintf("[step 8] sum_s w_s * beta_s (from the R2 decomposition)  = %s\n", fpr(u$recon)))
cat(sprintf("[step 8] difference (recon - pooled_additive)           = %s\n", fpr(u$recon - pooled8)))
cat(sprintf("[step 8] agencies appearing in >1 stack: %d of %d\n",
            nmulti, length(unique(est$agency.id))))

# =============================================================================
# persist the per-stack tables and an identity summary
# =============================================================================
write.csv(u$tab,  file.path(OUT_DIR, "fwl_decomp_unweighted.csv"), row.names = FALSE)
write.csv(wd$tab, file.path(OUT_DIR, "fwl_decomp_weighted.csv"),   row.names = FALSE)
summ <- data.frame(
  spec = c("unweighted (R2)", "weighted (R2, ebal)", "additive unit FE + time x stack",
           "unweighted, vhat over all rows"),
  pooled = c(u$pooled, wd$pooled, pooled8, u$pooled),
  recon  = c(u$recon,  wd$recon,  u$recon,  recon_all),
  diff   = c(u$recon - u$pooled, wd$recon - wd$pooled, u$recon - pooled8, recon_all - u$pooled),
  identity_holds = c(TRUE, TRUE, FALSE, FALSE))
write.csv(summ, file.path(OUT_DIR, "fwl_decomp_summary.csv"), row.names = FALSE)
cat("\nWrote output/fwl_decomp_{unweighted,weighted,summary}.csv\n")
