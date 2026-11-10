# =============================================================================
# taskA1_wald_diagnostics.R
#
# Diagnose the per-cohort joint pre-trend Wald tests (CS and SA). For each cohort
# report: # pre-treatment coefficients, Wald statistic, df used, numerical rank
# of the vcov block, its condition number, and the p-value. Decision rule: if
# rank < n_pre the p-value is NOT interpretable -> report NA (no pseudo-inverse
# p-value). Also report the ginv tolerance in use and the singular-value spectra
# for the 2002 and 2009 blocks.
# =============================================================================
suppressMessages({library(did); library(fixest); library(MASS)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987, dta$year.changed))

GINV_TOL <- sqrt(.Machine$double.eps)      # MASS::ginv default tolerance
cat(sprintf("ginv default tolerance = %.3e (drops singular values <= tol * max(sv))\n\n", GINV_TOL))

# full-rank-aware Wald: solve() if full rank, else NA (per decision rule)
wald_block <- function(a, V) {
  a <- as.numeric(a); n <- length(a)
  out <- list(n_pre = n, wald = NA_real_, df = NA_integer_, rank = NA_integer_,
              cond = NA_real_, p = NA_real_, ginv_p = NA_real_, sv = NULL)
  if (n == 0) return(out)
  V <- as.matrix(V)
  sv <- svd(V)$d; out$sv <- sv
  rk <- sum(sv > GINV_TOL * max(sv))        # same rule ginv uses
  out$rank <- rk
  out$cond <- if (min(sv) > 0) sv[1]/sv[length(sv)] else Inf
  if (rk < n) {                              # rank-deficient -> not interpretable
    # also compute the ginv-based number that was being (wrongly) reported
    Vi <- MASS::ginv(V); stat <- as.numeric(t(a) %*% Vi %*% a)
    out$wald <- stat; out$df <- rk; out$p <- NA_real_   # p = NA by rule
    out$ginv_p <- pchisq(stat, df = rk, lower.tail = FALSE)
    return(out)
  }
  Vi <- solve(V); stat <- as.numeric(t(a) %*% Vi %*% a)
  out$wald <- stat; out$df <- n; out$p <- pchisq(stat, df = n, lower.tail = FALSE)
  out$ginv_p <- out$p
  out
}

# ---- CS blocks (CORRECTED scaling: V_analytical / n; agency=unit => clustered) --
run_cs <- function(bp, tag) {
  set.seed(0)
  mp <- did::att_gt(yname="any.fatalities", tname="year", idname="agency.num",
                    gname="year.changed", xformla=~1, control_group="nevertreated",
                    clustervars="agency.num", base_period=bp, data=dta)
  n <- mp$DIDparams$n
  V <- as.matrix(mp$V_analytical) / n       # <-- correct scale: se = sqrt(diag(V)/n)
  att <- mp$att; G <- mp$group; Tt <- mp$t
  groups <- sort(unique(G[G>0]))
  ntr <- sapply(groups, function(g) length(unique(dta$agency.id[dta$year.changed==g])))
  rr <- list(); sp <- list()
  for (i in seq_along(groups)) { g <- groups[i]
    idx <- which(G==g & Tt < g & !is.na(att))
    w <- wald_block(att[idx], V[idx, idx, drop=FALSE])
    if (g %in% c(2002,2009)) sp[[paste0(tag,"_",g)]] <- w$sv
    rr[[length(rr)+1]] <- data.frame(estimator=tag, cohort=g, n_treated=ntr[i],
      n_pre=w$n_pre, rank=w$rank, cond=w$cond, wald=w$wald, df=w$df, p=w$p, ginv_p=w$ginv_p)
  }
  list(rows=do.call(rbind, rr), spectra=sp, mp=mp)
}
csV <- run_cs("varying",   "CS(varying)")
csU <- run_cs("universal", "CS(universal)")
cat("CS base_period as originally run in Task 4.2:", csV$mp$DIDparams$base_period, "\n\n")
rows <- c(list(csV$rows, csU$rows)); spectra <- c(csV$spectra, csU$spectra)
groups <- sort(unique(csV$mp$group[csV$mp$group>0]))
ntr <- sapply(groups, function(g) length(unique(dta$agency.id[dta$year.changed==g])))

# ---- SA blocks --------------------------------------------------------------
d_sa <- dta[!is.na(dta$any.fatalities) & (dta$year.changed==0 | dta$year.changed>=2001),]
d_sa$coh <- ifelse(d_sa$year.changed==0, 10000, d_sa$year.changed)
res <- fixest::feols(any.fatalities ~ sunab(coh, year) | agency.num + year,
                     data=d_sa, cluster=~agency.num)
ct <- summary(res, agg=FALSE)$coeftable; Vsa <- as.matrix(vcov(res))[rownames(ct), rownames(ct)]
mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee <- as.integer(vapply(mm, function(z) z[2], character(1)))
gg <- as.integer(vapply(mm, function(z) z[3], character(1)))
keep <- !is.na(ee) & !is.na(gg)
for (g in sort(unique(gg[keep]))) {
  idx <- which(gg==g & ee<0 & keep)
  w <- wald_block(ct[idx,1], Vsa[idx, idx, drop=FALSE])
  if (g %in% c(2002,2009)) spectra[[paste0("SA_",g)]] <- w$sv
  rows[[length(rows)+1]] <- data.frame(estimator="SA", cohort=g,
    n_treated=ntr[match(g,groups)], n_pre=w$n_pre, rank=w$rank, cond=w$cond,
    wald=w$wald, df=w$df, p=w$p, ginv_p=w$ginv_p)
}
tab <- do.call(rbind, rows); row.names(tab) <- NULL
tab$rank_deficient <- ifelse(!is.na(tab$n_pre) & tab$n_pre>0 & tab$rank < tab$n_pre, "YES", "")

options(width = 200)
cat("== Per-cohort pre-trend Wald diagnostics ==\n")
print(transform(tab, cond=signif(cond,3), wald=signif(wald,4),
                p=round(p,4), ginv_p=round(ginv_p,4)), row.names=FALSE)
cat(sprintf("\nrank-deficient blocks: CS %d/%d, SA %d/%d\n",
    sum(tab$rank_deficient=="YES" & tab$estimator=="CS"), sum(tab$estimator=="CS" & tab$n_pre>0),
    sum(tab$rank_deficient=="YES" & tab$estimator=="SA"), sum(tab$estimator=="SA" & tab$n_pre>0)))
cat("\nsingular-value spectra (2002, 2009 blocks):\n")
for (nm in names(spectra)) cat(sprintf("  %-9s: %s\n", nm, paste(signif(spectra[[nm]],3), collapse=", ")))

write.csv(tab, file.path(OUT_DIR, "taskA1_wald_diagnostics.csv"), row.names=FALSE)
cat("\nWrote output/taskA1_wald_diagnostics.csv\n")
