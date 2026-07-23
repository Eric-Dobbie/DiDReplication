# =============================================================================
# extract_atoms.R
#
# Data-agnostic extraction of the disaggregated "atoms" (component estimates)
# and their aggregation weights for three staggered-adoption DiD estimators:
#
#   1. Callaway & Sant'Anna (CS)   -- ATT(g,t)
#   2. Sun & Abraham        (SA)   -- CATT(g,e)
#   3. Stacked regression   (stacked) -- cohort-specific stack effects
#
# Every function takes column-name arguments explicitly; NO variable names from
# any particular dataset are hardcoded here. Dataset-specific preparation
# (recoding never-treated, choosing the estimation sample, etc.) lives in the
# driver script (run_extraction.R), not in these functions.
#
# Each extractor returns a tidy data frame with the common schema
#   estimator, atom_id, group, time, estimate, se, weight
# and carries the estimator's own package-reported pooled/overall estimate as
# attributes:  attr(df, "pooled_estimate"), attr(df, "pooled_se"),
#              attr(df, "pooled_label").
# The recombination check (Step 5) compares sum(weight * estimate) to that
# attribute.
# =============================================================================

suppressMessages({
  library(did)
  library(fixest)
  library(lfe)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Step 2 -- Callaway & Sant'Anna atoms
# -----------------------------------------------------------------------------
# yname, tname, idname, gname : column names (gname = first-treatment period,
#   with never-treated coded 0). xformla : covariate formula (default ~1).
# control_group, clustervars, ... : passed through to did::att_gt.
#
# Weight approach:  the "simple" aggregation in did::aggte weights each ATT(g,t)
# with t >= g by  pg_g / sum_{keepers} pg  , where pg_g = P(G = g) is the share
# of units in group g (see did:::compute.aggte, type == "simple":
#   simple.att <- sum(att[keepers] * pg[keepers]) / sum(pg[keepers]) ).
# We EXTRACT pg directly from the fitted object's internal one-row-per-unit
# analysis data (dp$data) rather than reconstructing group sizes from the raw
# input, because att_gt drops already-treated units and rows with missing
# outcomes -- so raw group counts do not match the counts did actually uses.
# This makes the recombination exact (to ~1e-16) instead of merely approximate.
# -----------------------------------------------------------------------------
extract_cs_atoms <- function(data, yname, tname, idname, gname,
                             xformla = ~1, control_group = "nevertreated",
                             clustervars = NULL, seed = 0, ...) {

  # att_gt uses a multiplier bootstrap for inference by default; seed it so the
  # reported SEs are reproducible (point estimates/weights are deterministic).
  set.seed(seed)
  mp <- did::att_gt(yname = yname, tname = tname, idname = idname, gname = gname,
                    xformla = xformla, control_group = control_group,
                    clustervars = clustervars, data = data, ...)

  agg <- did::aggte(mp, type = "simple", na.rm = TRUE)

  # --- exact simple-aggregation weights: pull pg from the fitted object -------
  dp     <- mp$DIDparams
  idata  <- dp$data
  gcol   <- dp$gname
  icol   <- dp$idname
  wcol   <- if (".w" %in% names(idata)) idata$.w else rep(1, nrow(idata))
  keepu  <- !duplicated(idata[[icol]])          # one row per unit
  ug     <- idata[[gcol]][keepu]                 # unit-level group
  uw     <- wcol[keepu]                          # unit-level weight (1 if none)

  atoms <- data.frame(group = mp$group, time = mp$t,
                      estimate = mp$att, se = mp$se,
                      stringsAsFactors = FALSE)
  atoms <- atoms[!is.na(atoms$estimate), , drop = FALSE]   # na.rm, as in did

  # pg_g = weighted share of units in group g (denominator cancels on normalize)
  pg <- sapply(atoms$group, function(g) sum(uw * (ug == g)) / sum(uw))
  keepers <- atoms$time >= atoms$group                     # post-treatment cells
  w <- rep(0, nrow(atoms))
  w[keepers] <- pg[keepers] / sum(pg[keepers])
  atoms$weight <- w

  out <- data.frame(
    estimator = "CS",
    atom_id   = sprintf("g%d_t%d", atoms$group, atoms$time),
    group     = atoms$group,
    time      = atoms$time,
    estimate  = atoms$estimate,
    se        = atoms$se,
    weight    = atoms$weight,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(agg$overall.att)
  attr(out, "pooled_se")       <- as.numeric(agg$overall.se)
  attr(out, "pooled_label")    <- "did::aggte(type='simple') overall ATT"
  out
}

# -----------------------------------------------------------------------------
# Step 3 -- Sun & Abraham atoms
# -----------------------------------------------------------------------------
# gname = cohort/first-treatment period; never-treated units must be flagged so
# they act as pure controls. Pass `never_value` = the code used for never-treated
# in gname (default 0); it is internally recoded to `ref_c` (a period far beyond
# the sample) so fixest::sunab treats those units as never treated.
#
# The atoms are the disaggregated cohort x relative-time coefficients CATT(g,e),
# read from summary(res, agg = FALSE). The overall ATT is summary(res,
# agg = "att").
#
# Weight approach:  fixest's agg = "att" (see fixest:::aggregate.fixest) sets
#   shares = colSums(sign(model_matrix[, e>=0 cols]))  (weighted by obs weights
#   if present), then shares/sum(shares). I.e. each CATT(g,e) with e >= 0 is
# weighted by the number of estimation-sample observations loading on that
# cohort x relative-time cell. We reconstruct exactly that from model.matrix().
# -----------------------------------------------------------------------------
extract_sa_atoms <- function(data, yname, idname, tname, gname,
                             never_value = 0, ref_c = 10000,
                             cluster = NULL, covars = NULL, ...) {

  d <- data
  d[["._coh"]] <- ifelse(d[[gname]] == never_value, ref_c, d[[gname]])

  # optional covariates enter as linear controls alongside the sunab term
  rhs_cov <- if (is.null(covars)) "" else paste0(paste(covars, collapse = " + "), " + ")
  fml <- stats::as.formula(sprintf(
    "%s ~ %ssunab(._coh, %s) | %s + %s", yname, rhs_cov, tname, idname, tname))
  cl  <- if (is.null(cluster)) stats::as.formula(paste0("~", idname)) else cluster
  res <- fixest::feols(fml, data = d, cluster = cl, ...)

  att <- fixest::coeftable(summary(res, agg = "att"))
  att_val <- att["ATT", 1]; att_se <- att["ATT", 2]

  ct  <- summary(res, agg = FALSE)$coeftable          # disaggregated CATT(g,e)
  cn  <- rownames(ct)
  est <- ct[, 1]; se <- ct[, 2]

  # robust parse of names like "year::-3:cohort::2009" (sub() backrefs are
  # unreliable with '::' separators here, so use regexec/regmatches)
  mm_parse <- regmatches(cn, regexec("::(-?[0-9]+):.*::([0-9]+)$", cn))
  e <- as.integer(vapply(mm_parse, function(z) z[2], character(1)))
  g <- as.integer(vapply(mm_parse, function(z) z[3], character(1)))

  # keep only the sunab cohort x period coefficients (covariate rows, if any,
  # do not parse to an (e, g) and are dropped)
  is_atom <- !is.na(e) & !is.na(g)
  e <- e[is_atom]; g <- g[is_atom]; est <- est[is_atom]; se <- se[is_atom]
  cn <- cn[is_atom]

  # exact fixest att weights: observation counts per cohort x period cell
  mmx    <- stats::model.matrix(res)
  shares <- colSums(abs(sign(mmx)))[cn]
  sel    <- (e >= 0) & !is.na(est)
  w      <- rep(0, length(est))
  w[sel] <- shares[sel] / sum(shares[sel])

  out <- data.frame(
    estimator = "SA",
    atom_id   = sprintf("g%d_e%d", g, e),
    group     = g,
    time      = e,                     # relative event time
    estimate  = as.numeric(est),
    se        = as.numeric(se),
    weight    = w,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(att_val)
  attr(out, "pooled_se")       <- as.numeric(att_se)
  attr(out, "pooled_label")    <- "fixest::sunab(agg='att') overall ATT"
  out
}

# -----------------------------------------------------------------------------
# FWL weight helper (used by the stacked extractor)
# -----------------------------------------------------------------------------
# Standard OLS/TWFE weighting result: the pooled coefficient on a treatment
# indicator equals a weighted average of stack-specific coefficients, with
# weights proportional to each stack's share of the total variance of the
# treatment indicator AFTER residualizing on the model's fixed effects
# (Frisch-Waugh-Lovell). We implement it explicitly:
#
#   1. residualize `treatname` on the fixed effects (idname + stacktime), via
#      OLS with no other covariates -- resid_treat.
#   2. weight_g = sum_g(resid_treat^2) / sum(resid_treat^2).
#
# IMPORTANT: the residualization MUST use the same rows the outcome regression
# uses. If the pooled/cohort regression drops rows with a missing outcome, those
# rows must be dropped here too, otherwise the weights are computed on a
# different sample than the coefficients and the FWL identity fails. `data`
# should therefore be the ESTIMATION sample (rows with non-missing outcome).
# -----------------------------------------------------------------------------
fwl_stack_weights <- function(data, treatname, idname, stacktime, stackid) {
  fml <- stats::as.formula(sprintf("%s ~ 1 | %s + %s", treatname, idname, stacktime))
  rt  <- lfe::felm(fml, data = data)$residuals[, 1]
  df  <- data.frame(stackid = data[[stackid]], rt2 = rt^2)
  w   <- df %>%
    dplyr::group_by(stackid) %>%
    dplyr::summarize(w_num = sum(rt2), .groups = "drop") %>%
    dplyr::mutate(weight = w_num / sum(w_num))
  stats::setNames(w$weight, as.character(w$stackid))
}

# -----------------------------------------------------------------------------
# Step 4 -- Stacked regression atoms
# -----------------------------------------------------------------------------
# Operates on an ALREADY-STACKED dataset -- i.e. one "sub-experiment" per
# treatment cohort, each combining that cohort's treated units with a clean
# control pool, stacked into one long frame with a stack id. (In the Payson &
# Parinandi data this is `stacked_fatal.csv`, built upstream with never-treated
# "clean controls" and a +/-4 balancing window; see run_extraction.R and the
# accompanying note. `build_stacks()` below reconstructs such a frame from a raw
# panel when one is not provided.)
#
#   stacktime : stack-specific time FE (e.g. year x cohort)   -> tname role
#   stackid   : sub-experiment id (e.g. treatment cohort)     -> gname role
#   treatname : 0/1 treatment indicator used in the pooled spec
#
# Atoms: stack-specific coefficients from the fully-interacted regression
#   y ~ treat:factor(stackid) | idname + stacktime         (cluster by idname).
# Weights: FWL variance shares from fwl_stack_weights(), computed on the
# estimation sample (resid_on = "estimation", the default) or, to reproduce a
# naive implementation, on all rows (resid_on = "all").
#
# Stacks whose interacted coefficient is NA/NaN (no identifying variation, e.g.
# a stack with no clean controls and no within-stack treatment contrast) are
# dropped, and the surviving weights are renormalized to sum to 1.
# -----------------------------------------------------------------------------
extract_stacked_atoms <- function(data, yname, idname, stacktime, stackid,
                                  treatname, cluster = idname,
                                  resid_on = c("estimation", "all"),
                                  window = NULL, control_group = "never_treated",
                                  ...) {
  resid_on <- match.arg(resid_on)

  # estimation sample = rows with a non-missing outcome (what felm actually uses)
  est <- data[!is.na(data[[yname]]), , drop = FALSE]

  # --- stack-specific (fully interacted) coefficients -----------------------
  fml <- stats::as.formula(sprintf(
    "%s ~ %s:factor(%s) | %s + %s | 0 | %s",
    yname, treatname, stackid, idname, stacktime, cluster))
  reg <- lfe::felm(fml, data = data)
  co  <- coef(reg); se_all <- reg$se
  nm  <- gsub(sprintf("%s:factor(%s)", treatname, stackid), "",
              names(co), fixed = TRUE)
  atoms <- data.frame(group = as.integer(nm),
                      estimate = as.numeric(co),
                      se = as.numeric(se_all[names(co)]),
                      stringsAsFactors = FALSE)

  # --- FWL variance-share weights -------------------------------------------
  wdat <- if (resid_on == "estimation") est else data
  wvec <- fwl_stack_weights(wdat, treatname, idname, stacktime, stackid)
  atoms$weight <- wvec[as.character(atoms$group)]

  # drop unidentified stacks, renormalize weights over survivors
  atoms <- atoms[!is.na(atoms$estimate) & !is.nan(atoms$estimate), , drop = FALSE]
  atoms$weight <- atoms$weight / sum(atoms$weight)

  # --- pooled target: the plain (unweighted, no-covariate) stacked coef ------
  # This is the spec the FWL decomposition targets (matches the residualization).
  pf <- stats::as.formula(sprintf("%s ~ %s | %s + %s | 0 | %s",
                                  yname, treatname, idname, stacktime, cluster))
  preg <- lfe::felm(pf, data = data)
  pooled <- coef(preg)[treatname]
  pooled_se <- preg$se[treatname]

  out <- data.frame(
    estimator = "stacked",
    atom_id   = sprintf("stack%d", atoms$group),
    group     = atoms$group,
    time      = NA_integer_,           # stack effect is not time-specific
    estimate  = atoms$estimate,
    se        = atoms$se,
    weight    = atoms$weight,
    stringsAsFactors = FALSE
  )
  attr(out, "pooled_estimate") <- as.numeric(pooled)
  attr(out, "pooled_se")       <- as.numeric(pooled_se)
  attr(out, "pooled_label")    <- sprintf(
    "pooled felm(%s ~ %s | %s + %s) [unweighted, no covariates]",
    yname, treatname, idname, stacktime)
  out
}

# -----------------------------------------------------------------------------
# build_stacks() -- reconstruct a stacked frame from a raw long panel
# -----------------------------------------------------------------------------
# Provided for the data-agnostic "build the stack" requirement. NOT used by the
# P&P driver (which reads their pre-built stacked_fatal.csv), but documents /
# implements the stack-construction logic:
#   * for each treatment cohort g, take units first-treated in g PLUS a clean
#     control pool
#       - control_group = "never_treated": units never treated (gname == never_value)
#       - control_group = "not_yet_treated": units with gname > g (+ never-treated)
#   * optionally restrict each sub-experiment to a symmetric +/-`window` around g
#   * define stack-specific time id `stacktime` = paste(time, g), a within-stack
#     treatment indicator, and stack id `stackid` = g.
# Returns the stacked data frame; feed it to extract_stacked_atoms().
# -----------------------------------------------------------------------------
build_stacks <- function(data, idname, tname, gname, yname,
                         never_value = 0, control_group = "never_treated",
                         window = NULL) {
  control_group <- match.arg(control_group, c("never_treated", "not_yet_treated"))
  cohorts <- sort(unique(data[[gname]][data[[gname]] != never_value]))
  never   <- data[data[[gname]] == never_value, , drop = FALSE]

  pieces <- lapply(cohorts, function(g) {
    treated <- data[data[[gname]] == g, , drop = FALSE]
    if (control_group == "never_treated") {
      ctrl <- never
    } else {
      nyt  <- data[data[[gname]] > g & data[[gname]] != never_value, , drop = FALSE]
      ctrl <- rbind(nyt, never)
    }
    sub <- rbind(treated, ctrl)
    if (!is.null(window)) {
      sub <- sub[abs(sub[[tname]] - g) <= window, , drop = FALSE]
    }
    sub[["._stackid"]]   <- g
    sub[["._stacktime"]] <- paste0(sub[[tname]], "_", g)
    sub[["._treat"]]     <- as.integer(sub[[gname]] == g & sub[[tname]] >= g)
    sub
  })
  do.call(rbind, pieces)
}
