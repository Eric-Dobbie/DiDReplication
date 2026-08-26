# =============================================================================
# stage_dn.R  --  Stage D (diagnostics) + Stage N (reporting-number reconcile)
#
# Stage D:
#   1. Effective number of cohorts  1/sum(w^2)  for R2, R2h, R5.
#   2. Wald test H0: all cohort effects equal, for CS and stacked, on both pools.
#   3. Leave-one-cohort-out for R2, R2h, R5 (flag sign flips); verify 2009 weight.
#   4. Weight concentration: top-1 and top-3 cohort weight share per estimator.
#   5. Cohorts with no clean controls under each pool.
# Stage N:
#   Confirm 278,324 file rows / 171,906 outcome-observed / 171,444 covariate-N;
#   reconciliation table.  (Flagging of mislabeled counts done separately by grep.)
# =============================================================================

suppressMessages({ library(did); library(fixest); library(lfe); library(dplyr) })
set.seed(0); options(stringsAsFactors = FALSE)

.here <- tryCatch({ a<-commandArgs(FALSE); f<-sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset="output")
logf <- file.path(OUT_DIR, "run.log")
lg <- function(...) { cat(..., "\n"); cat(..., "\n", file=logf, append=TRUE) }
lg("\n\n================ Stage D+N run.log ================")
lg("timestamp:", format(Sys.time(), tz="UTC", usetz=TRUE), " seed 0")

Y<-"any.fatalities"; TIME<-"year"; UNIT<-"agency.id"; UNITN<-"agency.num"
COHORT<-"year.changed"; S_TREAT<-"no.req"; S_TIME<-"year.cohort"; S_STACK<-"cohort"

# ---- prep dta (as run_extraction.R) ----------------------------------------
dta <- read.csv(file.path(DATA_DIR,"dta.csv"))
dta[[UNITN]] <- as.numeric(factor(dta[[UNIT]]))
dta[[COHORT]] <- ave(dta[[COHORT]], dta[[UNIT]],
  FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta[[COHORT]] <- ifelse(is.na(dta[[COHORT]]) & dta$no.req==0, 0,
                 ifelse(is.na(dta[[COHORT]]) & dta$no.req==1, 1987, dta[[COHORT]]))
ag_hist <- tapply(dta$no.req, dta$agency.id, function(x){u<-unique(x[!is.na(x)])
  if(length(u)==1 && u==0) "always0" else if(length(u)==1 && u==1) "always1" else "mixed"})
ag_nochange <- tapply(dta$change.type, dta$agency.id, function(ct) all(ct=="No Change"))
nochange_ids<- names(ag_nochange)[ag_nochange]
always0_ids <- intersect(names(ag_hist)[ag_hist=="always0"], nochange_ids)
always1_ids <- intersect(names(ag_hist)[ag_hist=="always1"], nochange_ids)

stacked <- read.csv(file.path(DATA_DIR,"stacked_fatal.csv"))
stacked$agency_stack <- paste(stacked[[UNIT]], stacked[[S_STACK]], sep="__")
stacked_h <- stacked[stacked$treat==1 | (stacked$treat==0 & stacked[[UNIT]] %in% always0_ids), ]

D <- list(); addD <- function(diagnostic, estimator, metric, value, note="")
  D[[length(D)+1L]] <<- data.frame(diagnostic=diagnostic, estimator=estimator,
    metric=metric, value=as.character(value), note=note)

# ---- helper: stacked cohort atoms + weights + full clustered vcov ----------
stacked_fit <- function(d) {
  est <- d[!is.na(d[[Y]]), ]
  reg <- lfe::felm(as.formula(sprintf(
    "%s ~ %s:factor(%s) | %s + %s | 0 | %s", Y, S_TREAT, S_STACK,
    "agency_stack", S_TIME, UNIT)), data=d)
  co <- coef(reg); V <- reg$clustervcv
  keep <- !is.na(co)
  co <- co[keep]; V <- V[keep, keep, drop=FALSE]
  grp <- as.integer(gsub(sprintf("%s:factor(%s)", S_TREAT, S_STACK), "",
                         names(co), fixed=TRUE))
  wv <- fwl_stack_weights(est, S_TREAT, "agency_stack", S_TIME, S_STACK)
  w  <- wv[as.character(grp)]; w <- w/sum(w)
  list(group=grp, beta=as.numeric(co), V=V, weight=as.numeric(w))
}
# ---- helper: CS fit -> cohort effects (post-avg) + weights + vcov ----------
cs_fit <- function(d) {
  mp <- did::att_gt(yname=Y, tname=TIME, idname=UNITN, gname=COHORT, xformla=~1,
                    control_group="nevertreated", clustervars=UNITN, data=d)
  keep <- !is.na(mp$att)
  att <- mp$att[keep]; g <- mp$group[keep]; t <- mp$t[keep]
  # analytical unit-clustered covariance of att(g,t) from the influence function
  # (mp$V_analytical is not the k x k att covariance): V = IF'IF / n^2, where the
  # 71 rows of inffunc are the clustering units. Matches mp$se to ~5% (bootstrap).
  IF <- as.matrix(mp$inffunc)[, keep, drop=FALSE]
  V <- crossprod(IF) / (mp$n^2)
  post <- t >= g
  gs <- sort(unique(g[post]))
  # cohort-average aggregation matrix A: row per cohort, avg its post cells
  A <- matrix(0, nrow=length(gs), ncol=length(att))
  for (i in seq_along(gs)) { sel <- post & g==gs[i]; A[i, sel] <- 1/sum(sel) }
  theta <- as.numeric(A %*% att); Vth <- A %*% V %*% t(A)
  # simple-aggregation cohort weights: pg_g * (#post cells g), normalized
  dp <- mp$DIDparams; idata <- dp$data
  keepu <- !duplicated(idata[[dp$idname]]); ug <- idata[[dp$gname]][keepu]
  pg <- sapply(gs, function(gg) mean(ug==gg))
  ncell <- sapply(gs, function(gg) sum(post & g==gg))
  W <- pg*ncell; W <- W/sum(W)
  list(group=gs, theta=theta, Vth=Vth, weight=W, overall=sum(W*theta))
}

wald_equal <- function(theta, V) {
  k <- length(theta); if (k < 2) return(c(stat=NA, df=NA, p=NA))
  C <- cbind(-1, diag(k-1))                     # theta_i - theta_1
  d <- C %*% theta; M <- C %*% V %*% t(C)
  stat <- as.numeric(t(d) %*% solve(M) %*% d)
  p <- pchisq(stat, df=k-1, lower.tail=FALSE)
  c(stat=stat, df=k-1, p=p)
}
eff_cohorts <- function(w) 1/sum(w^2)
loo <- function(group, beta, weight) {
  sapply(seq_along(group), function(k){
    w <- weight[-k]; b <- beta[-k]; sum(w*b)/sum(w) })
}
conc <- function(w){ s<-sort(w, decreasing=TRUE); c(top1=s[1], top3=sum(head(s,3))) }

# =============================================================================
# Fit everything
# =============================================================================
lg("\n-- fitting stacked (741, 41) and CS (41, 741) --")
S2  <- stacked_fit(stacked)     # R2
S2h <- stacked_fit(stacked_h)   # R2h
C41 <- cs_fit(dta)              # R5 pool
dta_d <- dta; dta_d[[COHORT]] <- ifelse(dta_d$agency.id %in% always1_ids, 0, dta_d[[COHORT]])
C741<- cs_fit(dta_d)           # R5d pool

# --- D1 effective cohorts ---
for (nm in c("R2","R2h","R5")) {
  w <- switch(nm, R2=S2$weight, R2h=S2h$weight, R5=C41$weight)
  addD("D1_effective_cohorts", nm, "eff_n_cohorts", round(eff_cohorts(w),4),
       sprintf("k=%d", length(w)))
  lg(sprintf("D1 %-4s effective cohorts = %.3f (of %d)", nm, eff_cohorts(w), length(w)))
}

# --- D2 Wald equal cohort effects ---
lg("\n-- D2 Wald: H0 all cohort effects equal --")
W_S2  <- wald_equal(S2$beta,  S2$V)
W_S2h <- wald_equal(S2h$beta, S2h$V)
W_C41 <- wald_equal(C41$theta, C41$Vth)
W_C741<- wald_equal(C741$theta, C741$Vth)
for (x in list(list("stacked","741(R2)",W_S2), list("stacked","41(R2h)",W_S2h),
               list("CS","41(R5)",W_C41), list("CS","741(R5d)",W_C741))) {
  s<-x[[3]]; rej <- is.finite(s["p"]) && s["p"]<0.05
  addD("D2_wald_equal_cohorts", sprintf("%s_%s",x[[1]],x[[2]]), "chi2_stat",
       round(s["stat"],3), sprintf("df=%d p=%.4g dispersion_exceeds_sampling=%s",
       s["df"], s["p"], rej))
  lg(sprintf("D2 %-8s pool %-8s chi2=%.2f df=%d p=%.4g reject=%s",
     x[[1]], x[[2]], s["stat"], s["df"], s["p"], rej))
}

# --- D3 leave-one-cohort-out ---
lg("\n-- D3 leave-one-cohort-out (pooled after dropping each cohort) --")
loo_report <- function(nm, group, beta, weight, full) {
  vals <- loo(group, beta, weight)
  for (i in order(group)) {
    flip <- sign(vals[i]) != sign(full)
    addD("D3_leave_one_cohort_out", nm, sprintf("drop_%d", group[i]),
         round(vals[i],5), sprintf("full=%.5f sign_flip=%s weight=%.4f",
         full, flip, weight[i]))
    if (flip) lg(sprintf("  FLAG %s drop %d -> %.5f (flips sign; full %.5f)",
                 nm, group[i], vals[i], full))
  }
  vals
}
full_R2  <- sum(S2$weight*S2$beta)
full_R2h <- sum(S2h$weight*S2h$beta)
full_R5  <- C41$overall
loo_report("R2",  S2$group,  S2$beta,  S2$weight,  full_R2)
loo_report("R2h", S2h$group, S2h$beta, S2h$weight, full_R2h)
loo_report("R5",  C41$group, C41$theta, C41$weight, full_R5)
w2009_R2  <- S2$weight[S2$group==2009]
w2009_R2h <- S2h$weight[S2h$group==2009]
lg(sprintf("D3 cohort-2009 stacked weight: R2=%.4f  R2h=%.4f (assumptions.md §9: ~0.37)",
   w2009_R2, w2009_R2h))
addD("D3_cohort2009_weight","R2","weight_2009", round(w2009_R2,4), "assumptions.md §9 ~0.37")
addD("D3_cohort2009_weight","R2h","weight_2009", round(w2009_R2h,4), "")

# --- D4 weight concentration ---
lg("\n-- D4 weight concentration --")
for (x in list(list("R2",S2$weight), list("R2h",S2h$weight), list("R5",C41$weight),
               list("R5d",C741$weight))) {
  cc <- conc(x[[2]])
  addD("D4_weight_concentration", x[[1]], "top1_share", round(cc["top1"],4), "")
  addD("D4_weight_concentration", x[[1]], "top3_share", round(cc["top3"],4), "")
  lg(sprintf("D4 %-4s top1=%.3f top3=%.3f", x[[1]], cc["top1"], cc["top3"]))
}

# --- D5 cohorts with no clean controls under each pool ---
lg("\n-- D5 cohorts with no clean controls --")
allcoh <- sort(unique(stacked[[S_STACK]]))
ctrl_by_cohort_741 <- sapply(allcoh, function(g) sum(stacked$treat==0 & stacked[[S_STACK]]==g))
ctrl_by_cohort_41  <- sapply(allcoh, function(g) sum(stacked_h$treat==0 & stacked_h[[S_STACK]]==g))
no_ctrl_741 <- allcoh[ctrl_by_cohort_741==0]
no_ctrl_41  <- allcoh[ctrl_by_cohort_41==0]
# also: identified cohorts (atoms) vs eligible
lost_741 <- setdiff(allcoh, S2$group)
lost_41  <- setdiff(allcoh, S2h$group)
lg("all stacked cohorts:", paste(allcoh, collapse=","))
lg("741-pool cohorts with 0 control rows:", paste(no_ctrl_741, collapse=","))
lg("41-pool  cohorts with 0 control rows:", paste(no_ctrl_41, collapse=","))
lg("741-pool cohorts not identified (no atom):", paste(lost_741, collapse=","))
lg("41-pool  cohorts not identified (no atom):", paste(lost_41, collapse=","))
addD("D5_no_clean_controls","stacked_741","cohorts_zero_control_rows", paste(no_ctrl_741,collapse=";"),"")
addD("D5_no_clean_controls","stacked_41","cohorts_zero_control_rows", paste(no_ctrl_41,collapse=";"),"")
addD("D5_no_clean_controls","stacked_741","cohorts_unidentified", paste(lost_741,collapse=";"),"treated-only or collinear")
addD("D5_no_clean_controls","stacked_41","cohorts_unidentified", paste(lost_41,collapse=";"),"")

diag <- do.call(rbind, D)
write.csv(diag, file.path(OUT_DIR,"diagnostics.csv"), row.names=FALSE)
lg("\nWrote diagnostics.csv (", nrow(diag), "rows )")

# =============================================================================
# Stage N -- reporting-number reconciliation
# =============================================================================
lg("\n================ Stage N ================")
file_rows <- nrow(stacked)
outcome_obs <- sum(!is.na(stacked[[Y]]))
covars <- c("log.pop","log.med.inc","pct.white","pct.white.officers.imputed")
complete_cov <- outcome_obs - sum(!is.na(stacked[[Y]]) &
                  !complete.cases(stacked[, covars]))
# felm N under plain and covariate specs (confirm)
m_plain <- lfe::felm(any.fatalities ~ no.req | agency.id + year.cohort | 0 | agency.id,
                     data=stacked)
m_cov <- lfe::felm(any.fatalities ~ no.req + log.pop + log.med.inc + pct.white +
                   pct.white.officers.imputed | agency.id + year.cohort | 0 | agency.id,
                   data=stacked)
N_plain <- m_plain$N; N_cov <- m_cov$N
control_only <- 741*21*11
treated_obs <- outcome_obs - control_only
lg(sprintf("file rows                : %d (expect 278324)", file_rows))
lg(sprintf("outcome-observed rows    : %d (expect 171906)", outcome_obs))
lg(sprintf("  = control 741*21*11    : %d  + treated %d", control_only, treated_obs))
lg(sprintf("felm N plain             : %d (expect 171906)", N_plain))
lg(sprintf("felm N covariate (m3/m4) : %d (expect 171444)", N_cov))
lg(sprintf("dropped for missing cov  : %d (expect 462)", outcome_obs - N_cov))

nrec <- data.frame(
  quantity=c("file_rows","outcome_observed_rows","control_rows(741x21x11)",
             "treated_outcome_rows","felm_N_plain","felm_N_covariate",
             "rows_dropped_missing_pctwhiteofficers"),
  value=c(file_rows, outcome_obs, control_only, treated_obs, N_plain, N_cov,
          outcome_obs-N_cov),
  expected=c(278324,171906,171171,735,171906,171444,462),
  kind=c("file_rows","file_rows(outcome-observed)","file_rows(control only)",
         "file_rows(treated)","estimation_N","estimation_N","dropped_rows"),
  stringsAsFactors=FALSE)
nrec$is_estimation_N <- nrec$kind=="estimation_N"
nrec$matches_expected <- nrec$value==nrec$expected

# append the windowed Cengiz row counts (stacked_dimensionality.csv) so the
# single table also carries the 40-48k values the slides put next to 171,444
dimf <- file.path(OUT_DIR,"stacked_dimensionality.csv")
if (file.exists(dimf)) {
  dm <- read.csv(dimf)
  dmrows <- data.frame(
    quantity=sprintf("cengiz_window_rows[%s,k=%d]", dm$control, dm$k),
    value=dm$rows, expected=dm$rows,
    kind="row_count_windowed", is_estimation_N=FALSE, matches_expected=TRUE,
    stringsAsFactors=FALSE)
  nrec <- rbind(nrec, dmrows)
}
write.csv(nrec, file.path(OUT_DIR,"n_reconciliation.csv"), row.names=FALSE)
lg("\nWrote n_reconciliation.csv (", nrow(nrec), "rows )")
lg(paste(capture.output(print(nrec[nrec$kind!="row_count_windowed",], row.names=FALSE)), collapse="\n"))

# labeling audit of output/: find files whose column/label "Rows"|"N"|"Observations"
# attaches to BOTH a windowed row count and an estimation N (the slide error).
lg("\n-- Stage N labeling audit of output/ --")
est_N_vals <- c("171444","171,444","171906","171,906")
audit <- data.frame()
for (f in list.files(OUT_DIR, pattern="\\.(csv|tex|md)$", full.names=TRUE)) {
  txt <- tryCatch(readLines(f, warn=FALSE), error=function(e) character(0))
  has_rowslabel <- any(grepl("\\b[Rr]ows\\b|Observations|Num\\.? ?[Oo]bs", txt))
  has_estN <- any(sapply(est_N_vals, function(v) any(grepl(v, txt, fixed=TRUE))))
  # collision = same file labels windowed row-counts "Rows" AND carries an estimation N
  has_windowed <- any(grepl("\\b4[0-8][0-9]{3}\\b|33,?470|36,?428|14,?21[0-9]", txt))
  if (has_rowslabel || has_estN) audit <- rbind(audit, data.frame(
    file=basename(f), labels_rows=has_rowslabel, has_estimation_N=has_estN,
    has_windowed_rowcount=has_windowed,
    collision=has_rowslabel && has_estN && has_windowed, stringsAsFactors=FALSE))
}
write.csv(audit, file.path(OUT_DIR,"n_label_audit.csv"), row.names=FALSE)
lg(paste(capture.output(print(audit, row.names=FALSE)), collapse="\n"))
lg("collisions found in output/:", sum(audit$collision))

lg("\n================ Stage D+N sessionInfo ================")
lg(paste(capture.output(print(sessionInfo())), collapse="\n"))
cat("\nDONE Stage D+N\n")
