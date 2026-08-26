# =============================================================================
# stage_l_ladder.R  --  Stage L: harmonized ladder (new rows only)
#
# Existing rows R0-R5 are tabulated from output files (not recomputed).
# New rows computed here, reusing R/extract_atoms.R:
#   R2h  Stacked, interacted FE, controls restricted to the 41 always-no.req==0
#   R3h  Cengiz +/-4, controls = 41 always-0
#   R5d  CS with the FULL 741 pool (always-1 recoded gname=0)  [DIAGNOSTIC ONLY]
#   R5b  CS simple aggregation weights applied to stacked cohort atoms, x2
#        (once on R2's 741-pool atoms, once on R2h's 41-pool atoms)
#   R6   CS ATT(g) and ATT(e) on the 41-unit pool
#
# Writes output/ladder.csv and output/atoms_harmonized.csv; appends to run.log.
# =============================================================================

suppressMessages({ library(did); library(fixest); library(lfe); library(dplyr) })
set.seed(0)
options(stringsAsFactors = FALSE)

.here <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=","",a[grep("^--file=",a)])
  if (length(f)) dirname(normalizePath(f)) else "R" }, error=function(e) "R")
source(file.path(.here, "extract_atoms.R"))

DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUT_DIR, "run.log")
lg <- function(...) { cat(..., "\n"); cat(..., "\n", file = logf, append = TRUE) }
lg("\n\n================ Stage L run.log ================")
lg("timestamp:", format(Sys.time(), tz="UTC", usetz=TRUE), " seed: 0")
lg("R:", R.version.string, " did:", as.character(packageVersion("did")),
   " fixest:", as.character(packageVersion("fixest")), " lfe:", as.character(packageVersion("lfe")))

Zc <- qnorm(0.975)
supports <- function(est, se) is.finite(est) && is.finite(se) && est < 0 && abs(est/se) > Zc

# ---- column names (as in run_extraction.R) ---------------------------------
Y<-"any.fatalities"; TIME<-"year"; UNIT<-"agency.id"; UNITN<-"agency.num"
COHORT<-"year.changed"; S_TREAT<-"no.req"; S_TIME<-"year.cohort"; S_STACK<-"cohort"

# =============================================================================
# Prep dta exactly as run_extraction.R (first-treatment collapse + gname recode)
# =============================================================================
dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta[[UNITN]] <- as.numeric(factor(dta[[UNIT]]))
dta[[COHORT]] <- ave(dta[[COHORT]], dta[[UNIT]],
  FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta[[COHORT]] <- ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]]==0, 0,
                 ifelse(is.na(dta[[COHORT]]) & dta[["no.req"]]==1, 1987, dta[[COHORT]]))

# treatment-history sets among the 741 non-changers
ag_hist <- tapply(dta$no.req, dta$agency.id, function(x){u<-unique(x[!is.na(x)])
  if(length(u)==1 && u==0) "always0" else if(length(u)==1 && u==1) "always1" else "mixed"})
ag_nochange <- tapply(dta$change.type, dta$agency.id, function(ct) all(ct=="No Change"))
nochange_ids <- names(ag_nochange)[ag_nochange]
always0_ids <- intersect(names(ag_hist)[ag_hist=="always0"], nochange_ids)  # 41
always1_ids <- intersect(names(ag_hist)[ag_hist=="always1"], nochange_ids)  # 700
lg("always0 (41-pool):", length(always0_ids), " always1:", length(always1_ids))
stopifnot(length(always0_ids)==41L, length(always1_ids)==700L)

# =============================================================================
# Load stacked
# =============================================================================
stacked <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))
STACKUNIT <- "agency_stack"
stacked[[STACKUNIT]] <- paste(stacked[[UNIT]], stacked[[S_STACK]], sep="__")

stacked_est_summary <- function(d, label) {
  e <- d[!is.na(d[[Y]]), ]
  data.frame(label=label,
    n_obs=nrow(e),
    n_units=length(unique(e[[UNIT]])),
    n_control_units=length(unique(e[[UNIT]][e$treat==0])),
    stringsAsFactors=FALSE)
}

# ----------------------------------------------------------------------------
# R2 baseline (full 741 pool, interacted FE) -- computed live to form R2h - R2
# ----------------------------------------------------------------------------
lg("\n-- R2 (baseline, 741 pool, interacted FE) --")
st_R2 <- extract_stacked_atoms(stacked, yname=Y, idname=STACKUNIT,
          stacktime=S_TIME, stackid=S_STACK, treatname=S_TREAT, cluster=UNIT,
          resid_on="estimation")
R2_est <- attr(st_R2,"pooled_estimate"); R2_se <- attr(st_R2,"pooled_se")
R2_sum <- stacked_est_summary(stacked, "R2")
R2_cohorts <- sort(st_R2$group)
lg(sprintf("R2 pooled = %+.6f (se %.6f)  cohorts[%d]: %s",
   R2_est, R2_se, length(R2_cohorts), paste(R2_cohorts, collapse=",")))

# ----------------------------------------------------------------------------
# R2h: restrict controls to the 41 always-0 (treated rows unchanged)
# ----------------------------------------------------------------------------
lg("\n-- R2h (41 pool, interacted FE) --")
stacked_h <- stacked[stacked$treat==1 |
                     (stacked$treat==0 & stacked[[UNIT]] %in% always0_ids), ]
st_R2h <- extract_stacked_atoms(stacked_h, yname=Y, idname=STACKUNIT,
           stacktime=S_TIME, stackid=S_STACK, treatname=S_TREAT, cluster=UNIT,
           resid_on="estimation")
R2h_est <- attr(st_R2h,"pooled_estimate"); R2h_se <- attr(st_R2h,"pooled_se")
R2h_sum <- stacked_est_summary(stacked_h, "R2h")
R2h_cohorts <- sort(st_R2h$group)
cohorts_lost <- setdiff(R2_cohorts, R2h_cohorts)
lg(sprintf("R2h pooled = %+.6f (se %.6f)  cohorts retained[%d]: %s",
   R2h_est, R2h_se, length(R2h_cohorts), paste(R2h_cohorts, collapse=",")))
lg("R2h cohorts lost vs R2:", if(length(cohorts_lost)) paste(cohorts_lost,collapse=",") else "(none)")
lg(sprintf("R2h - R2 = %+.6f  (sign flips: %s)", R2h_est - R2_est,
   sign(R2h_est) != sign(R2_est)))

# ----------------------------------------------------------------------------
# R3h: Cengiz +/-4 with 41-pool controls (rebuild windowed stack)
# ----------------------------------------------------------------------------
lg("\n-- R3h (Cengiz +/-4, 41 pool) --")
K<-4; FIRST_OUT<-2000; LAST_OUT<-2020
REVERSIBLE <- c("memphis tennessee","portsmouth new hampshire")
build_cengiz <- function(control_ids) {
  d <- dta[!dta$agency.id %in% REVERSIBLE, ]
  allcoh <- sort(unique(d$year.changed[!is.na(d$year.changed) & d$year.changed>=2000]))
  coh <- allcoh[allcoh-K>=FIRST_OUT & allcoh+K<=LAST_OUT]
  nochange <- d[d$agency.id %in% control_ids, ]
  parts <- list()
  for (g in coh) {
    tr <- d[!is.na(d$year.changed) & d$year.changed==g & d$year>=g-K & d$year<=g+K, ]
    if (nrow(tr)==0) next
    tr$cohort<-g; tr$treat<-1L
    ct <- nochange[nochange$year>=g-K & nochange$year<=g+K, ]
    ct$cohort<-g; ct$treat<-0L
    parts[[length(parts)+1L]] <- rbind(tr, ct)
  }
  cz <- do.call(rbind, parts)
  cz$stackunit <- paste0(cz$agency.id,"__",cz$cohort)
  cz$stacktime <- paste0(cz$year,"__",cz$cohort)
  cz
}
cz_h <- build_cengiz(always0_ids)
st_R3h <- extract_stacked_atoms(cz_h, yname=Y, idname="stackunit",
           stacktime="stacktime", stackid="cohort", treatname=S_TREAT,
           cluster="agency.id", resid_on="estimation")
R3h_est <- attr(st_R3h,"pooled_estimate"); R3h_se <- attr(st_R3h,"pooled_se")
R3h_cohorts <- sort(st_R3h$group)
cz_h_est <- cz_h[!is.na(cz_h[[Y]]), ]
lg(sprintf("R3h pooled = %+.6f (se %.6f)  cohorts survived[%d]: %s",
   R3h_est, R3h_se, length(R3h_cohorts), paste(R3h_cohorts, collapse=",")))
lg("R3h eligible cohorts (window+outcome):",
   paste(sort(unique(cz_h$cohort)), collapse=","))

# ----------------------------------------------------------------------------
# R5d: CS with FULL 741 pool  [DIAGNOSTIC ONLY -- asserts a false treat history]
# ----------------------------------------------------------------------------
lg("\n-- R5d (CS, 741 pool, DIAGNOSTIC ONLY) --")
dta_d <- dta
dta_d[[COHORT]] <- ifelse(dta_d$agency.id %in% always1_ids, 0, dta_d[[COHORT]])
cs_d <- extract_cs_atoms(dta_d, yname=Y, tname=TIME, idname=UNITN, gname=COHORT,
          xformla=~1, control_group="nevertreated", clustervars=UNITN)
R5d_est <- attr(cs_d,"pooled_estimate"); R5d_se <- attr(cs_d,"pooled_se")
d5d_est <- dta_d[!is.na(dta_d[[Y]]) & (dta_d[[COHORT]]==0 | dta_d[[COHORT]]>=2001), ]
R5d_nobs <- sum(!is.na(dta_d[[Y]]))  # att_gt effective outcome-observed rows
R5d_ncontrol <- length(unique(dta_d$agency.id[dta_d[[COHORT]]==0]))
R5d_ncohort <- length(unique(cs_d$group))
lg(sprintf("R5d overall ATT = %+.6f (se %.6f)  control units(gname=0): %d  cohorts: %d",
   R5d_est, R5d_se, R5d_ncontrol, R5d_ncohort))

# ----------------------------------------------------------------------------
# R5 baseline (41 pool CS) computed live for R5d diff, R5b weights, R6
# ----------------------------------------------------------------------------
lg("\n-- R5 baseline (CS 41 pool) live --")
cs <- extract_cs_atoms(dta, yname=Y, tname=TIME, idname=UNITN, gname=COHORT,
        xformla=~1, control_group="nevertreated", clustervars=UNITN)
R5_est <- attr(cs,"pooled_estimate"); R5_se <- attr(cs,"pooled_se")
R5_ncontrol <- length(unique(dta$agency.id[dta[[COHORT]]==0]))
R5_ncohort <- length(unique(cs$group))
lg(sprintf("R5 overall ATT = %+.6f (se %.6f)  control units: %d  cohorts: %d",
   R5_est, R5_se, R5_ncontrol, R5_ncohort))
lg(sprintf("R5d - R5 = %+.6f", R5d_est - R5_est))

# CS simple cohort weights (sum atom weights within group)
cs_cohort_w <- tapply(cs$weight, cs$group, sum)
cs_cohort_w <- cs_cohort_w[cs_cohort_w > 0]

# ----------------------------------------------------------------------------
# R5b: apply CS simple cohort weights to stacked cohort atoms (aggregate_atoms)
# ----------------------------------------------------------------------------
aggregate_atoms <- function(atoms, weights_by_group) {
  g <- as.character(atoms$group)
  w <- weights_by_group[g]
  keep <- !is.na(w) & is.finite(atoms$estimate)
  w <- w[keep]/sum(w[keep])
  list(value=sum(w*atoms$estimate[keep]),
       matched=atoms$group[keep], w=w)
}
lg("\n-- R5b (CS-simple weights on stacked atoms) --")
r5b_R2  <- aggregate_atoms(st_R2,  cs_cohort_w)
r5b_R2h <- aggregate_atoms(st_R2h, cs_cohort_w)
lg(sprintf("R5b on R2  atoms (741 pool): %+.6f  [cohorts matched: %s]",
   r5b_R2$value,  paste(sort(r5b_R2$matched),  collapse=",")))
lg(sprintf("R5b on R2h atoms (41 pool):  %+.6f  [cohorts matched: %s]",
   r5b_R2h$value, paste(sort(r5b_R2h$matched), collapse=",")))
lg(sprintf("R5 (CS pooled target) = %+.6f", R5_est))
move_R2  <- abs(r5b_R2$value  - R5_est) < abs(R2_est  - R5_est)
move_R2h <- abs(r5b_R2h$value - R5_est) < abs(R2h_est - R5_est)
lg(sprintf("R5b-R2h moves toward R5: %s | R5b-R2 moves toward R5: %s",
   move_R2h, move_R2))
lg("FLAG: pool-driven" , if(move_R2h && !move_R2) "TRIGGERED" else "not (this pattern)")
lg("FLAG: aggregation-driven", if(move_R2h && move_R2) "TRIGGERED" else "not (this pattern)")

# ----------------------------------------------------------------------------
# R6: CS ATT(g) and ATT(e) on the 41-unit pool
# ----------------------------------------------------------------------------
lg("\n-- R6 (CS ATT(g), ATT(e), 41 pool) --")
mp <- did::att_gt(yname=Y, tname=TIME, idname=UNITN, gname=COHORT,
        xformla=~1, control_group="nevertreated", clustervars=UNITN, data=dta)
agg_g <- did::aggte(mp, type="group", na.rm=TRUE)
agg_e <- did::aggte(mp, type="dynamic", na.rm=TRUE)
attg <- data.frame(kind="ATT(g)", key=agg_g$egt, estimate=agg_g$att.egt,
                   se=agg_g$se.egt)
atte <- data.frame(kind="ATT(e)", key=agg_e$egt, estimate=agg_e$att.egt,
                   se=agg_e$se.egt)
attg$ci_low<-attg$estimate-Zc*attg$se; attg$ci_high<-attg$estimate+Zc*attg$se
atte$ci_low<-atte$estimate-Zc*atte$se; atte$ci_high<-atte$estimate+Zc*atte$se
share_g_neg <- mean(attg$estimate < 0)
share_g_sig <- mean(abs(attg$estimate/attg$se) > Zc)
lg(sprintf("R6: %d cohort ATT(g); share negative = %.3f; share individually sig(5%%) = %.3f",
   nrow(attg), share_g_neg, share_g_sig))
lg(paste(capture.output(print(attg, row.names=FALSE, digits=4)), collapse="\n"))

# =============================================================================
# atoms_harmonized.csv  (R2h, R3h, R6)
# =============================================================================
mk <- function(atoms, spec) data.frame(spec=spec, atom_id=atoms$atom_id,
  group=atoms$group, time=atoms$time, estimate=atoms$estimate, se=atoms$se,
  weight=atoms$weight, stringsAsFactors=FALSE)
r6_atoms <- rbind(
  data.frame(spec="R6_ATTg", atom_id=paste0("g",attg$key), group=attg$key,
             time=NA, estimate=attg$estimate, se=attg$se, weight=NA),
  data.frame(spec="R6_ATTe", atom_id=paste0("e",atte$key), group=NA,
             time=atte$key, estimate=atte$estimate, se=atte$se, weight=NA))
atoms_h <- rbind(mk(st_R2h,"R2h"), mk(st_R3h,"R3h"), r6_atoms)
write.csv(atoms_h, file.path(OUT_DIR,"atoms_harmonized.csv"), row.names=FALSE)
lg("\nWrote atoms_harmonized.csv (", nrow(atoms_h), "rows )")

# =============================================================================
# ladder.csv  (existing tabulated + new)
# =============================================================================
row <- function(id,label,pool,nctrl,dim,est,se,nobs,nunits,ncoh,src) {
  ci_l <- if(is.finite(se)) est-Zc*se else NA
  ci_h <- if(is.finite(se)) est+Zc*se else NA
  data.frame(id=id,label=label,control_pool=pool,n_controls=nctrl,
    changed_dimension=dim,estimate=est,se=se,ci_low=ci_l,ci_high=ci_h,
    n_obs=nobs,n_units=nunits,n_cohorts=ncoh,source=src,
    supports_pp=supports(est,se),stringsAsFactors=FALSE)
}
L <- list()
# --- existing (tabulated from output files) ---
L[[length(L)+1]] <- row("R0","P&P Table 3 m4: full-panel stack, shared FE, cov+ebal wts",
  "741 No-Change",741,"reference (paper spec)",-0.0908090808214521,0.0043288708766928,
  171444,NA,NA,"output/cengiz_control_groups_summary.csv; assumptions.md §8,§14")
L[[length(L)+1]] <- row("R1","R0 plain (no cov, no wts), shared agency FE",
  "741 No-Change",741,"drop cov+wts",-0.104397217091685,0.0382690468295203,
  171906,NA,NA,"output/cengiz_control_groups_summary.csv; assumptions.md §8")
L[[length(L)+1]] <- row("R2","R1 with interacted agency×stack FE",
  "741 No-Change",741,"interacted FE",R2_est,R2_se,
  R2_sum$n_obs,R2_sum$n_units,length(R2_cohorts),
  "output/recombination_check.csv (est); computed here (se,N)")
L[[length(L)+1]] <- row("R3","Cengiz ±4, never-treated","741 No-Change",741,
  "windowed ±4",-0.0615125201022484,0.0692438762505653,NA,NA,NA,
  "output/cengiz_control_groups_summary.csv")
L[[length(L)+1]] <- row("R4","Cengiz ±4, not-yet-treated","741+not-yet",NA,
  "not-yet controls",-0.0616300138652742,0.0692119598826385,NA,NA,NA,
  "output/cengiz_control_groups_summary.csv")
L[[length(L)+1]] <- row("R5","CS=SA never-treated (41), no cov, overall ATT",
  "41 always-0",41,"CS estimator + 41 pool",R5_est,R5_se,
  NA,NA,R5_ncohort,"output/cs_variants_summary.csv (est); computed here")
# --- new ---
L[[length(L)+1]] <- row("R2h","Stacked interacted FE, harmonized 41-pool",
  "41 always-0",41,"control pool 741→41",R2h_est,R2h_se,
  R2h_sum$n_obs,R2h_sum$n_units,length(R2h_cohorts),"computed here")
L[[length(L)+1]] <- row("R3h","Cengiz ±4, harmonized 41-pool","41 always-0",41,
  "control pool 741→41",R3h_est,R3h_se,nrow(cz_h_est),
  length(unique(cz_h_est$agency.id)),length(R3h_cohorts),"computed here")
L[[length(L)+1]] <- row("R5d","CS FULL 741 pool [DIAGNOSTIC ONLY — false history]",
  "741 (always-1 recoded g=0)",741,"CS + 741 pool",R5d_est,R5d_se,
  R5d_nobs,NA,R5d_ncohort,"computed here (DECOMPOSITION DEVICE, not a candidate spec)")
L[[length(L)+1]] <- row("R5b_R2","CS-simple weights × R2 atoms (741 pool)",
  "741 No-Change",741,"CS weights on stacked atoms",r5b_R2$value,NA,
  R2_sum$n_obs,NA,length(r5b_R2$matched),"computed here")
L[[length(L)+1]] <- row("R5b_R2h","CS-simple weights × R2h atoms (41 pool)",
  "41 always-0",41,"CS weights on stacked atoms",r5b_R2h$value,NA,
  R2h_sum$n_obs,NA,length(r5b_R2h$matched),"computed here")
L[[length(L)+1]] <- row("R6","CS ATT(g)/ATT(e) disaggregated (41 pool)","41 always-0",41,
  "disaggregated targets",R5_est,R5_se,NA,NA,nrow(attg),
  "computed here (see atoms_harmonized.csv; overall=R5)")

ladder <- do.call(rbind, L)
write.csv(ladder, file.path(OUT_DIR,"ladder.csv"), row.names=FALSE)
lg("\nWrote ladder.csv (", nrow(ladder), "rows )")
lg(paste(capture.output(print(ladder[,c("id","control_pool","n_controls","estimate","se","n_obs","n_cohorts","supports_pp")], row.names=FALSE, digits=5)), collapse="\n"))

# stash R5b flag + R6 shares for the report
lg(sprintf("\nR6 share_negative=%.3f share_sig=%.3f", share_g_neg, share_g_sig))
lg("\n================ Stage L sessionInfo ================")
lg(paste(capture.output(print(sessionInfo())), collapse="\n"))
cat("\nDONE Stage L\n")
