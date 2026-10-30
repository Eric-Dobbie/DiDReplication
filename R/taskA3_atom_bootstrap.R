# =============================================================================
# taskA3_atom_bootstrap.R
#
# Is the Task-1 atom channel (-0.117) a real identification difference or small-
# cohort sampling noise? For the 11 cohorts identified under both estimators,
# compare ATT_CS(g) and beta_stacked(g) with a PAIRED agency-cluster bootstrap
# that resamples agencies and recomputes BOTH estimators on each draw (the two
# estimates share data, so an analytic combination of separate SEs is wrong).
#
# Fast per-draw estimators (verified below against the committed atoms):
#   ATT_CS(g)      = mean_{t>=g} { [Ybar_{g,t}-Ybar_{g,g-1}] - [Ybar_{nt,t}-Ybar_{nt,g-1}] }
#                    never-treated controls, base g-1 (post ATT is base-invariant).
#   beta_stk(g)    = feols(any.fatalities ~ no.req | agency + year) on cohort g's
#                    shipped stack (within one stack this IS the interacted-FE atom).
# Aggregation weights are held at their full-sample values (WCS simple, w_s FWL):
#   the question is atom-noise, not weight-noise. Atom channel = sum_g WCS(g)*
#   (beta_stk(g) - ATT_CS(g)).
# Writes output/taskA3_atom_bootstrap.csv (per-cohort) and prints the atom-channel
# interval. B draws, seeded.
# =============================================================================
suppressMessages({library(fixest)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
B <- as.integer(Sys.getenv("A3_B", unset = "500"))

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$yc <- ifelse(is.na(dta$year.changed) & dta$no.req==0, 0,
          ifelse(is.na(dta$year.changed) & dta$no.req==1, 1987, dta$year.changed))
st <- read.csv(file.path(DATA_DIR, "stacked_fatal.csv"))

cohorts <- c(2002,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
# full-sample weights (CS simple cohort weights + stacked FWL), from committed inputs
ci <- read.csv(file.path(OUT_DIR, "task1_cohort_inputs.csv"))
WCS <- setNames(ci$WCS_restr, ci$cohort)[as.character(cohorts)]

# ---- fast CS cohort effect on a panel frame d (cols: agency.id,id,year,yc,any.fatalities)
att_cs <- function(d, g) {
  base <- g - 1
  nt <- d[d$yc == 0, ]; tr <- d[d$yc == g, ]
  if (!nrow(tr)) return(NA_real_)
  ybar <- function(df, yr) {                       # mean outcome at year yr by picking rows
    v <- df$any.fatalities[df$year == yr]; mean(v, na.rm=TRUE) }
  ntb <- ybar(nt, base); trb <- ybar(tr, base)
  ts <- g:2020
  mean(sapply(ts, function(t) (ybar(tr,t)-trb) - (ybar(nt,t)-ntb)))
}
# ---- fast stacked atom on cohort g's stack subset (resampled) ----------------
beta_stk <- function(s) {
  m <- tryCatch(fixest::feols(any.fatalities ~ no.req | agency.id2 + year,
                              data = s, warn=FALSE, notes=FALSE), error=function(e) NULL)
  if (is.null(m)) return(NA_real_)
  cf <- tryCatch(coef(m), error=function(e) NULL)
  if (is.null(cf) || !("no.req" %in% names(cf))) return(NA_real_)
  v <- unname(cf["no.req"]); if (length(v)!=1 || is.na(v)) NA_real_ else v
}

# ---- full-sample point estimates (verify vs committed) ----------------------
ntr <- sapply(cohorts, function(g) length(unique(dta$agency.id[dta$yc==g])))
A_CS0 <- sapply(cohorts, function(g) att_cs(dta, g))
st$agency.id2 <- st$agency.id
A_ST0 <- sapply(cohorts, function(g){ s <- st[st$cohort==g,]; beta_stk(s) })
cat("== full-sample check vs committed atoms ==\n")
chk <- data.frame(cohort=cohorts, A_CS_fast=round(A_CS0,4), A_CS_committed=round(ci$A_CS,4),
                  A_ST_fast=round(A_ST0,4), A_ST_committed=round(ci$A_stacked,4))
print(chk, row.names=FALSE)
atom_channel0 <- sum(WCS * (A_ST0 - A_CS0))
cat(sprintf("\nfull-sample atom channel = %.5f (Task 1 reported -0.117)\n\n", atom_channel0))

# ---- paired agency-cluster bootstrap (vectorized) ---------------------------
ag <- unique(dta$agency.id)
OBS <- as.character(2000:2020)
# never-treated outcome matrix [agency x year]
mkmat <- function(ids) {
  m <- matrix(NA_real_, length(ids), length(OBS), dimnames=list(ids, OBS))
  sub <- dta[dta$agency.id %in% ids, ]
  m[cbind(as.character(sub$agency.id), as.character(sub$year))[sub$year>=2000,]] <-
    sub$any.fatalities[sub$year>=2000]
  m
}
nt_ids <- unique(dta$agency.id[dta$yc==0]); Ynt <- mkmat(nt_ids)
tr_ids <- lapply(cohorts, function(g) unique(dta$agency.id[dta$yc==g])); names(tr_ids)<-cohorts
Ytr <- lapply(cohorts, function(g) mkmat(tr_ids[[as.character(g)]])); names(Ytr)<-cohorts
# stacked: per-cohort stack df + row indices by agency
stg_df <- lapply(cohorts, function(g) st[st$cohort==g, c("no.req","year","any.fatalities","agency.id")])
names(stg_df)<-cohorts
agidx <- lapply(cohorts, function(g) split(seq_len(nrow(stg_df[[as.character(g)]])),
                                           stg_df[[as.character(g)]]$agency.id)); names(agidx)<-cohorts

wmean <- function(M, cnt) {                         # weighted column means, weights=counts
  w <- cnt[rownames(M)]; w[is.na(w)] <- 0
  if (sum(w)==0) return(setNames(rep(NA_real_,ncol(M)), colnames(M)))
  colSums(M*w, na.rm=TRUE)/colSums((!is.na(M))*w)
}
att_cs_fast <- function(ntY, trY, g) {
  bk <- as.character(g-1); if (is.na(ntY[bk]) || is.na(trY[bk])) return(NA_real_)
  tt <- as.character(g:2020); mean((trY[tt]-trY[bk]) - (ntY[tt]-ntY[bk]), na.rm=TRUE)
}
set.seed(1)
BCS <- matrix(NA_real_, B, length(cohorts)); BST <- BCS
BCH <- rep(NA_real_, B); complete <- logical(B)
for (b in 1:B) {
  draw <- sample(ag, length(ag), replace=TRUE)
  cnt <- table(draw); cntv <- setNames(as.integer(cnt), names(cnt))
  ntY <- wmean(Ynt, cntv)
  a_cs <- numeric(length(cohorts)); a_st <- numeric(length(cohorts))
  for (j in seq_along(cohorts)) { g <- cohorts[j]
    a_cs[j] <- att_cs_fast(ntY, wmean(Ytr[[j]], cntv), g)
    il <- agidx[[j]][draw]; present <- !vapply(il, is.null, logical(1))
    if (!any(present)) { a_st[j] <- NA_real_; next }
    rows <- unlist(il[present], use.names=FALSE)
    ids  <- rep.int(which(present), lengths(il[present]))
    s <- stg_df[[j]][rows,]; s$agency.id2 <- ids
    a_st[j] <- beta_stk(s)
  }
  BCS[b,] <- a_cs; BST[b,] <- a_st
  ok <- !is.na(a_cs) & !is.na(a_st); complete[b] <- all(ok)
  BCH[b] <- if (any(ok)) sum(WCS[ok]*(a_st[ok]-a_cs[ok])) else NA_real_   # partial sum; dropped cohorts are low-weight (0.003-0.12)
}

# ---- per-cohort summary -----------------------------------------------------
res <- data.frame(cohort=cohorts, n_treated=ntr,
  A_CS=round(A_CS0,4), se_A_CS=round(apply(BCS,2,sd,na.rm=TRUE),4),
  A_ST=round(A_ST0,4), se_A_ST=round(apply(BST,2,sd,na.rm=TRUE),4),
  diff=round(A_CS0-A_ST0,4))
dmat <- BCS - BST
res$se_diff <- round(apply(dmat,2,sd,na.rm=TRUE),4)
res$p_diff <- round(mapply(function(j) 2*min(mean(dmat[,j]<=0,na.rm=TRUE),
                                             mean(dmat[,j]>=0,na.rm=TRUE)), seq_along(cohorts)),3)
res$eff_draws <- apply(!is.na(dmat),2,sum)
res$distinguishable <- ifelse(res$p_diff < 0.05, "YES","")
options(width=200); print(res, row.names=FALSE)
cat(sprintf("\ncohorts with difference distinguishable from 0 (p<0.05): %d of %d\n",
            sum(res$distinguishable=="YES"), nrow(res)))
qc  <- quantile(BCH, c(.025,.5,.975), na.rm=TRUE)
BCHc <- BCH[complete]; qcc <- quantile(BCHc, c(.025,.5,.975), na.rm=TRUE)
cat(sprintf("atom channel (all draws, partial sum over estimable cohorts): point %.4f | mean %.4f sd %.4f | 95%% CI [%.4f, %.4f]\n",
            atom_channel0, mean(BCH,na.rm=TRUE), sd(BCH,na.rm=TRUE), qc[1], qc[3]))
cat(sprintf("atom channel (complete draws only, n=%d of %d): mean %.4f sd %.4f | 95%% CI [%.4f, %.4f]\n",
            sum(complete), B, mean(BCHc,na.rm=TRUE), sd(BCHc,na.rm=TRUE), qcc[1], qcc[3]))
write.csv(res, file.path(OUT_DIR,"taskA3_atom_bootstrap.csv"), row.names=FALSE)
write.csv(data.frame(draw=1:B, atom_channel=BCH), file.path(OUT_DIR,"taskA3_atom_channel_draws.csv"), row.names=FALSE)
cat("\nWrote output/taskA3_atom_bootstrap.csv and taskA3_atom_channel_draws.csv\n")
