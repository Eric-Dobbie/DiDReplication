# =============================================================================
# design_sweep_core.R  --  DETERMINISTIC design sweep over panel geometry
#
# NOT a simulation. There is NO outcome variable, NO DGP, NO sampling of data,
# and NO replications. Aggregation weights are EXACT deterministic functions of
# the design (panel length T, adoption-period vector g, cohort sizes n_g, control
# pool N_c). This file defines design_weights(); every figure/validation script
# is this function enumerated over design configurations.
#
# Per cohort it returns: pre_g, post_g, k = pre_g - 1; CS weights (simple, group,
# dynamic, calendar); stacked weight w_s EXACT (two-way within-stack demeaning of
# the treatment indicator -> V_s/V) and the closed-form approximation
# N_s*T_s*p(1-p)*q(1-q); k/n_g; and identification status per estimator.
#
# Periods are indexed 1..T. A cohort adopts at period gp in 2..T (needs >=1 pre).
# pre_g = gp-1 ; post_g = T-gp+1 (periods t>=gp) ; k = pre_g-1 (joint pre-test
# restrictions, from B1). Control pool N_c is common to all stacks (never-treated).
# =============================================================================

# ---- CS cohort weights (exact; match did::aggte, verified vs P&P below) ------
cs_weights <- function(gp, n_g, T) {
  G <- length(gp); post <- T - gp + 1
  simple <- n_g * post;                         simple <- simple/sum(simple)
  group  <- n_g;                                group  <- group /sum(group)
  ming <- min(gp)
  # dynamic: average over event times e=0..(T-ming); at e, cohorts with gp<=T-e share by n_g
  Edyn <- 0:(T - ming); wdyn <- numeric(G)
  for (e in Edyn) { elig <- which(gp <= T - e); den <- sum(n_g[elig])
    wdyn[elig] <- wdyn[elig] + (n_g[elig]/den) }
  wdyn <- (wdyn/length(Edyn)); wdyn <- wdyn/sum(wdyn)
  # calendar: average over calendar periods t=ming..T; at t, cohorts with gp<=t share by n_g
  Tcal <- ming:T; wcal <- numeric(G)
  for (t in Tcal) { elig <- which(gp <= t); den <- sum(n_g[elig])
    wcal[elig] <- wcal[elig] + (n_g[elig]/den) }
  wcal <- (wcal/length(Tcal)); wcal <- wcal/sum(wcal)
  list(simple=simple, group=group, dynamic=wdyn, calendar=wcal)
}

# ---- stacked V_s EXACT: within-stack two-way demean of the treatment indicator
# balanced stack: units 1..(n_g+N_c) x periods 1..T; D[i,t]=1 iff i<=n_g & t>=gp.
# Two-way FE residual (exact for balanced) = D - rowmean - colmean + grandmean.
Vs_exact <- function(gp, n_g, N_c, T) {
  Nu <- n_g + N_c
  if (N_c == 0 && n_g <= 1) return(0)           # no contrast -> collinear
  if (gp > T || gp < 2) return(0)               # no post or no pre -> absorbed
  D <- matrix(0, Nu, T); if (n_g >= 1) D[1:n_g, gp:T] <- 1
  R <- D - rowMeans(D) - rep(colMeans(D), each=Nu) + mean(D)
  sum(R^2)
}
# ---- closed-form approximation ----------------------------------------------
Vs_approx <- function(gp, n_g, N_c, T) {
  Nu <- n_g + N_c; p <- n_g/Nu; q <- (T - gp + 1)/T
  Nu * T * p*(1-p) * q*(1-q)
}

# ---- master: all weights + identification for a design -----------------------
design_weights <- function(T, gp, n_g, N_c) {
  stopifnot(length(gp)==length(n_g))
  G <- length(gp); pre <- gp - 1; post <- T - gp + 1; k <- pre - 1
  cs <- cs_weights(gp, n_g, T)
  Vs  <- mapply(Vs_exact,  gp, n_g, MoreArgs=list(N_c=N_c, T=T))
  Vsa <- mapply(Vs_approx, gp, n_g, MoreArgs=list(N_c=N_c, T=T))
  ws  <- if (sum(Vs)>0)  Vs/sum(Vs)   else rep(NA_real_, G)
  wsa <- if (sum(Vsa)>0) Vsa/sum(Vsa) else rep(NA_real_, G)
  # identification per cohort
  id_cs <- ifelse(pre < 1, "collinear:no_pre", ifelse(n_g==1, "single_treated", "ok"))
  id_st <- ifelse(N_c==0 & n_g<=1, "no_contrast",
            ifelse(gp>T | gp<2, "collinear:FE",
             ifelse(Vs==0, "collinear:FE", ifelse(n_g==1, "single_treated", "ok"))))
  ws[Vs==0] <- NA_real_
  data.frame(cohort=seq_len(G), gp=gp, n_g=n_g, pre=pre, post=post, k=k, k_over_n=k/n_g,
    w_simple=cs$simple, w_group=cs$group, w_dynamic=cs$dynamic, w_calendar=cs$calendar,
    V_s=Vs, w_stacked=ws, V_s_approx=Vsa, w_stacked_approx=wsa,
    id_CS=id_cs, id_stacked=id_st, stringsAsFactors=FALSE)
}

# =============================================================================
# SELF-TEST (runs when sourced directly): CS weights vs committed P&P weights,
# and stacked exact demean vs an felm residualization on one synthetic stack.
# =============================================================================
if (sys.nframe()==0 || identical(environment(), globalenv())) {
 .selftest <- function() {
  OUT <- Sys.getenv("DIDREP_OUT", unset="output")
  # P&P CS cohorts (period index, 2000=1): years 2001..2020
  yrs <- c(2001,2002,2003,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
  gp  <- yrs - 2000 + 1; n_g <- c(1,5,1,2,1,7,1,3,1,3,3,1,1); T <- 21
  cs <- cs_weights(gp, n_g, T)
  ref <- read.csv(file.path(OUT,"task2_cohort_weights.csv"))
  cat("CS weight self-test vs task2_cohort_weights.csv (max abs diff):\n")
  for (nm in c("simple","group","dynamic","calendar"))
    cat(sprintf("  %-9s %.2e\n", nm, max(abs(cs[[nm]] - ref[[nm]]))))
  # stacked exact demean vs felm on a synthetic balanced stack
  suppressMessages(library(lfe))
  gps<-5; ngs<-3; Ncs<-40; Ts<-21; Nu<-ngs+Ncs
  d <- expand.grid(i=1:Nu, t=1:Ts); d$D <- as.integer(d$i<=ngs & d$t>=gps)
  r <- felm(D ~ 1 | factor(i) + factor(t), data=d)$residuals[,1]
  cat(sprintf("stacked V_s self-test: demean %.6f vs felm %.6f (diff %.2e)\n",
      Vs_exact(gps,ngs,Ncs,Ts), sum(r^2), abs(Vs_exact(gps,ngs,Ncs,Ts)-sum(r^2))))
 }
 if (!exists(".sourced_quiet")) .selftest()
}
