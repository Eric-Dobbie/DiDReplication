# =============================================================================
# design_sweep_validation.R  --  closed form vs exact V_s (deterministic).
# Derivation: for a BALANCED stack with a control pool of all-zero controls common
# to the stack, the two-way-demeaned treatment residual factorizes,
#   resid_{it} = (a_i - pbar)(b_t - qbar),  a_i=1[treated], b_t=1[post],
# so V_s = [N p(1-p)] * [T q(1-q)] EXACTLY. This script confirms that across the
# sweep (deviation ~ machine precision) and measures where it degrades:
#   (1) panel imbalance (missing cells), (2) stack-specific vs common control pool,
#   (3) balancing/sampling weights. Then applies it to P&P's 11 stacks.
# =============================================================================
suppressMessages({library(lfe)})
.sourced_quiet <- TRUE; source("R/design_sweep_core.R")
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); DATA <- Sys.getenv("DIDREP_DATA", unset="data")

# exact V_s via felm on an arbitrary (possibly unbalanced/weighted) stack frame
Vs_felm <- function(d, w=NULL) {
  if (is.null(w)) r <- felm(D ~ 1 | iid + tid, data=d)$residuals[,1]
  else            r <- felm(D ~ 1 | iid + tid, data=d, weights=w)$residuals[,1]
  if (is.null(w)) sum(r^2) else sum(w*r^2)
}
build_stack <- function(gp, n_g, N_c, T) {
  d <- expand.grid(i=1:(n_g+N_c), t=1:T); d$D <- as.integer(d$i<=n_g & d$t>=gp)
  d$iid <- factor(d$i); d$tid <- factor(d$t); d
}

# ---- (0) sweep: balanced common-pool -> exact should equal approx ------------
set.seed(3); dev <- c()
for (it in 1:200) {
  T <- sample(12:30,1); gp <- sample(2:T,1); n_g <- sample(1:8,1); N_c <- sample(20:60,1)
  ex <- Vs_exact(gp,n_g,N_c,T); ap <- Vs_approx(gp,n_g,N_c,T)
  if (ex>0) dev <- c(dev, abs(ex-ap)/ex)
}
cat(sprintf("(0) BALANCED common-pool: max relative |exact-approx|/exact over 200 designs = %.2e  => closed form is EXACT\n", max(dev)))

# ---- (1) imbalance: drop a share of cells, exact via felm vs approx ----------
cat("\n(1) panel imbalance (random missing cells):\n")
for (miss in c(0,0.05,0.10,0.20,0.40)) {
  set.seed(5); rr <- c()
  for (it in 1:60) {
    T<-21; gp<-sample(3:18,1); n_g<-sample(2:8,1); N_c<-40
    d <- build_stack(gp,n_g,N_c,T)
    if (miss>0) d <- d[sample(nrow(d), round(nrow(d)*(1-miss))),]
    ex <- Vs_felm(d); ap <- Vs_approx(gp,n_g,N_c,T)   # approx uses nominal balanced p,q
    if (ex>0) rr <- c(rr, abs(ex-ap)/ex)
  }
  cat(sprintf("   missing share %.2f: median rel deviation %.3f, max %.3f\n", miss, median(rr), max(rr)))
}

# ---- (2) weights: downweight controls (balancing-weight analogue) ------------
cat("\n(2) balancing weights (controls downweighted):\n")
for (wc in c(1, 0.5, 0.1)) {
  set.seed(9); rr<-c()
  for (it in 1:60) { T<-21; gp<-sample(3:18,1); n_g<-sample(2:8,1); N_c<-40
    d <- build_stack(gp,n_g,N_c,T); w <- ifelse(d$i<=n_g, 1, wc)
    ex <- Vs_felm(d, w); ap <- Vs_approx(gp,n_g,N_c,T)
    if (ex>0) rr<-c(rr, abs(ex-ap)/ex) }
  cat(sprintf("   control weight %.2f: median rel deviation %.3f\n", wc, median(rr)))
}

# ---- apply to P&P: exact V_s (fwl_decomp) vs closed form ---------------------
cat("\n(P&P) exact V_s/V (from fwl_decomp) vs closed-form prediction:\n")
fw <- read.csv(file.path(OUT,"fwl_decomp_unweighted.csv")); fw <- fw[fw$w_s>0 & !is.na(fw$beta_s),]
st <- read.csv(file.path(DATA,"stacked_fatal.csv")); est <- st[!is.na(st$any.fatalities),]
Tpp <- length(unique(est$year))                                   # 21
pp <- do.call(rbind, lapply(fw$stack, function(g){
  s <- est[est$cohort==g,]; n_g <- length(unique(s$agency.id[s$treat==1])); N_c <- length(unique(s$agency.id[s$treat==0]))
  post <- length(unique(s$year[s$year>=g]))
  data.frame(stack=g, n_g=n_g, N_c=N_c, post=post, Vs_exact=fw$V_s[fw$stack==g],
             Vs_cf=Vs_approx(gp=g-2000+1 + (Tpp - length(unique(s$year))), n_g, N_c, Tpp)) }))
# simpler closed form using post share directly (T=21):
pp$p <- pp$n_g/(pp$n_g+pp$N_c); pp$q <- pp$post/Tpp
pp$Vs_cf <- (pp$n_g+pp$N_c)*Tpp*pp$p*(1-pp$p)*pp$q*(1-pp$q)
pp$w_exact <- pp$Vs_exact/sum(pp$Vs_exact); pp$w_cf <- pp$Vs_cf/sum(pp$Vs_cf)
print(transform(pp[,c("stack","n_g","N_c","post","Vs_exact","Vs_cf","w_exact","w_cf")],
      Vs_exact=round(Vs_exact,3), Vs_cf=round(Vs_cf,3), w_exact=round(w_exact,3), w_cf=round(w_cf,3)), row.names=FALSE)
cat(sprintf("\nP&P w_exact vs w_cf: correlation %.3f, L1 %.3f\n",
    cor(pp$w_exact,pp$w_cf), sum(abs(pp$w_exact-pp$w_cf))))
r_exact <- pp$Vs_exact[pp$stack==2009]/pp$Vs_exact[pp$stack==2013]
r_cf    <- pp$Vs_cf[pp$stack==2009]/pp$Vs_cf[pp$stack==2013]
cat(sprintf("2009-vs-2013 V_s ratio: closed-form %.2f vs exact %.2f\n", r_cf, r_exact))
write.csv(pp, file.path(OUT,"design_validation_pp.csv"), row.names=FALSE)
cat("\nWrote output/design_validation_pp.csv\n")
