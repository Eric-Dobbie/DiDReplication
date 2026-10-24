# =============================================================================
# design_sweep_checkability.R  --  weight-vs-checkability (analytic MDE) across
# the design space. Sweep version of weight_checkability.R / Figure C.
# REVISION (supersedes commit 55ddc68): units, pools, and headline corrected.
#
# DETERMINISTIC. No outcome variable, no DGP, no sampling of data. For every
# design configuration (panel length T, adoption-period vector g, cohort sizes
# n_g, never-treated count N_c) and every aggregation scheme we compute, per
# cohort: its aggregation weight (exact, design_sweep_core.R); whether its joint
# pre-trend test is CHECKABLE (k>=1 pre-leads AND k/n_g<=1); and, if checkable,
# the analytic MINIMUM DETECTABLE pre-trend at 80% power / alpha 0.05, in sigma
# units.
#
# UNITS (fix 1). The treatment effect is a LEVEL, so the comparison is level-on-
# level. We report two level measures, NOT the per-period slope:
#   - JUMP MDE  : smallest detectable uniform level pre-trend  (a = rep(1,k)):
#                 MDE_jump/sigma = sqrt(lambda*(k)/[k/(k+1)]) * sqrt(1/n_g+1/n_c).
#   - TREND-BIAS: slope MDE converted to implied mean post bias = slope * h,
#                 h = mean post event time = (T - gp)/2 (balanced).
# "Resolvable" = checkable AND MDE_jump/sigma <= 0.303 ( = |beta|=0.10 at
# sigma=0.3297, a P&P-scale fatal-encounters effect). Jump is the headline level
# measure (a persistent level pre-trend biases DiD by its full size regardless of
# horizon, so trend-bias -> 0 as horizon -> 0 is avoided).
#
# POOLS (fix 4) come from each configuration's geometry, not fixed at 41/741:
#   - CS family (simple/group/dynamic/calendar) + SA: NEVER-TREATED count N_c.
#   - stacked (and CS-NYT): NOT-YET-TREATED + never-treated at g =
#       N_c + sum_{g': gp'>gp} n_{g'}.  (In the data the 741 pool also folds in
#       700 ALWAYS-treated/no-requirement units, which are not valid not-yet-
#       treated controls; a geometry-honest pool excludes them -- reported.)
# The stacked weight V_s is computed with each stack's own pool.
#
# HEADLINE (fix 7): the weight-weighted MDE distribution, since "resolvable
# share" is degenerate (0 almost everywhere) after fix 1.
#
# P&P ANCHOR (fix 5): generated with the SAME grid functions on P&P's timing
# vector + CS cohort sizes + N_c=41 (so it is comparable and can be placed in the
# swept distribution). The committed real-data point (741 pool, from
# weight_checkability_cohorts.csv) is shown only as a separate reference marker --
# no percentile is computed for it.
# =============================================================================
suppressMessages({library(ggplot2)})
.sourced_quiet <- TRUE
source("R/design_sweep_core.R")
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); FIG <- Sys.getenv("DIDREP_FIG", unset="figures")

SIGMA   <- 0.3297
BETA_PP <- 0.10
THR     <- BETA_PP/SIGMA                 # resolvable line on MDE/sigma (0.303)
SCHEMES <- c("simple","group","dynamic","calendar","stacked")
CS_FAM  <- c("simple","group","dynamic","calendar")

# ---- analytic MDE/sigma: lambda*(k), denom per shape, precomputed per k -------
lambda_star <- function(k){ if (k<1) return(NA_real_); crit <- qchisq(0.95, k)
  uniroot(function(l) pchisq(crit, df=k, ncp=l, lower.tail=FALSE)-0.80, c(1e-6, 400), tol=1e-10)$root }
KMAX <- 40; kk <- 1:KMAX
den_jump  <- kk/(kk+1)                                   # a'(I+11')^-1 a, a=1
den_trend <- sapply(kk, function(k){ a<-1:k; sum(a^2) - sum(a)^2/(k+1) })
lam  <- sapply(kk, lambda_star)
Cjump  <- sqrt(lam/den_jump)                             # MDE_jump/sigma  = Cjump[k]*sqrt(1/n_g+1/n_c)
Ctrend <- sqrt(lam/den_trend)                            # slope MDE/sigma = Ctrend[k]*sqrt(1/n_g+1/n_c)
mde_jump_sig  <- function(k,n_g,n_c){ ok<-!is.na(k)&k>=1&k<=KMAX&n_c>0&n_g>0
  o<-rep(NA_real_,length(k)); o[ok]<-Cjump[k[ok]]*sqrt(1/n_g[ok]+1/n_c[ok]); o }
mde_tbias_sig <- function(k,n_g,n_c,h){ ok<-!is.na(k)&k>=1&k<=KMAX&n_c>0&n_g>0
  o<-rep(NA_real_,length(k)); o[ok]<-Ctrend[k[ok]]*sqrt(1/n_g[ok]+1/n_c[ok])*h[ok]; o }

# ---- geometry-derived control pools for a design -----------------------------
pools_of <- function(gp, n_g, N_c){           # returns n_c per cohort, per family
  nyt <- sapply(seq_along(gp), function(i) sum(n_g[gp > gp[i]]))   # strictly-later adopters
  list(cs = rep(N_c, length(gp)), stacked = N_c + nyt) }

# ---- stacked weights with per-cohort pool ------------------------------------
stacked_w <- function(gp, n_g, T, nc_vec){
  Vs <- mapply(function(g,ng,nc) Vs_exact(g, ng, nc, T), gp, n_g, nc_vec)
  w <- if (sum(Vs)>0) Vs/sum(Vs) else rep(NA_real_, length(gp)); w[Vs==0] <- NA_real_; w }

# ---- binding reason + shares for one design, one scheme ----------------------
#   reasons: no_pre (k<1) | size (k/n_g>1) | underpowered (jump MDE>THR) | resolvable
classify <- function(k, n_g, n_c, w, size_strict=FALSE){
  mj <- mde_jump_sig(k, n_g, n_c)
  size <- !is.na(k) & k>=1 & (if (size_strict) k/n_g >= 1 else k/n_g > 1)
  reason <- ifelse(is.na(k)|k<1, "no_pre", ifelse(size, "size",
             ifelse(is.finite(mj) & mj>THR, "underpowered", "resolvable")))
  list(reason=reason, mj=mj) }
shares <- function(reason, mj, w){
  ident <- !is.na(w); tot <- sum(w[ident]); if (tot<=0) return(NULL)
  sh <- function(r) sum(w[ident & reason==r])/tot
  ok <- ident & reason %in% c("underpowered","resolvable") & is.finite(mj)
  wmed <- if (any(ok)){ o<-order(mj[ok]); cw<-cumsum(w[ok][o])/sum(w[ok]); mj[ok][o][which(cw>=0.5)[1]] } else NA_real_
  c(no_pre=sh("no_pre"), size=sh("size"), underpowered=sh("underpowered"),
    resolvable=sh("resolvable"), wmed_jump=wmed, n_ident=sum(ident)) }

# ---- design generators (same universe as design_sweep_figureC.R) -------------
place <- function(T, G, type){ slots <- 2:T
  if (type=="early") head(slots,G) else if (type=="late") tail(slots,G)
  else if (type=="spread") unique(round(seq(2,T,length.out=G))) else sort(sample(slots,G)) }
sizes <- function(G, type){ if (type=="equal") rep(3L,G)
  else { s<-rep(1L,G); s[1]<-8L; if(G>=2) s[2]<-5L; if(G>=3) s[3]<-3L; s } }

grid_rows <- list(); coh_rows <- list()
addcfg <- function(T,G,pl,sz,Nc,arm,cfgid){
  gp <- place(T,G,pl); if (length(gp)<G) return(invisible())
  ng <- sizes(length(gp), sz); k <- gp-2; h <- (T-gp)/2
  cw <- cs_weights(gp, ng, T); pl_nc <- pools_of(gp, ng, Nc)
  ws <- stacked_w(gp, ng, T, pl_nc$stacked)
  wlist <- list(simple=cw$simple, group=cw$group, dynamic=cw$dynamic, calendar=cw$calendar, stacked=ws)
  for (s in SCHEMES){
    nc <- if (s=="stacked") pl_nc$stacked else pl_nc$cs
    cl <- classify(k, ng, nc, wlist[[s]]); sh <- shares(cl$reason, cl$mj, wlist[[s]])
    if (is.null(sh)) next
    grid_rows[[length(grid_rows)+1L]] <<- data.frame(cfg=cfgid, T=T, G=length(gp), placement=pl,
      size=sz, N_c=Nc, arm=arm, scheme=s, no_pre=sh["no_pre"], size_fail=sh["size"],
      underpowered=sh["underpowered"], resolvable=sh["resolvable"], wmed_jump=sh["wmed_jump"],
      n_ident=sh["n_ident"], row.names=NULL)
    if (arm=="grid" && Nc==41)
      coh_rows[[length(coh_rows)+1L]] <<- data.frame(cfg=cfgid, T=T, placement=pl, size=sz,
        scheme=s, cohort=gp, k=k, n_g=ng, n_c=nc, weight=wlist[[s]], reason=cl$reason,
        mde_jump_sig=cl$mj, row.names=NULL) }
}
cfg <- 0L
for (T in c(15,21,30)) for (G in c(3,5,8,11)) for (pl in c("early","late","spread"))
  for (sz in c("equal","skewed")) for (Nc in c(10,41,120)) if (G<=T-1){ cfg<-cfg+1L; addcfg(T,G,pl,sz,Nc,"grid",cfg) }
set.seed(11)
for (T in 21) for (G in c(5,8,11)) for (sz in c("equal","skewed")) for (r in 1:300){ cfg<-cfg+1L; addcfg(T,G,"random",sz,41,"random",cfg) }
grid <- do.call(rbind, grid_rows); coh <- do.call(rbind, coh_rows)
write.csv(grid, file.path(OUT,"design_checkability_grid.csv"), row.names=FALSE)
write.csv(coh,  file.path(OUT,"design_checkability_cohorts_grid.csv"), row.names=FALSE)
cat(sprintf("swept configs: %d (grid %d x {N_c} + random %d); scheme-rows %d; cohort-rows %d\n",
    length(unique(grid$cfg)), sum(grid$arm=="grid" & grid$scheme=="simple"),
    sum(grid$arm=="random" & grid$scheme=="simple"), nrow(grid), nrow(coh)))
cat(sprintf("[thresholds] resolvable: MDE_jump/sigma <= %.3f (|beta|=%.2f, sigma=%.4f). Pools geometry-derived (CS: N_c; stacked: N_c + not-yet-treated).\n", THR, BETA_PP, SIGMA))

# ---- P&P anchor via the SAME grid functions ----------------------------------
yrs <- c(2001,2002,2003,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
gpp <- yrs-2000+1; ngp <- c(1,5,1,2,1,7,1,3,1,3,3,1,1); Tpp <- 21; Ncpp <- 41
kpp <- gpp-2; hpp <- (Tpp-gpp)/2
cwp <- cs_weights(gpp, ngp, Tpp); plp <- pools_of(gpp, ngp, Ncpp)
wsp <- stacked_w(gpp, ngp, Tpp, plp$stacked)
wlp <- list(simple=cwp$simple, group=cwp$group, dynamic=cwp$dynamic, calendar=cwp$calendar, stacked=wsp)
pp_rows <- list()
for (s in SCHEMES){ nc <- if (s=="stacked") plp$stacked else plp$cs
  cl <- classify(kpp, ngp, nc, wlp[[s]]); sh <- shares(cl$reason, cl$mj, wlp[[s]])
  pp_rows[[s]] <- data.frame(scheme=s, no_pre=sh["no_pre"], size_fail=sh["size"],
    underpowered=sh["underpowered"], resolvable=sh["resolvable"], wmed_jump=sh["wmed_jump"], row.names=NULL) }
pp <- do.call(rbind, pp_rows)
cat("\n== P&P anchor (grid functions; N_c=41; geometry pools) -- weight share by binding reason ==\n")
print(transform(pp, no_pre=round(no_pre,3), size_fail=round(size_fail,3), underpowered=round(underpowered,3),
  resolvable=round(resolvable,3), wmed_jump=round(wmed_jump,3)), row.names=FALSE)
cat(sprintf("stacked pool per cohort (geometry): %s\n", paste(sprintf("%d:%d",yrs,plp$stacked),collapse=" ")))

# ---- stacked 2009 under geometry vs committed pools, with/without Memphis -----
i9 <- which(yrs==2009); nc9_geo <- plp$stacked[i9]
f <- function(ng,nc) c(jump=mde_jump_sig(8,ng,nc), tbias=mde_tbias_sig(8,ng,nc,hpp[i9]))
cat("\n== stacked 2009: MDE_jump/sigma and trend-bias/sigma (k=8) ==\n")
cat(sprintf("  geometry pool  n_c=%d, Memphis IN  n_g=8: jump %.3f, tbias %.3f ; k/n_g=%.2f -> %s\n",
    nc9_geo, f(8,nc9_geo)["jump"], f(8,nc9_geo)["tbias"], 8/8, ifelse(8/8<=1,"checkable(<=1)","uncheckable")))
cat(sprintf("  geometry pool  n_c=%d, Memphis OUT n_g=7: jump %.3f, tbias %.3f ; k/n_g=%.2f -> %s\n",
    nc9_geo, f(7,nc9_geo)["jump"], f(7,nc9_geo)["tbias"], 8/7, ifelse(8/7<=1,"checkable","uncheckable(size)")))
cat(sprintf("  committed pool n_c=741, Memphis IN n_g=8: jump %.3f, tbias %.3f (reference marker)\n", f(8,741)["jump"], f(8,741)["tbias"]))
cat("  resolvable line THR =", round(THR,3), "sigma -> none of the above is resolvable\n")

# ---- gate sensitivity: k/n_g <= 1 vs < 1 (which cohorts flip) -----------------
nflip <- sum(coh$k==coh$n_g & coh$k>=1)
cat(sprintf("\n== gate sensitivity (size k/n_g<=1 vs <1) ==\n  grid cohort-cells with k/n_g exactly 1 (flip): %d of %d (%.1f%%)\n",
    nflip, nrow(coh), 100*nflip/nrow(coh)))
flipw <- tapply(coh$weight[coh$k==coh$n_g & coh$k>=1], coh$scheme[coh$k==coh$n_g & coh$k>=1], sum, na.rm=TRUE)
cat("  grid aggregation weight sitting exactly on the k/n_g=1 boundary, by scheme:\n"); print(round(flipw,3))
cat(sprintf("  P&P: stacked 2009 has k=8, n_g=8 -> k/n_g=1.00: CHECKABLE under <=1, UNCHECKABLE under <1 (flips). weight=%.3f\n",
    wlp$stacked[i9]))
cat("       all other P&P cohorts have k/n_g != 1 (no flip).\n")

# ---- binding-reason summary: grid median vs P&P ------------------------------
med <- aggregate(cbind(no_pre,size_fail,underpowered,resolvable,wmed_jump)~scheme, grid, median, na.rm=TRUE)
med <- med[match(SCHEMES,med$scheme),]
bind <- data.frame(scheme=SCHEMES,
  grid_no_pre=med$no_pre, grid_size=med$size_fail, grid_underpow=med$underpowered,
  grid_resolvable=med$resolvable, grid_wmed=med$wmed_jump,
  pp_no_pre=pp$no_pre, pp_size=pp$size_fail, pp_underpow=pp$underpowered,
  pp_resolvable=pp$resolvable, pp_wmed=pp$wmed_jump)
write.csv(bind, file.path(OUT,"design_checkability_binding.csv"), row.names=FALSE)
cat("\n== binding-reason weight shares: grid median vs P&P ==\n")
print(transform(bind, grid_no_pre=round(grid_no_pre,2), grid_size=round(grid_size,2), grid_underpow=round(grid_underpow,2),
  grid_resolvable=round(grid_resolvable,2), grid_wmed=round(grid_wmed,2),
  pp_no_pre=round(pp_no_pre,2), pp_size=round(pp_size,2), pp_underpow=round(pp_underpow,2),
  pp_resolvable=round(pp_resolvable,2), pp_wmed=round(pp_wmed,2)), row.names=FALSE)

# ---- headline table: weighted-median jump MDE/sigma, grid vs P&P -------------
summ <- bind[,c("scheme","grid_wmed","grid_resolvable","pp_wmed","pp_resolvable")]
write.csv(summ, file.path(OUT,"tab_design_checkability.csv"), row.names=FALSE)
fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "--")
L <- c("\\begin{tabular}{lrrrr}","\\toprule",
  "& \\multicolumn{2}{c}{median across designs} & \\multicolumn{2}{c}{P\\&P (grid fns)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Scheme & wtd.\\ med.\\ jump MDE/$\\sigma$ & resolvable & wtd.\\ med.\\ jump MDE/$\\sigma$ & resolvable \\\\","\\midrule",
  apply(summ,1,function(r) sprintf("%s & %s & %s & %s & %s \\\\", r[["scheme"]],
    fmt(as.numeric(r[["grid_wmed"]])), fmt(as.numeric(r[["grid_resolvable"]])),
    fmt(as.numeric(r[["pp_wmed"]])), fmt(as.numeric(r[["pp_resolvable"]])))),
  "\\bottomrule","\\end{tabular}")
writeLines(L, file.path(OUT,"tab_design_checkability.tex"))

# ---- P&P percentile (grid-function anchor only) ------------------------------
cat("\nP&P (grid-fn anchor) percentile of weighted-median jump MDE/sigma:\n")
for (s in SCHEMES){ v<-grid$wmed_jump[grid$scheme==s]; v<-v[is.finite(v)]; p<-pp$wmed_jump[pp$scheme==s]
  cat(sprintf("  %-9s wtd.med %.3f sigma -> %.0fth pct (resolvable share %.3f; grid-median resolvable %.3f)\n",
    s, p, if(is.finite(p)&&length(v)) 100*mean(v<=p) else NA, pp$resolvable[pp$scheme==s], med$resolvable[med$scheme==s])) }

# =============================================================================
# FIGURES  (headline = weighted-median jump MDE distribution; + scatter)
# =============================================================================
scheme_lv <- factor(SCHEMES, levels=SCHEMES)
g3d <- grid[is.finite(grid$wmed_jump),]; g3d$scheme <- factor(g3d$scheme, levels=SCHEMES)
ppl <- data.frame(scheme=scheme_lv, v=pp$wmed_jump)
g3 <- ggplot(g3d, aes(wmed_jump)) +
  geom_histogram(bins=30, fill="#009E73", color="white") +
  geom_vline(xintercept=THR, linetype="dashed", color="grey45") +
  geom_vline(data=ppl[is.finite(ppl$v),], aes(xintercept=v), color="#D55E00", linewidth=0.9) +
  geom_text(data=ppl[is.finite(ppl$v),], aes(x=v, y=Inf, label=sprintf("P&P %.2f",v)),
            color="#D55E00", vjust=1.4, hjust=-0.08, size=2.6) +
  scale_x_log10() + facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("weight-weighted median JUMP MDE / sigma among checkable cohorts (dashed = resolvable line %.2f)", THR),
       y="design configurations") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_checkability_wmed.pdf"), g3, width=12, height=3.6, device=cairo_pdf)
ggsave(file.path(FIG,"design_checkability_wmed.png"), g3, width=12, height=3.6, dpi=200)

# scatter: weight vs jump MDE/sigma, grid cohorts + P&P (grid-fn) overlaid
sc <- coh[coh$reason %in% c("resolvable","underpowered") & is.finite(coh$mde_jump_sig) & !is.na(coh$weight) & coh$weight>0,]
sc$scheme <- factor(sc$scheme, levels=SCHEMES)
unck <- coh[coh$reason %in% c("no_pre","size") & !is.na(coh$weight) & coh$weight>0,]; unck$scheme <- factor(unck$scheme, levels=SCHEMES)
xmax <- max(sc$mde_jump_sig, na.rm=TRUE); x_unck <- xmax*1.9; unck$mde_jump_sig <- x_unck
ppc <- do.call(rbind, lapply(SCHEMES, function(s){ nc<-if(s=="stacked") plp$stacked else plp$cs
  mj<-mde_jump_sig(kpp,ngp,nc); rs<-ifelse(is.na(kpp)|kpp<1,"no_pre",ifelse(kpp/ngp>1,"size",ifelse(mj>THR,"underpowered","resolvable")))
  data.frame(scheme=s, cohort=yrs, weight=wlp[[s]], mde_jump_sig=ifelse(rs%in%c("no_pre","size"),x_unck,mj), reason=rs) }))
ppc <- ppc[!is.na(ppc$weight) & ppc$weight>0,]; ppc$scheme <- factor(ppc$scheme, levels=SCHEMES)
g2 <- ggplot() +
  annotate("rect", xmin=x_unck/1.3, xmax=x_unck*1.3, ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.5) +
  annotate("text", x=x_unck, y=Inf, label="uncheckable", angle=90, vjust=1.2, hjust=1, size=2.6, color="grey35") +
  geom_vline(xintercept=THR, linetype="dashed", color="grey45") +
  geom_vline(xintercept=1, linetype="dotted", color="grey65") +
  geom_point(data=sc,   aes(mde_jump_sig, weight), color="#9ecae1", alpha=0.30, size=0.9) +
  geom_point(data=unck, aes(mde_jump_sig, weight), color="grey75", alpha=0.18, size=0.9, shape=4) +
  geom_point(data=ppc,  aes(mde_jump_sig, weight), color="#D55E00", size=2.0) +
  ggrepel::geom_text_repel(data=ppc, aes(mde_jump_sig, weight, label=cohort), color="#D55E00", size=2.2, max.overlaps=40, seed=1) +
  scale_x_log10() + facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("JUMP MDE / sigma   (dashed = resolvable line %.2f; dotted = one sigma; right strip = uncheckable)", THR),
       y="cohort aggregation weight") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_checkability_scatter.pdf"), g2, width=13, height=3.8, device=cairo_pdf)
ggsave(file.path(FIG,"design_checkability_scatter.png"), g2, width=13, height=3.8, dpi=200)
for (f in c("design_checkability_wmed","design_checkability_scatter"))
  file.copy(file.path(FIG,paste0(f,".png")), file.path(OUT,paste0(f,".png")), overwrite=TRUE)

cat("\nWrote design_checkability_grid.csv, _cohorts_grid.csv, _binding.csv, tab_design_checkability.{csv,tex},\n")
cat("  figures/design_checkability_{wmed,scatter}.*\n")
