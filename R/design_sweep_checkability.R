# =============================================================================
# design_sweep_checkability.R  --  weight-vs-checkability (analytic MDE) across
# the design space. The sweep version of weight_checkability.R / Figure C.
#
# DETERMINISTIC. No outcome variable, no DGP, no sampling of data, no
# replications. For every design configuration (panel length T, adoption-period
# vector g, cohort sizes n_g) and every aggregation scheme we compute, per cohort:
#   - its aggregation weight  (exact, from design_sweep_core.R);
#   - whether its joint pre-trend test is CHECKABLE: k>=1 pre-leads AND k/n_g<=1
#     (the B1/B2 rank/size gate); and
#   - if checkable, the analytic MINIMUM DETECTABLE pre-trend (MDE) at 80% power,
#     alpha 0.05, EXPRESSED IN UNITS OF sigma (the outcome's noise SD).
#
# Why sigma units: MDE = sqrt(lambda*(k)/denom) * sigma * sqrt(1/n_g + 1/n_c), so
# MDE/sigma = sqrt(lambda*(k)/denom) * sqrt(1/n_g + 1/n_c) is a PURE function of
# the geometry (k, n_g, n_c) -- sigma cancels; the sweep needs no outcomes. To
# connect to weight_checkability.R's effect-relative axis we mark the line
# MDE = a P&P-scale pre-trend |beta|=0.10 at the committed sigma=0.3297, i.e.
# MDE/sigma = 0.10/0.3297 = 0.303 ("resolvable" at/below it).
#
# CONTROL POOLS are scheme-specific and realistic, as in the data: the CS family
# (simple/group/dynamic/calendar) and SA use the NEVER-TREATED pool (~41); the
# stacked estimator uses the much larger NOT-YET-TREATED pool (~741). The pool is
# the first-order driver of the resolvable line (1/n_c floor), so it is carried,
# not swept away. The sweep varies the geometry (T, placement, sizes); within a
# config all schemes share one geometry and differ only in weights + pool.
#
# The P&P ANCHOR is taken EXACTLY from the committed weight_checkability_cohorts
# .csv (real per-scheme n_g, cohort membership, pools -- including the Memphis
# 2009 n_g=8-vs-7 and the 2002 no-clean-controls stack), not re-derived; the grid
# shows the surrounding geometric landscape. Balanced grid abstractions from the
# real data (equal cell counts, a single n_g per cohort across schemes) are noted.
# =============================================================================
suppressMessages({library(ggplot2); library(tidyr)})
.sourced_quiet <- TRUE
source("R/design_sweep_core.R")
OUT <- Sys.getenv("DIDREP_OUT", unset="output"); FIG <- Sys.getenv("DIDREP_FIG", unset="figures")

NC_CS   <- 41        # never-treated pool  (CS family + SA)
NC_ST   <- 741       # not-yet-treated pool (stacked)
SIGMA   <- 0.3297    # committed never-treated two-way-FE residual SD (weight_checkability.R)
BETA_PP <- 0.10      # a P&P-scale fatal-encounters pre-trend, for the "resolvable" line
THR_EFF <- BETA_PP/SIGMA   # MDE/sigma below which a 0.10 pre-trend is detectable (0.303)
THR_SIG <- 1.0             # secondary: detect a pre-trend as large as one sigma
SCHEMES <- c("simple","group","dynamic","calendar","stacked")
POOL    <- c(simple=NC_CS, group=NC_CS, dynamic=NC_CS, calendar=NC_CS, stacked=NC_ST)

# ---- analytic MDE/sigma: lambda*(k) and denom precomputed per k --------------
lambda_star <- function(k){ if (k<1) return(NA_real_); crit <- qchisq(0.95, k)
  uniroot(function(l) pchisq(crit, df=k, ncp=l, lower.tail=FALSE)-0.80, c(1e-6, 400), tol=1e-10)$root }
denom_shape <- function(k, shape){ a <- if (shape=="trend") (1:k) else rep(1,k)
  M <- solve(diag(k) + matrix(1,k,k)); as.numeric(t(a) %*% M %*% a) }     # a'(I+11')^-1 a
KMAX <- 40; kseq <- 1:KMAX
Ctrend <- sqrt(sapply(kseq,lambda_star)/sapply(kseq,denom_shape,shape="trend"))
mde_sig <- function(k, n_g, n_c){                                        # trend shape, sigma units
  ok <- !is.na(k) & k>=1 & k<=KMAX & n_c>0 & n_g>0
  out <- rep(NA_real_, length(k)); out[ok] <- Ctrend[k[ok]]*sqrt(1/n_g[ok] + 1/n_c)
  out }

# ---- per-scheme weighted checkability shares for one design ------------------
scheme_shares <- function(k, n_g, w, n_c){
  ident <- !is.na(w); tot <- sum(w[ident])
  if (tot<=0) return(c(uncheckable=NA,underpow_eff=NA,good_eff=NA,good_sig=NA,wmed_mde=NA,n_ident=sum(ident)))
  gate_ok <- ident & !is.na(k) & k>=1 & (k/n_g <= 1)
  ms <- mde_sig(k, n_g, n_c)
  unck   <- ident & !gate_ok
  good_e <- gate_ok & is.finite(ms) & ms <= THR_EFF
  good_s <- gate_ok & is.finite(ms) & ms <= THR_SIG
  under_e<- gate_ok & !good_e
  okm <- gate_ok & is.finite(ms)
  wmed <- if (any(okm)){ o<-order(ms[okm]); cw<-cumsum(w[okm][o])/sum(w[okm]); ms[okm][o][which(cw>=0.5)[1]] } else NA_real_
  c(uncheckable=sum(w[unck])/tot, underpow_eff=sum(w[under_e])/tot, good_eff=sum(w[good_e])/tot,
    good_sig=sum(w[good_s])/tot, wmed_mde=wmed, n_ident=sum(ident)) }

# ---- design generators (identical universe to design_sweep_figureC.R) --------
place <- function(T, G, type){ slots <- 2:T
  if (type=="early") head(slots,G) else if (type=="late") tail(slots,G)
  else if (type=="spread") unique(round(seq(2,T,length.out=G))) else sort(sample(slots,G)) }
sizes <- function(G, type){ if (type=="equal") rep(3L,G)
  else { s<-rep(1L,G); s[1]<-8L; if(G>=2) s[2]<-5L; if(G>=3) s[3]<-3L; s } }

grid_rows <- list(); coh_rows <- list()
addcfg <- function(T,G,pl,sz,arm,cfgid){
  gp <- place(T,G,pl); if (length(gp)<G) return(invisible())
  ng <- sizes(length(gp), sz)
  dw <- design_weights(T, gp, ng, NC_ST)       # CS weights pool-independent; stacked w_s at NC_ST
  wcols <- list(simple=dw$w_simple, group=dw$w_group, dynamic=dw$w_dynamic,
                calendar=dw$w_calendar, stacked=dw$w_stacked)
  for (s in SCHEMES){ sh <- scheme_shares(dw$k, dw$n_g, wcols[[s]], POOL[[s]])
    grid_rows[[length(grid_rows)+1L]] <<- data.frame(cfg=cfgid, T=T, G=length(gp),
      placement=pl, size=sz, arm=arm, scheme=s,
      uncheckable=sh["uncheckable"], underpow_eff=sh["underpow_eff"], good_eff=sh["good_eff"],
      good_sig=sh["good_sig"], wmed_mde=sh["wmed_mde"], n_ident=sh["n_ident"], row.names=NULL) }
  if (arm=="grid"){
    for (s in SCHEMES){ ms <- mde_sig(dw$k, dw$n_g, POOL[[s]])
      gate_ok <- !is.na(dw$k) & dw$k>=1 & (dw$k/dw$n_g <= 1)
      coh_rows[[length(coh_rows)+1L]] <<- data.frame(cfg=cfgid, T=T, placement=pl, size=sz,
        scheme=s, cohort=dw$cohort, gp=dw$gp, k=dw$k, n_g=dw$n_g, n_c=POOL[[s]],
        weight=wcols[[s]], gate=ifelse(gate_ok,"ok",ifelse(is.na(dw$k)|dw$k<1,"no_pre","size")),
        mde_sigma=ifelse(gate_ok, ms, NA_real_), row.names=NULL) } }
}
cfg <- 0L
for (T in c(15,21,30)) for (G in c(3,5,8,11)) for (pl in c("early","late","spread"))
  for (sz in c("equal","skewed")) if (G <= T-1){ cfg<-cfg+1L; addcfg(T,G,pl,sz,"grid",cfg) }
set.seed(11)
for (T in c(21)) for (G in c(5,8,11)) for (sz in c("equal","skewed")) for (r in 1:300){
  cfg<-cfg+1L; addcfg(T,G,"random",sz,"random",cfg) }
grid <- do.call(rbind, grid_rows); coh <- do.call(rbind, coh_rows)
write.csv(grid, file.path(OUT,"design_checkability_grid.csv"), row.names=FALSE)
write.csv(coh,  file.path(OUT,"design_checkability_cohorts_grid.csv"), row.names=FALSE)
cat(sprintf("swept configurations: %d (grid %d + random %d); scheme-rows %d; grid cohort-rows %d\n",
    length(unique(grid$cfg)), sum(grid$arm=="grid")/length(SCHEMES),
    sum(grid$arm=="random")/length(SCHEMES), nrow(grid), nrow(coh)))
cat(sprintf("[pools] CS family + SA: n_c=%d ; stacked: n_c=%d\n", NC_CS, NC_ST))
cat(sprintf("[thresholds] resolvable (good_eff): MDE/sigma <= %.3f  (|beta|=%.2f at sigma=%.4f); good_sig: MDE/sigma <= 1\n",
    THR_EFF, BETA_PP, SIGMA))

# ---- P&P anchor: EXACT, from committed weight_checkability_cohorts.csv --------
wc <- read.csv(file.path(OUT,"weight_checkability_cohorts.csv"))
wc$mde_sigma <- wc$mde_trend/SIGMA
# dynamic/calendar not in that file -> same geometry/pool as CS, weights from task2
w_cs <- read.csv(file.path(OUT,"task2_cohort_weights.csv"))
cs_geo <- unique(wc[wc$scheme=="CS simple", c("cohort","n_g","k","gate","mde_sigma")])
anchor_scheme <- function(label){
  if (label %in% c("dynamic","calendar")){
    d <- merge(cs_geo, w_cs[,c("cohort",label)], by="cohort"); w <- d[[label]]
    k<-d$k; n_g<-d$n_g; gate<-d$gate; ms<-d$mde_sigma
  } else {
    key <- c(simple="CS simple", group="CS group", stacked="stacked")[label]
    d <- wc[wc$scheme==key,]; w<-d$weight; k<-d$k; n_g<-d$n_g; gate<-d$gate; ms<-d$mde_sigma }
  gate_ok <- gate=="ok"; tot<-sum(w)
  good_e <- gate_ok & is.finite(ms) & ms<=THR_EFF
  okm <- gate_ok & is.finite(ms)
  wmed <- if (any(okm)){ o<-order(ms[okm]); cw<-cumsum(w[okm][o])/sum(w[okm]); ms[okm][o][which(cw>=0.5)[1]] } else NA_real_
  data.frame(scheme=label, uncheckable=sum(w[!gate_ok])/tot, good_eff=sum(w[good_e])/tot,
             wmed_mde=wmed, row.names=NULL) }
pp_df <- do.call(rbind, lapply(SCHEMES, anchor_scheme))
cat("\n== P&P anchor (exact, from weight_checkability_cohorts.csv; real per-scheme pools) ==\n")
print(transform(pp_df, uncheckable=round(uncheckable,3), good_eff=round(good_eff,3), wmed_mde=round(wmed_mde,3)), row.names=FALSE)

# ---- validate the sweep's MDE/sigma reproduces the committed per-cohort MDE --
cat("\nvalidation vs weight_checkability_cohorts.csv:\n")
for (rc in list(c("CS simple",2002,NC_CS), c("stacked",2009,NC_ST))){
  r <- wc[wc$scheme==rc[[1]] & wc$cohort==as.integer(rc[[2]]),][1,]
  got <- mde_sig(r$k, r$n_g, as.numeric(rc[[3]])); want <- r$mde_sigma
  cat(sprintf("  %-9s %d: sweep %.4f vs committed %.4f (diff %.2e)\n", rc[[1]], as.integer(rc[[2]]), got, want, abs(got-want))) }

# ---- P&P percentile in the swept distribution (per scheme) -------------------
pctile <- function(s, col, p){ v<-grid[[col]][grid$scheme==s]; v<-v[is.finite(v)]
  if (!is.finite(p)||!length(v)) NA_real_ else mean(v<=p) }
cat("\nP&P vs swept distribution (per scheme):\n")
for (s in SCHEMES){ p<-pp_df[pp_df$scheme==s,]
  cat(sprintf("  %-9s resolvable %.3f (%.0fth pct) | uncheckable %.3f (%.0fth pct) | wtd.med MDE/sig %.3f (%.0fth pct)\n",
    s, p$good_eff, 100*pctile(s,"good_eff",p$good_eff), p$uncheckable, 100*pctile(s,"uncheckable",p$uncheckable),
    p$wmed_mde, 100*pctile(s,"wmed_mde",p$wmed_mde))) }

# ---- summary table (median across configs + P&P), csv + tex -------------------
med <- aggregate(cbind(good_eff,uncheckable,wmed_mde)~scheme, grid, median, na.rm=TRUE)
med <- med[match(SCHEMES, med$scheme),]
summ <- data.frame(scheme=SCHEMES, grid_good_eff=med$good_eff, grid_uncheckable=med$uncheckable,
  grid_wmed_mde=med$wmed_mde, pp_good_eff=pp_df$good_eff, pp_uncheckable=pp_df$uncheckable, pp_wmed_mde=pp_df$wmed_mde)
write.csv(summ, file.path(OUT,"tab_design_checkability.csv"), row.names=FALSE)
cat("\n== summary: median across designs vs P&P ==\n")
print(transform(summ, grid_good_eff=round(grid_good_eff,3), grid_uncheckable=round(grid_uncheckable,3),
  grid_wmed_mde=round(grid_wmed_mde,3), pp_good_eff=round(pp_good_eff,3),
  pp_uncheckable=round(pp_uncheckable,3), pp_wmed_mde=round(pp_wmed_mde,3)), row.names=FALSE)
fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "--")
L <- c("\\begin{tabular}{lrrrrrr}","\\toprule",
  "& \\multicolumn{3}{c}{median across designs} & \\multicolumn{3}{c}{P\\&P (actual)} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Scheme & resolvable & uncheckable & wtd.\\ med.\\ MDE/$\\sigma$ & resolvable & uncheckable & wtd.\\ med.\\ MDE/$\\sigma$ \\\\",
  "\\midrule",
  apply(summ,1,function(r) sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", r[["scheme"]],
    fmt(as.numeric(r[["grid_good_eff"]])), fmt(as.numeric(r[["grid_uncheckable"]])), fmt(as.numeric(r[["grid_wmed_mde"]])),
    fmt(as.numeric(r[["pp_good_eff"]])), fmt(as.numeric(r[["pp_uncheckable"]])), fmt(as.numeric(r[["pp_wmed_mde"]])))),
  "\\bottomrule","\\end{tabular}")
writeLines(L, file.path(OUT,"tab_design_checkability.tex"))

# =============================================================================
# FIGURES
# =============================================================================
scheme_lv <- factor(SCHEMES, levels=SCHEMES)
# ---- Figure 1: distribution of the resolvable-weight share, P&P marked -------
g1d <- grid[is.finite(grid$good_eff),]; g1d$scheme <- factor(g1d$scheme, levels=SCHEMES)
ppl <- data.frame(scheme=scheme_lv, share=pp_df$good_eff)
g1 <- ggplot(g1d, aes(good_eff)) +
  geom_histogram(breaks=seq(0,1,0.05), fill="#0072B2", color="white") +
  geom_vline(data=ppl, aes(xintercept=share), color="#D55E00", linewidth=0.9) +
  geom_text(data=ppl, aes(x=share, y=Inf, label=sprintf("P&P %.2f", share)),
            color="#D55E00", vjust=1.4, hjust=-0.08, size=2.7) +
  facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("share of aggregation weight on CHECKABLE & RESOLVABLE cohorts (gate ok, MDE <= %.2f sigma)", THR_EFF),
       y="design configurations") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_checkability_share.pdf"), g1, width=12, height=3.6, device=cairo_pdf)
ggsave(file.path(FIG,"design_checkability_share.png"), g1, width=12, height=3.6, dpi=200)

# ---- Figure 2: the comparison itself, swept -- weight vs MDE/sigma -----------
sc <- coh[coh$gate=="ok" & is.finite(coh$mde_sigma) & !is.na(coh$weight) & coh$weight>0,]
sc$scheme <- factor(sc$scheme, levels=SCHEMES)
unck <- coh[coh$gate!="ok" & !is.na(coh$weight) & coh$weight>0,]; unck$scheme <- factor(unck$scheme, levels=SCHEMES)
xmax <- max(sc$mde_sigma, na.rm=TRUE); x_unck <- xmax*1.9; unck$mde_sigma <- x_unck
# P&P cohorts, real pools, from committed table
ppc <- rbind(
  transform(wc[wc$scheme %in% c("CS simple","CS group","stacked"),
    c("scheme","cohort","weight","gate","mde_sigma")],
    scheme=c("CS simple"="simple","CS group"="group","stacked"="stacked")[scheme]),
  { d<-merge(cs_geo, w_cs[,c("cohort","dynamic","calendar")], by="cohort")
    rbind(data.frame(scheme="dynamic",  cohort=d$cohort, weight=d$dynamic,  gate=d$gate, mde_sigma=d$mde_sigma),
          data.frame(scheme="calendar", cohort=d$cohort, weight=d$calendar, gate=d$gate, mde_sigma=d$mde_sigma)) })
ppc <- ppc[!is.na(ppc$weight) & ppc$weight>0,]
ppc$mde_sigma[ppc$gate!="ok"] <- x_unck; ppc$scheme <- factor(ppc$scheme, levels=SCHEMES)
g2 <- ggplot() +
  annotate("rect", xmin=x_unck/1.3, xmax=x_unck*1.3, ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.5) +
  annotate("text", x=x_unck, y=Inf, label="uncheckable", angle=90, vjust=1.2, hjust=1, size=2.6, color="grey35") +
  geom_vline(xintercept=THR_EFF, linetype="dashed", color="grey45") +
  geom_vline(xintercept=1, linetype="dotted", color="grey65") +
  geom_point(data=sc,   aes(mde_sigma, weight), color="#9ecae1", alpha=0.30, size=0.9) +
  geom_point(data=unck, aes(mde_sigma, weight), color="grey75", alpha=0.20, size=0.9, shape=4) +
  geom_point(data=ppc,  aes(mde_sigma, weight), color="#D55E00", size=2.0) +
  ggrepel::geom_text_repel(data=ppc, aes(mde_sigma, weight, label=cohort),
       color="#D55E00", size=2.2, max.overlaps=40, seed=1) +
  scale_x_log10() + facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("MDE / sigma   (dashed = |beta|=%.2f pre-trend; dotted = one sigma; right strip = uncheckable)", BETA_PP),
       y="cohort aggregation weight") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_checkability_scatter.pdf"), g2, width=13, height=3.8, device=cairo_pdf)
ggsave(file.path(FIG,"design_checkability_scatter.png"), g2, width=13, height=3.8, dpi=200)

# ---- Figure 3: distribution of weighted-median MDE/sigma, P&P marked ---------
g3d <- grid[is.finite(grid$wmed_mde),]; g3d$scheme <- factor(g3d$scheme, levels=SCHEMES)
ppl3 <- data.frame(scheme=scheme_lv, v=pp_df$wmed_mde)
g3 <- ggplot(g3d, aes(wmed_mde)) +
  geom_histogram(bins=30, fill="#009E73", color="white") +
  geom_vline(xintercept=THR_EFF, linetype="dashed", color="grey45") +
  geom_vline(data=ppl3[is.finite(ppl3$v),], aes(xintercept=v), color="#D55E00", linewidth=0.9) +
  scale_x_log10() + facet_wrap(~scheme, ncol=5) +
  labs(x=sprintf("weight-weighted median MDE / sigma among checkable cohorts (dashed = %.2f)", THR_EFF),
       y="design configurations") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"design_checkability_wmed.pdf"), g3, width=12, height=3.6, device=cairo_pdf)
ggsave(file.path(FIG,"design_checkability_wmed.png"), g3, width=12, height=3.6, dpi=200)
for (f in c("design_checkability_share","design_checkability_scatter","design_checkability_wmed"))
  file.copy(file.path(FIG,paste0(f,".png")), file.path(OUT,paste0(f,".png")), overwrite=TRUE)

cat("\nWrote design_checkability_grid.csv, _cohorts_grid.csv, tab_design_checkability.{csv,tex},\n")
cat("  figures/design_checkability_{share,scatter,wmed}.*\n")
