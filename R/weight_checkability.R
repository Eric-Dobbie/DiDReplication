# =============================================================================
# weight_checkability.R  --  weight vs. checkability (analytic MDE) diagnostic
#
# Per cohort g, per scheme {CS group, CS simple, SA, stacked}: the cohort's
# aggregation weight vs an analytic minimum detectable pre-trend violation (MDE),
# with binary checkability gates. R5 settings reused (never-treated controls, no
# covariates, universal base period, no binning); cohorts/sample unchanged.
# NB: estimator_decomposition_instructions.md is not present on disk this session;
# the R5 conventions are taken from run_extraction.R and the A-chain (as the task
# restates them). Stacked uses the corrected unit x stack + year x cohort FE (R2).
#
# Weights reused where available: CS simple/group from task2_cohort_weights.csv;
# stacked FWL w_s from fwl_decomp_unweighted.csv (sum to 1). SA computed here.
#
# UNITS FIX (supersedes the first version): the treatment effect is a LEVEL, so a
# level-on-level comparison is used. The "trend" MDE from mde_shape(a=1:k) is a
# per-period SLOPE; comparing it to the level effect understated the detectable
# violation. We now report (i) the JUMP MDE -- the detectable uniform level shift
# (a=rep(1,k)) -- as the level-comparable measure, and (ii) the slope MDE
# CONVERTED to an implied mean post-period bias = slope * h, h = mean post event
# time = (2020 - g)/2. Both are compared to |aggregate effect|. See stacked 2009.
# =============================================================================
suppressMessages({library(did); library(fixest); library(lfe); library(ggplot2); library(ggrepel)})
DATA <- Sys.getenv("DIDREP_DATA", unset="data"); OUT <- Sys.getenv("DIDREP_OUT", unset="output")
FIG  <- Sys.getenv("DIDREP_FIG",  unset="figures")
FIRST_OUT <- 2000

dta <- read.csv(file.path(DATA,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$yc <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta$gname <- ifelse(is.na(dta$yc)&dta$no.req==0,0, ifelse(is.na(dta$yc)&dta$no.req==1,1987,dta$yc))
never_ids <- unique(dta$agency.id[dta$gname==0])
st <- read.csv(file.path(DATA,"stacked_fatal.csv")); st$agency_stack <- paste(st$agency.id, st$cohort, sep="__")

# ---- sigma: residual SD, unit+year FE on untreated (never-treated) obs -------
uu <- dta[dta$gname==0 & !is.na(dta$any.fatalities), ]
sig <- sd(felm(any.fatalities ~ 1 | agency.num + year, data=uu)$residuals[,1])
pbar <- mean(uu$any.fatalities); sig_bin <- sqrt(pbar*(1-pbar))
cat(sprintf("sigma (unit+year FE residual SD, untreated) = %.4f ; sqrt(p(1-p)) check = %.4f (pbar=%.3f)\n", sig, sig_bin, pbar))

# ---- CS fit (universal base) for the 2002 sanity check -----------------------
set.seed(0)
mp <- did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="gname",
                  xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period="universal",data=dta)
nU <- mp$DIDparams$n; Vcs <- as.matrix(mp$V_analytical)/nU
i2002 <- which(mp$group==2002 & mp$t==2000)                 # 2002's single placebo lead (e=-2)
se_attgt_2002 <- if(length(i2002)) sqrt(Vcs[i2002,i2002]) else NA_real_

# ---- weights -----------------------------------------------------------------
w_cs <- read.csv(file.path(OUT,"task2_cohort_weights.csv"))     # simple, group, ...
fw   <- read.csv(file.path(OUT,"fwl_decomp_unweighted.csv")); fw <- fw[fw$w_s>0 & !is.na(fw$beta_s),]
# SA cohort weights (att aggregation obs-count shares), same sample/recode as extract_sa_atoms
d_sa <- dta[!is.na(dta$any.fatalities) & (dta$gname==0 | dta$gname>=2001), ]; d_sa$coh <- ifelse(d_sa$gname==0,10000,d_sa$gname)
res <- fixest::feols(any.fatalities ~ sunab(coh, year) | agency.num + year, data=d_sa, cluster=~agency.num)
ct <- summary(res, agg=FALSE)$coeftable; mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee <- as.integer(vapply(mm,function(z)z[2],character(1))); gg <- as.integer(vapply(mm,function(z)z[3],character(1)))
mmx <- stats::model.matrix(res); sh <- colSums(abs(sign(mmx)))[rownames(ct)]
wpost <- ifelse(!is.na(ee)&ee>=0, sh, 0); wpost[is.na(wpost)] <- 0
wSA <- tapply(wpost, gg, sum); wSA <- wSA/sum(wpost)
att_SA <- as.numeric(fixest::coeftable(summary(res, agg="att"))["ATT",1])
agg <- read.csv(file.path(OUT,"task2_aggregation_sensitivity.csv"))
AGG <- c(group=agg$overall[agg$scheme=="group"], simple=agg$overall[agg$scheme=="simple"],
         SA=att_SA, stacked=-0.098915)

# ---- lambda* (power 0.80, alpha 0.05) and MDE for a shape a -------------------
lambda_star <- function(k){ if(k<1) return(NA_real_); crit<-qchisq(0.95,k)
  uniroot(function(l) pchisq(crit, df=k, ncp=l, lower.tail=FALSE)-0.80, c(1e-6, 200))$root }
mde_shape <- function(k, n_g, n_c, sigma, shape){
  if (is.na(k)||k<1||n_c==0) return(NA_real_)
  a <- if (shape=="trend") (1:k) else rep(1,k)
  M <- solve(diag(k) + matrix(1,k,k)); denom <- as.numeric(t(a)%*%M%*%a)
  sqrt(lambda_star(k)/denom) * sigma * sqrt(1/n_g + 1/n_c)
}

# ---- cohort universe + per-scheme rows ---------------------------------------
cs_cohorts <- sort(w_cs$cohort)                                    # 13
n_g_cs <- sapply(cs_cohorts, function(g) length(unique(dta$agency.id[which(dta$yc==g)])))  # which() avoids NA-index
n_g_st <- sapply(fw$stack,   function(g) length(unique(st$agency.id[st$cohort==g & st$treat==1])))
n_c_st <- sapply(fw$stack,   function(g) length(unique(st$agency.id[st$cohort==g & st$treat==0])))
k_of  <- function(g) (g - FIRST_OUT) - 1                           # pre yrs observed - 1 (base)

mkrow <- function(g, scheme, weight, n_g, n_c){
  aggkey <- c("CS group"="group","CS simple"="simple","SA"="SA","stacked"="stacked")[scheme]
  eff <- abs(unname(AGG[aggkey]))
  k <- k_of(g); h <- (2020 - g)/2                       # mean post event time (slope -> level)
  no_pre <- k < 1; no_controls <- n_c==0; size <- (!is.na(k) & n_g>0 & (k/n_g > 1))
  gate_fail <- no_pre | no_controls | size
  gate <- if (no_pre) "no_pre" else if (no_controls) "no_controls" else if (size) "size" else "ok"
  slope <- if (gate_fail) NA_real_ else mde_shape(k,n_g,n_c,sig,"trend")   # per-period SLOPE MDE
  mj    <- if (gate_fail) NA_real_ else mde_shape(k,n_g,n_c,sig,"jump")    # uniform LEVEL (jump) MDE
  tbias <- slope * h                                                       # implied mean post bias (LEVEL)
  data.frame(scheme=scheme, cohort=g, weight=weight, n_g=n_g, n_c=n_c, k=k, k_over_n=k/n_g,
    no_pre=no_pre, no_controls=no_controls, size=size, gate=gate, h_post=h,
    mde_trend_slope=slope, mde_trend_bias=tbias, mde_jump=mj,
    mde_rel_tbias=tbias/eff, mde_rel_jump=mj/eff, stringsAsFactors=FALSE, row.names=NULL)
}
rows <- list()
for (i in seq_along(cs_cohorts)){ g<-cs_cohorts[i]
  rows[[length(rows)+1]] <- mkrow(g,"CS group",  w_cs$group [w_cs$cohort==g], n_g_cs[i], 41)
  rows[[length(rows)+1]] <- mkrow(g,"CS simple", w_cs$simple[w_cs$cohort==g], n_g_cs[i], 41)
  rows[[length(rows)+1]] <- mkrow(g,"SA",        as.numeric(wSA[as.character(g)]), n_g_cs[i], 41)
}
for (i in seq_along(fw$stack)){ g<-fw$stack[i]
  rows[[length(rows)+1]] <- mkrow(g,"stacked", fw$w_s[fw$stack==g], n_g_st[i], n_c_st[i]) }
D <- do.call(rbind, rows); D$weight[is.na(D$weight)] <- 0
write.csv(D, file.path(OUT,"weight_checkability_cohorts.csv"), row.names=FALSE)

# ---- sanity check: CS 2002 analytic lead SE vs att_gt ------------------------
se_analytic_2002 <- sig*sqrt((1/n_g_cs[cs_cohorts==2002] + 1/41)*2)   # k=1: Var=2*sig^2*(1/n_g+1/n_c)
cat(sprintf("\nSANITY (CS 2002 single lead): analytic SE %.4f vs att_gt SE %.4f ; ratio att_gt/analytic = %.2f\n",
    se_analytic_2002, se_attgt_2002, se_attgt_2002/se_analytic_2002))

# ---- summary table (level-on-level: JUMP MDE vs |effect|) ---------------------
summ <- do.call(rbind, lapply(unique(D$scheme), function(s){ d<-D[D$scheme==s,]
  unck <- d$gate!="ok"; mj1 <- !unck & !is.na(d$mde_rel_jump) & d$mde_rel_jump>1
  res  <- !unck & !is.na(d$mde_rel_jump) & d$mde_rel_jump<=1
  data.frame(scheme=s, w_uncheckable=sum(d$weight[unck]), w_jump_gt_eff=sum(d$weight[mj1]),
             w_resolvable=sum(d$weight[res]), w_hardunion=sum(d$weight[unck|mj1])) }))
options(width=170); cat("\n== summary: weight on hard-to-check cohorts (LEVEL / jump MDE) ==\n"); print(transform(summ,
  w_uncheckable=round(w_uncheckable,3), w_jump_gt_eff=round(w_jump_gt_eff,3),
  w_resolvable=round(w_resolvable,3), w_hardunion=round(w_hardunion,3)), row.names=FALSE)
write.csv(summ, file.path(OUT,"tab_weight_checkability.csv"), row.names=FALSE)
L <- c("\\begin{tabular}{lrrrr}","\\toprule",
  "Scheme & Uncheckable & jump MDE$>|$eff$|$ & Resolvable & Hard (union) \\\\","\\midrule",
  apply(summ,1,function(r) sprintf("%s & %.3f & %.3f & %.3f & %.3f \\\\", r[1], as.numeric(r[2]),
        as.numeric(r[3]), as.numeric(r[4]), as.numeric(r[5]))),
  "\\bottomrule","\\end{tabular}")
writeLines(L, file.path(OUT,"tab_weight_checkability.tex"))

# ---- plot --------------------------------------------------------------------
make_plot <- function(relcol, file, xlab){
  D$mrel <- D[[relcol]]
  fin <- D[is.finite(D$mrel) & D$gate=="ok", ]
  xmax <- max(fin$mrel, na.rm=TRUE); x_unck <- xmax*1.8; strip_lo <- xmax*1.3; strip_hi <- xmax*2.4
  unc <- D[D$gate!="ok", ]; unc$mrel <- x_unck
  pal <- c(no_pre="#D55E00", no_controls="#CC79A7", size="#E69F00")
  g <- ggplot() +
    annotate("rect", xmin=strip_lo, xmax=strip_hi, ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.5) +
    annotate("text", x=sqrt(strip_lo*strip_hi), y=Inf, label="Uncheckable", angle=90, vjust=1.2, hjust=1, size=2.7, color="grey35") +
    geom_vline(xintercept=1, linetype="dashed", color="grey50") +
    geom_point(data=fin, aes(mrel, weight), color="#0072B2", size=2.4) +
    geom_point(data=unc, aes(mrel, weight, shape=gate), color="black", size=2.4) +
    ggrepel::geom_text_repel(data=rbind(fin[,c("scheme","cohort","weight","mrel")], unc[,c("scheme","cohort","weight","mrel")]),
       aes(mrel, weight, label=cohort), size=2.4, max.overlaps=30, seed=1) +
    scale_shape_manual(values=c(no_pre=4, no_controls=17, size=15), name="failed gate") +
    scale_x_log10() + facet_wrap(~scheme, scales="free_y", ncol=2) +
    labs(x=xlab, y="cohort aggregation weight") +
    theme_minimal(base_size=11) + theme(legend.position="top", panel.grid.minor=element_blank())
  ggsave(file.path(FIG, paste0(file,".pdf")), g, width=10, height=7, device=cairo_pdf)
  ggsave(file.path(FIG, paste0(file,".png")), g, width=10, height=7, dpi=200)
}
# primary: JUMP MDE (level, a level-on-level comparison with the effect)
make_plot("mde_rel_jump","fig_weight_checkability",
  "relative jump MDE  (detectable uniform level pre-trend / |aggregate effect|; dashed = MDE equals effect)")
# secondary: slope MDE converted to implied mean post-period bias (also a level)
make_plot("mde_rel_tbias","fig_weight_checkability_bias",
  "relative trend-bias MDE  (detectable slope x mean post horizon / |aggregate effect|; dashed = equals effect)")
for (f in c("fig_weight_checkability","fig_weight_checkability_bias"))
  for (e in c(".pdf",".png")) file.copy(file.path(FIG,paste0(f,e)), file.path(OUT,paste0(f,e)), overwrite=TRUE)

# ---- stacked 2009 under both level measures (the headline correction) --------
r09 <- D[D$scheme=="stacked" & D$cohort==2009,][1,]
cat(sprintf("\nstacked 2009 (n_g=%d, n_c=%d, k=%d, h=%.1f): jump MDE %.4f (%.3f sigma, %.2f x |eff|) ; trend-bias MDE %.4f (%.3f sigma, %.2f x |eff|)\n",
  r09$n_g, r09$n_c, r09$k, r09$h_post, r09$mde_jump, r09$mde_jump/sig, r09$mde_rel_jump,
  r09$mde_trend_bias, r09$mde_trend_bias/sig, r09$mde_rel_tbias))
cat(sprintf("  (slope MDE was %.4f = %.3f sigma; comparing THAT to the level effect was the units bug)\n",
  r09$mde_trend_slope, r09$mde_trend_slope/sig))

cat("\n== per-cohort table (level measures) ==\n"); print(transform(D[,c("scheme","cohort","weight","n_g","n_c","k","k_over_n","gate","mde_jump","mde_rel_jump","mde_rel_tbias")],
  weight=round(weight,3), k_over_n=round(k_over_n,2), mde_jump=round(mde_jump,3),
  mde_rel_jump=round(mde_rel_jump,2), mde_rel_tbias=round(mde_rel_tbias,2)), row.names=FALSE)
cat("\nWrote weight_checkability_cohorts.csv, tab_weight_checkability.{csv,tex}, fig_weight_checkability{,_bias}.*\n")
