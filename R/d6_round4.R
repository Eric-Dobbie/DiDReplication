# =============================================================================
# d6_round4.R  --  D6 round-4 follow-ups.
#
# 1. Exact CS-imputation gap at event time 0 (s=1):
#    Delta^Imp_{g,0} := tau_imp(g,g) - tau_CS-NYT(g,g), per cohort, verified to
#    machine precision (2001 should be 0). Decompose each cohort's AVERAGE gap
#    into (a) own last-pre-period block bias and (b) inherited adjustment-cohort
#    term (exact remainder). The Liu (2025) Props 2.1/2.2 weights
#    N_k / (sum_{j>=k} N_j + N_inf) are constructed from cohort sizes and
#    reported as the theoretical attribution of (b) across adjustment cohorts.
#    NOTE: Liu (2025) is not on hand; the TOTAL decomposition is exact by
#    construction, and the per-cohort weight attribution is a reconstruction of
#    the stated formula -- flagged as such. This becomes the MAIN gap figure;
#    the round-3 linear-slope predicted-vs-actual plot is moved to the appendix.
#
# 2. mde_g for ALL cohorts via slope-only restricted wild cluster bootstrap
#    (Webb, B=999): mde = Zc * SE_WCR, SE_WCR = SD of restricted bootstrap
#    slopes. Add n_treated_clusters. Redraw the main figure on the new x-axis
#    and report which cohorts change rank vs the CR1 mde.
#
# 3. RI diagnostic: pre-period outcome mean/variance and agency size (officers,
#    population) of each treated cohort vs. the never-treated placebo pool; and
#    joint RI with placebo cohorts drawn from SIZE-MATCHED never-treated
#    agencies, reported alongside WCR. (Presented, not adjudicated.)
# =============================================================================
suppressMessages({library(fixest); library(did); library(ggplot2); library(ggrepel); library(tidyr)})
DATA <- Sys.getenv("DIDREP_DATA", unset="data"); OUT <- Sys.getenv("DIDREP_OUT", unset="output")
FIG  <- Sys.getenv("DIDREP_FIG",  unset="figures")
FIRST_OUT <- 2000; LAST <- 2020; Zc <- qnorm(0.975)+qnorm(0.80)
WEBB <- c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
REV <- c("revere massachusetts","portsmouth new hampshire","memphis tennessee",
         "fall river massachusetts","ramapo new york","springfield massachusetts")
code <- function(d){ d$agency.num <- as.numeric(factor(d$agency.id))
  d$yc <- ave(d$year.changed, d$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
  d$gname <- ifelse(is.na(d$yc)&d$no.req==0,0, ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc)); d }
d <- code(read.csv(file.path(DATA,"dta.csv")))
d <- code(d[!(d$agency.id %in% REV), ])
dw <- d[d$year>=FIRST_OUT & d$year<=LAST & !is.na(d$any.fatalities), ]
never_ids <- unique(dw$agency.id[dw$gname==0]); N_inf <- length(never_ids)
nochange_ids <- unique(dw$agency.id[dw$change.type=="No Change"])
cohorts <- sort(unique(dw$gname[dw$gname %in% 2001:2020]))
N_k <- sapply(cohorts, function(g) length(unique(dw$agency.id[dw$gname==g]))); names(N_k)<-cohorts

# ============================================================================
# ITEM 1
# ============================================================================
treated <- dw$gname %in% 2001:2020 & dw$year>=dw$gname
untreated <- (dw$gname==0) | (dw$gname %in% 2001:2020 & dw$year<dw$gname)
s1 <- feols(any.fatalities ~ 1 | agency.num + year, dw[untreated,], fixef.rm="none")
fe <- fixef(s1); A <- fe$agency.num; G <- fe$year
yhat0 <- function(ag,yr) unname(A[as.character(ag)] + G[as.character(yr)])
set.seed(0)
mp <- suppressWarnings(did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="gname",
   xformla=~1, control_group="notyettreated", base_period="varying", data=dw))
att_cs <- function(g,t){ i<-which(mp$group==g & mp$t==t); if(length(i)) mp$att[i] else NA_real_ }

# per-cohort s=1 gap + average-gap own/inherited decomposition
r1 <- lapply(cohorts, function(g){
  tg <- dw[dw$gname==g & dw$year==g, ]; pg <- dw[dw$gname==g & dw$year==(g-1), ]
  tau_imp_0 <- mean(tg$any.fatalities - yhat0(tg$agency.num, tg$year))
  cs_0  <- att_cs(g,g); gap0 <- tau_imp_0 - cs_0
  bb0   <- if(nrow(pg)) mean(pg$any.fatalities - yhat0(pg$agency.num, pg$year)) else NA_real_
  # cohort-average gap over all post periods (event times >=0)
  post <- dw[dw$gname==g & dw$year>=g, ]
  tau_imp_avg <- mean(post$any.fatalities - yhat0(post$agency.num, post$year))
  cs_avg <- mean(sapply(g:LAST, function(t) att_cs(g,t)), na.rm=TRUE)
  gap_avg <- tau_imp_avg - cs_avg
  data.frame(cohort=g, n_g=N_k[as.character(g)], tau_imp_0=tau_imp_0, cs_nyt_0=cs_0, gap_0=gap0,
    own_block=bb0, inherited_0=gap0-bb0, gap_avg=gap_avg, row.names=NULL) })
D1 <- do.call(rbind, r1)
write.csv(D1, file.path(OUT,"d6_gap_identity_s1.csv"), row.names=FALSE)
cat("=== ITEM 1: exact gap at s=1 (event time 0) ===\n")
print(transform(D1, tau_imp_0=round(tau_imp_0,4), cs_nyt_0=round(cs_nyt_0,4), gap_0=round(gap_0,4),
  own_block=round(own_block,4), inherited_0=round(inherited_0,4), gap_avg=round(gap_avg,4)), row.names=FALSE)
cat(sprintf("\nVERIFY 2001 gap_0 = %.2e (should be 0 to 1e-8)\n", D1$gap_0[D1$cohort==2001]))
cat(sprintf("VERIFY own+inherited = gap_0 exactly: max|own+inh - gap0| = %.2e\n",
    max(abs(D1$own_block + D1$inherited_0 - D1$gap_0))))

# Liu (2025) Props 2.1/2.2 weights  w_k = N_k / (sum_{j>=k} N_j + N_inf)
liu_w <- sapply(seq_along(cohorts), function(i){ k<-cohorts[i]; N_k[i] / (sum(N_k[cohorts>=k]) + N_inf) })
names(liu_w) <- cohorts
LIU <- data.frame(adjustment_cohort=cohorts, N_k=as.integer(N_k),
  denom=sapply(cohorts, function(k) sum(N_k[cohorts>=k])+N_inf), weight=round(liu_w,5), row.names=NULL)
write.csv(LIU, file.path(OUT,"d6_liu_weights.csv"), row.names=FALSE)
cat("\nLiu Props 2.1/2.2 adjustment-cohort weights  N_k/(sum_{j>=k}N_j + N_inf):  [reconstructed]\n")
print(LIU, row.names=FALSE)

# MAIN gap figure: own block bias + inherited, stacked, per cohort (s=1)
long <- pivot_longer(D1[,c("cohort","own_block","inherited_0")], cols=c("own_block","inherited_0"),
                     names_to="component", values_to="value")
long$component <- factor(long$component, levels=c("own_block","inherited_0"),
                         labels=c("own last-pre block bias","inherited (adjustment cohorts)"))
pg <- ggplot(long, aes(factor(cohort), value, fill=component)) +
  geom_col() +
  geom_point(data=D1, aes(factor(cohort), gap_0), inherit.aes=FALSE, shape=18, size=2.6, color="black") +
  geom_hline(yintercept=0, color="grey60") +
  scale_fill_manual(values=c("own last-pre block bias"="#0072B2","inherited (adjustment cohorts)"="#E69F00"), name=NULL) +
  labs(x="cohort", y=expression(hat(Delta)^Imp [g*",0"]*"   (imputation - CS-NYT at event time 0)"),
       title="CS-imputation gap at s=1: own block bias + inherited adjustment-cohort term",
       subtitle="diamond = total gap tau_imp(g,g) - tau_CS-NYT(g,g); 2001 = 0 by construction (single pre-period)",
       caption="Inherited term is attributed across later (adjustment) cohorts by the Liu (2025) weights N_k/(sum_{j>=k}N_j+N_inf); see d6_liu_weights.csv. Total is exact.") +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank(), legend.position="top",
    plot.caption=element_text(size=7, hjust=0, color="grey40"))
ggsave(file.path(FIG,"d6_gap_decomposition.pdf"), pg, width=9, height=5, device=cairo_pdf)
ggsave(file.path(FIG,"d6_gap_decomposition.png"), pg, width=9, height=5, dpi=200)
file.copy(file.path(FIG,"d6_gap_decomposition.png"), file.path(OUT,"d6_gap_decomposition.png"), overwrite=TRUE)

# move round-3 linear-slope plot to appendix name (if present)
for (ext in c("pdf","png")){ old<-file.path(FIG,paste0("d6_gap_predicted_vs_actual.",ext))
  if (file.exists(old)) file.copy(old, file.path(FIG,paste0("d6_appendix_gap_linear_slope.",ext)), overwrite=TRUE) }
if (file.exists(file.path(OUT,"d6_gap_predicted_vs_actual.png")))
  file.copy(file.path(OUT,"d6_gap_predicted_vs_actual.png"), file.path(OUT,"d6_appendix_gap_linear_slope.png"), overwrite=TRUE)

# ============================================================================
# ITEM 2: mde_g via slope-only restricted WCR for all cohorts
# ============================================================================
cr1_vcov <- function(X, resid, cl, XtXinv){ Gs<-unique(cl); p<-ncol(X); meat<-matrix(0,p,p)
  for(g in Gs){ ix<-cl==g; Xg<-X[ix,,drop=FALSE]; sc<-crossprod(Xg,resid[ix]); meat<-meat+tcrossprod(sc) }
  nG<-length(Gs); n<-nrow(X); adj<-(nG/(nG-1))*((n-1)/(n-p)); XtXinv %*% (adj*meat) %*% XtXinv }
wcr_slope_se <- function(s, B=999, seed=1){
  X<-model.matrix(~year*D,s); y<-s$any.fatalities; cl<-s$agency.num; k<-which(colnames(X)=="year:D")
  if(length(unique(cl))<3 || length(unique(s$D))<2) return(NA_real_)
  XtXinv<-tryCatch(solve(crossprod(X)),error=function(e) NULL); if(is.null(XtXinv)) return(NA_real_)
  Xr<-X[,-k,drop=FALSE]; br<-solve(crossprod(Xr))%*%crossprod(Xr,y); yhat_r<-Xr%*%br; ur<-y-yhat_r
  Gs<-unique(cl); set.seed(seed); bs<-numeric(B)
  for(b in 1:B){ wg<-sample(WEBB,length(Gs),replace=TRUE); names(wg)<-Gs
    ystar<-yhat_r+wg[as.character(cl)]*ur; bs[b]<-(XtXinv%*%crossprod(X,ystar))[k] }
  sd(bs) }
comp_ids <- function(est, g){ switch(est,
  CS=never_ids, stacked_cengiz=never_ids, imputation=union(never_ids, unique(dw$agency.id[dw$gname>g & dw$gname<3000])),
  stacked=nochange_ids) }
pre_sub <- function(g, tids, cids){ pre<-FIRST_OUT:(g-1)
  s<-dw[dw$year %in% pre & dw$agency.id %in% c(tids,cids),]; s$D<-as.integer(s$agency.id %in% tids); s }

TAB <- read.csv(file.path(OUT,"weight_checkability_table.csv"))
v2 <- TAB[TAB$sample=="droppers" & TAB$estimator %in% c("CS","stacked","stacked_cengiz","imputation"), ]
v2$n_treated_clusters <- v2$n_g
v2$mde_wcr <- NA_real_
for (i in seq_len(nrow(v2))){ g<-v2$cohort[i]; est<-v2$estimator[i]
  if (g==2000) next
  if ((g - FIRST_OUT) < 3) next                      # need >=3 pre-periods for a slope
  tids <- unique(dw$agency.id[dw$gname==g])
  s <- pre_sub(g, tids, comp_ids(est,g))
  se <- tryCatch(wcr_slope_se(s), error=function(e) NA_real_); v2$mde_wcr[i] <- Zc*se }
write.csv(v2[,c("estimator","cohort","n_g","n_treated_clusters","pre_g","weight","mde_g","mde_wcr","testable_g")],
          file.path(OUT,"d6_mde_wcr.csv"), row.names=FALSE)
cat("\n=== ITEM 2: mde_g CR1 vs WCR (droppers), + n_treated_clusters ===\n")
print(transform(v2[,c("estimator","cohort","n_treated_clusters","pre_g","mde_g","mde_wcr","testable_g")],
  mde_g=round(mde_g,4), mde_wcr=round(mde_wcr,4)), row.names=FALSE)
# rank changes per estimator (testable cohorts only)
cat("\nrank changes (testable cohorts, by increasing mde):\n")
for (est in c("CS","stacked","stacked_cengiz","imputation")){ p<-v2[v2$estimator==est & v2$testable_g & is.finite(v2$mde_wcr) & is.finite(v2$mde_g),]
  if(nrow(p)<2) next
  p$rank_cr1<-rank(p$mde_g); p$rank_wcr<-rank(p$mde_wcr); ch<-p[p$rank_cr1!=p$rank_wcr,]
  cat(sprintf("  %s: %s\n", est, if(nrow(ch)) paste(sprintf("%d (CR1 #%d->WCR #%d)",ch$cohort,ch$rank_cr1,ch$rank_wcr),collapse=", ") else "no rank changes")) }

# redraw main figure (v2) on the WCR mde x-axis
mk_v2wcr <- function(){
  sub <- v2; base_lab <- c(CS="A: CS (simple)", stacked="B: stacked (R2 FE)",
    stacked_cengiz="B2: stacked Cengiz +/-4 (NT)", imputation="C: imputation (BJS)")
  ord <- c("CS","stacked","stacked_cengiz","imputation")
  spear <- function(w,pg) if(length(w)>=3) suppressWarnings(cor(w,pg,method="spearman")) else NA
  stat <- sapply(ord, function(e){ p<-sub[sub$estimator==e,]
    sprintf("%s\nuntestable wt=%.2f  rho(w,pre_g)=%.2f", base_lab[e],
      sum(p$weight[!p$testable_g],na.rm=TRUE), spear(p$weight,p$pre_g)) })
  sub$estimator <- factor(sub$estimator, levels=ord, labels=stat)
  fin <- sub[sub$testable_g & is.finite(sub$mde_wcr),]; xmax<-max(fin$mde_wcr,na.rm=TRUE); x_un<-xmax*1.9
  sub$xpos <- ifelse(sub$testable_g & is.finite(sub$mde_wcr), sub$mde_wcr, x_un)
  g <- ggplot(sub, aes(xpos, weight)) +
    annotate("rect", xmin=xmax*1.35, xmax=xmax*2.6, ymin=-Inf, ymax=Inf, fill="grey88", alpha=0.6) +
    geom_point(aes(size=n_g, fill=testable_g), shape=21, color="grey25", alpha=0.9) +
    ggrepel::geom_text_repel(aes(label=cohort), size=2.3, max.overlaps=30, seed=1) +
    annotate("text", x=sqrt(xmax*1.35*xmax*2.6), y=Inf, label="not testable", angle=90, vjust=1.1, hjust=1, size=2.6, color="grey35") +
    scale_fill_manual(values=c(`TRUE`="#0072B2",`FALSE`="#D55E00"), name="testable") +
    scale_size_area(max_size=8, name=expression(n[g])) + scale_x_log10() +
    facet_wrap(~estimator, scales="free_y", nrow=1) +
    labs(x="minimum detectable pre-trend slope  mde_g  (WILD CLUSTER BOOTSTRAP SE; log; right strip = not testable)",
         y="cohort aggregation weight", title="Weight vs. checkability (droppers), x = WCR-based mde_g",
         caption="Panel B2 (Cengiz +/-4): cohorts with pre_g<4 (2001-2003) are excluded by the window, not down-weighted.") +
    theme_minimal(base_size=11) + theme(legend.position="right", panel.grid.minor=element_blank(),
      strip.text=element_text(size=8, lineheight=1.05), plot.caption=element_text(size=7,hjust=0,color="grey40"))
  ggsave(file.path(FIG,"fig_weight_checkability_d6_v2_wcrmde.pdf"), g, width=15, height=4.2, device=cairo_pdf)
  ggsave(file.path(FIG,"fig_weight_checkability_d6_v2_wcrmde.png"), g, width=15, height=4.2, dpi=200)
  file.copy(file.path(FIG,"fig_weight_checkability_d6_v2_wcrmde.png"), file.path(OUT,"fig_weight_checkability_d6_v2_wcrmde.png"), overwrite=TRUE) }
mk_v2wcr()

# ============================================================================
# ITEM 3: RI diagnostic + size-matched joint RI
# ============================================================================
sizevar <- if ("num.officers" %in% names(dw)) "num.officers" else "pop"
pre_stats <- function(ids, g){ pre<-FIRST_OUT:(g-1)
  s<-dw[dw$year %in% pre & dw$agency.id %in% ids,]
  c(mean_y=mean(s$any.fatalities), var_y=var(s$any.fatalities),
    size=mean(s[[sizevar]], na.rm=TRUE)) }
cat(sprintf("\n=== ITEM 3: pre-period diagnostics (size var = %s) ===\n", sizevar))
diag_rows <- list()
for (g in c(2009,2013)){ tids<-unique(dw$agency.id[dw$gname==g])
  ts<-pre_stats(tids,g); ns<-pre_stats(never_ids,g)
  diag_rows[[length(diag_rows)+1L]] <- data.frame(cohort=g, group="treated", mean_y=ts["mean_y"], var_y=ts["var_y"], size=ts["size"], row.names=NULL)
  diag_rows[[length(diag_rows)+1L]] <- data.frame(cohort=g, group="never-treated pool", mean_y=ns["mean_y"], var_y=ns["var_y"], size=ns["size"], row.names=NULL) }
DIAG <- do.call(rbind, diag_rows); write.csv(DIAG, file.path(OUT,"d6_ri_diagnostics.csv"), row.names=FALSE)
print(transform(DIAG, mean_y=round(mean_y,3), var_y=round(var_y,4), size=round(size,0)), row.names=FALSE)

# size-matched joint RI: placebo treated = nearest never-treated agencies by size
agsize <- tapply(dw[[sizevar]], dw$agency.id, function(x) mean(x,na.rm=TRUE))
ri_joint_matched <- function(g, n_g, B=499, seed=3, matched=TRUE){
  pre<-FIRST_OUT:(g-1); base<-dw[dw$year %in% pre & dw$agency.id %in% never_ids,]
  ids<-unique(base$agency.id); tids<-unique(dw$agency.id[dw$gname==g])
  tsz<-mean(agsize[tids],na.rm=TRUE)
  # observed joint Wald (treated vs never)
  so<-base; so<-rbind(so, transform(dw[dw$year %in% pre & dw$agency.id %in% tids,]))
  so$D<-as.integer(so$agency.id %in% tids)
  Wobs<-tryCatch(unname(fixest::wald(feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,so,cluster=~agency.num),"year::",print=FALSE)$stat),error=function(e)NA)
  set.seed(seed); Wb<-numeric(B)
  pool_order <- names(sort(abs(agsize[ids]-tsz)))        # size-ranked never-treated
  for(b in 1:B){
    if(matched){ cand<-head(pool_order, max(n_g*3, n_g+5)); tr<-sample(cand, n_g) }
    else tr<-sample(ids, n_g)
    dd<-base; dd$D<-as.integer(dd$agency.id %in% tr)
    mj<-tryCatch(feols(any.fatalities~i(year,D,ref=g-1)|agency.num+year,dd),error=function(e)NULL)
    Wb[b]<- if(is.null(mj)) NA else tryCatch(unname(fixest::wald(mj,"year::",print=FALSE)$stat),error=function(e)NA) }
  (1+sum(Wb>=Wobs,na.rm=TRUE))/(1+sum(!is.na(Wb))) }
cat("\n=== ITEM 3: joint RI, unmatched vs size-matched never-treated placebos (alongside WCR) ===\n")
rob <- read.csv(file.path(OUT,"d6_pretests_2009_2013_robust.csv"))
ri_rows <- list()
for (g in c(2009,2013)){ n_g<-length(unique(dw$agency.id[dw$gname==g]))
  wcr_j <- mean(rob$joint_p_WCR[rob$cohort==g]);
  ri_u <- ri_joint_matched(g,n_g,matched=FALSE); ri_m <- ri_joint_matched(g,n_g,matched=TRUE)
  ri_rows[[length(ri_rows)+1L]] <- data.frame(cohort=g, joint_p_WCR_mean=round(wcr_j,3),
    joint_p_RI_unmatched=round(ri_u,3), joint_p_RI_sizematched=round(ri_m,3), row.names=NULL) }
RI <- do.call(rbind, ri_rows); write.csv(RI, file.path(OUT,"d6_ri_sizematched.csv"), row.names=FALSE)
print(RI, row.names=FALSE)

cat("\nWrote d6_gap_identity_s1.csv, d6_liu_weights.csv, d6_gap_decomposition.*, d6_mde_wcr.csv,\n")
cat("  fig_weight_checkability_d6_v2_wcrmde.*, d6_ri_diagnostics.csv, d6_ri_sizematched.csv\n")
