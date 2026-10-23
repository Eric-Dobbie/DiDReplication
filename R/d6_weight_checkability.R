# =============================================================================
# d6_weight_checkability.R  --  weight vs. checkability, one panel per estimator.
#
# Core claim (on P&P): under every modern estimator, the cohorts carrying the
# most aggregation weight are the cohorts whose parallel-trends assumption is
# least checkable. x = checkability (mde_g, and a companion with pre_g);
# y = cohort-level aggregation weight; one small-multiple panel per estimator.
#
# Estimators (R5: never-treated controls, no covariates, universal base, no bin):
#   A  CS      : did::att_gt + aggte(type="simple"); cohort weight = sum of atom
#               weights over that cohort's post-periods (exact linear weights,
#               recovered by differentiating overall.att wrt each att(g,t) atom).
#   B  stacked : corrected agency x stack + year x cohort FE (R2); w_s = V_s/V
#               (Prop. 1), reused from fwl_decomp.R (deck) / rebuilt (droppers).
#   C  imputation (BJS/Gardner): Y(0) imputed from agency+year FE on untreated
#               obs; ATT = mean tau over treated obs; simple weights w_it=1/N_1
#               so cohort weight = treated-obs share. Requires ABSORBING
#               treatment -> run on the droppers sample only.
#
# did2s / didimputation are NOT installable here (CRAN unreachable through the
# proxy; R 4.3.3 predates the binaries). Panel C is a from-scratch BJS/Gardner
# imputation, which under R5 (no covariates) targets the identical estimand; the
# ATT it returns is used for validation 2 in place of did_imputation's.
#
# Treatment reversals: imputation needs absorbing treatment. The DROPPERS sample
# excludes the 6 "other-direction" (Adopted) in-window switchers -- revere 2002,
# portsmouth 2003, memphis 2004, fall river 2012, ramapo 2018, springfield MA
# 2018 -- which also removes the two absorbed reversals (memphis 2009,
# portsmouth 2016, the second switch of memphis/portsmouth). Two figure versions:
# (i) A/B on the deck sample, C on droppers; (ii) all three on droppers.
#
# Cohort identity: stacked cohort 2002 (zero clean controls) is NOT CS cohort
# 2002 (41 never-treated controls). Points are labelled by estimator-specific
# cohort; 2002 in Panel B is footnoted.
#
# FILE NAMES: the brief asks for fig_weight_checkability.{pdf,png}; that name is
# already a committed deliverable of weight_checkability.R (a different analysis).
# To avoid clobbering it, D6 figures use a d6_ prefix. FLAGGED in the summary.
# =============================================================================
suppressMessages({library(did); library(fixest); library(lfe); library(ggplot2); library(ggrepel)})
DATA <- Sys.getenv("DIDREP_DATA", unset="data"); OUT <- Sys.getenv("DIDREP_OUT", unset="output")
FIG  <- Sys.getenv("DIDREP_FIG",  unset="figures")
FIRST_OUT <- 2000; LAST <- 2020; WINDOW <- 4
Zc <- qnorm(0.975) + qnorm(0.80)                      # 1.95996 + 0.84162 = 2.80158

REVERSERS <- c("revere massachusetts","portsmouth new hampshire","memphis tennessee",
               "fall river massachusetts","ramapo new york","springfield massachusetts")

code <- function(d){
  d$agency.num <- as.numeric(factor(d$agency.id))
  d$yc <- ave(d$year.changed, d$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
  d$gname <- ifelse(is.na(d$yc)&d$no.req==0,0, ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc))
  d }
dta  <- code(read.csv(file.path(DATA,"dta.csv")))
dta_drop <- code(dta[!(dta$agency.id %in% REVERSERS), ])     # absorbing-droppers sample

# ============================================================================
# checkability: pre_g, post_g, mde_g (estimator-specific comparison group)
# ============================================================================
pre_of  <- function(g) g - FIRST_OUT               # observed untreated pre-periods
post_of <- function(g) LAST - g + 1
never_ids <- function(d) unique(d$agency.id[d$gname==0])
nochange_ids <- function(d) unique(d$agency.id[d$change.type=="No Change"])
notyet_ids <- function(d,g) unique(d$agency.id[d$gname>g & d$gname<3000])

# minimum detectable linear pre-trend slope: Zc * clustered SE of the differential
# (treated - comparison) pre-period time trend. Pre-periods only; never post.
mde_g <- function(treated_ids, comp_ids, g, dsrc){
  pre <- FIRST_OUT:(g-1)
  if (length(pre) < 3) return(NA_real_)                       # need >=3 pre for residual df
  if (!length(treated_ids) || !length(comp_ids)) return(NA_real_)
  sub <- dsrc[dsrc$year %in% pre & !is.na(dsrc$any.fatalities) &
              dsrc$agency.id %in% c(treated_ids, comp_ids), ]
  sub$D <- as.integer(sub$agency.id %in% treated_ids)
  if (length(unique(sub$D))<2 || length(unique(sub$year))<3) return(NA_real_)
  m <- tryCatch(feols(any.fatalities ~ year*D, data=sub, cluster=~agency.num), error=function(e) NULL)
  if (is.null(m)) return(NA_real_)
  ct <- tryCatch(coeftable(m), error=function(e) NULL)
  if (is.null(ct) || !("year:D" %in% rownames(ct))) return(NA_real_)
  Zc * ct["year:D","Std. Error"] }

# testable flag + mde, for one cohort under one estimator
check_cohort <- function(g, estimator, treated_ids, dsrc, has_controls=TRUE){
  pg <- pre_of(g)
  comp <- switch(estimator,
    CS         = never_ids(dsrc),
    imputation = union(never_ids(dsrc), notyet_ids(dsrc,g)),
    stacked    = if (has_controls) nochange_ids(dsrc) else character(0))
  cg <- switch(estimator, CS="never-treated", imputation="never + not-yet-treated at g",
               stacked=if(has_controls) "stack clean controls (No Change)" else "none")
  md <- if (!has_controls) NA_real_ else mde_g(treated_ids, comp, g, dsrc)
  testable <- (pg>=3) && has_controls && is.finite(md) && length(treated_ids)>0
  list(pre_g=pg, post_g=post_of(g), mde_g=md, testable=testable, comparison=cg) }

# ============================================================================
# Panel A: CS simple -- exact cohort weights via att_gt + aggte gradient
# ============================================================================
cs_panel <- function(d, sample, cg="nevertreated"){
  set.seed(0)
  mp <- suppressWarnings(did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",
        gname="gname", xformla=~1, control_group=cg, clustervars="agency.num",
        base_period="universal", data=d))
  base <- aggte(mp, type="simple")$overall.att
  n <- length(mp$att); w <- numeric(n); dd <- 1e-4
  for (k in seq_len(n)){ m2<-mp; m2$att[k]<-m2$att[k]+dd; w[k]<-(aggte(m2,type="simple")$overall.att-base)/dd }
  gg <- mp$group
  wg  <- tapply(w, gg, sum); tau <- tapply(w*mp$att, gg, sum)/wg
  cohorts <- as.integer(names(wg))
  rows <- lapply(seq_along(cohorts), function(i){ g<-cohorts[i]
    tids <- unique(d$agency.id[d$gname==g & !is.na(d$any.fatalities)])
    ck <- check_cohort(g, "CS", tids, d)
    data.frame(estimator="CS", sample=sample, cohort=g, n_g=length(tids),
      pre_g=ck$pre_g, post_g=ck$post_g, weight=as.numeric(wg[i]), mde_g=ck$mde_g,
      testable_g=ck$testable, comparison_group=ck$comparison,
      atom_or_stack_estimate=as.numeric(tau[i]), row.names=NULL) })
  list(tab=do.call(rbind, rows), overall=base) }

# ============================================================================
# Panel B: stacked (R2 FE), w_s = V_s/V. Reuse fwl_decomp logic.
# ============================================================================
decomp_stack <- function(d){                     # d = reconstructed stack, estimation sample
  d$agency_stack <- paste(d$agency.id, d$cohort, sep="__")
  stacks <- sort(unique(d$cohort))
  r1 <- felm(no.req ~ 1 | agency_stack + year.cohort, data=d); d$vhat <- as.numeric(r1$residuals[,1])
  Vs <- sapply(stacks, function(s){ ix<-d$cohort==s; sum(d$vhat[ix]^2) }); V<-sum(Vs); ws<-Vs/V
  betas <- sapply(stacks, function(s){ ds<-d[d$cohort==s,]
    m<-tryCatch(felm(any.fatalities ~ no.req | agency.id + year, data=ds), error=function(e) NULL)
    if(is.null(m)) return(NA_real_); b<-coef(m)["no.req"]; if(is.na(b)) NA_real_ else unname(b) })
  ng <- sapply(stacks, function(s) length(unique(d$agency.id[d$cohort==s & d$treat==1])))
  nc <- sapply(stacks, function(s) length(unique(d$agency.id[d$cohort==s & d$treat==0])))
  data.frame(stack=stacks, V_s=Vs, w_s=ws, beta_s=betas, n_g=ng, n_c=nc, row.names=NULL) }

build_stack <- function(d){
  cohorts <- sort(unique(d$year.changed[!is.na(d$year.changed) & d$year.changed>=FIRST_OUT]))
  nochange <- d[d$change.type=="No Change", ]; parts <- list()
  for (g in cohorts){ tr<-d[!is.na(d$year.changed)&d$year.changed==g,]; tr$cohort<-g; tr$treat<-1L
    parts[[length(parts)+1L]]<-tr
    if (g-WINDOW>=FIRST_OUT){ ct<-nochange; ct$cohort<-g; ct$treat<-0L; parts[[length(parts)+1L]]<-ct } }
  recon <- do.call(rbind, parts)
  recon$year.cohort <- as.numeric(paste0(recon$year, recon$cohort)); recon }

stacked_panel <- function(d, sample, deck_csv=NULL, stacked_fatal=NULL){
  if (!is.null(deck_csv)){                       # deck: reuse committed decomposition
    fw <- read.csv(deck_csv)
    ng <- sapply(fw$stack, function(g) length(unique(stacked_fatal$agency.id[stacked_fatal$cohort==g & stacked_fatal$treat==1])))
    nc <- sapply(fw$stack, function(g) length(unique(stacked_fatal$agency.id[stacked_fatal$cohort==g & stacked_fatal$treat==0])))
    dec <- data.frame(stack=fw$stack, w_s=fw$w_s, beta_s=fw$beta_s, n_g=ng, n_c=nc)
    dsrc <- dta
  } else {                                       # droppers: rebuild + decompose
    rec <- build_stack(d); est <- rec[!is.na(rec$any.fatalities), ]
    dec <- decomp_stack(est); dsrc <- d }
  overall <- sum(ifelse(dec$w_s>0 & !is.na(dec$beta_s), dec$w_s*dec$beta_s, 0))
  rows <- lapply(seq_len(nrow(dec)), function(i){ g<-dec$stack[i]
    tids <- unique(dsrc$agency.id[!is.na(dsrc$year.changed) & dsrc$year.changed==g])
    ck <- check_cohort(g, "stacked", tids, dsrc, has_controls=(dec$n_c[i]>0))
    data.frame(estimator="stacked", sample=sample, cohort=g, n_g=dec$n_g[i],
      pre_g=ck$pre_g, post_g=ck$post_g, weight=dec$w_s[i], mde_g=ck$mde_g,
      testable_g=ck$testable, comparison_group=ck$comparison,
      atom_or_stack_estimate=dec$beta_s[i], row.names=NULL) })
  list(tab=do.call(rbind, rows), overall=overall) }

# ============================================================================
# Panel B2: Cengiz +/-4 stacks with NEVER-TREATED controls (robustness on B)
#   balanced +/-4 window, cohorts 2004..2016 (full window in 2000-2020),
#   controls = never-treated (gname==0, 41), Cengiz FE = agency x stack +
#   year x stack (= the R2 interactive FE). V_s/V via the same decomposition.
# ============================================================================
cengiz_panel <- function(d, sample){
  K <- WINDOW; never <- unique(d$agency.id[d$gname==0]); nt <- d[d$agency.id %in% never, ]
  cohorts <- sort(unique(d$year.changed[!is.na(d$year.changed) & d$year.changed>=FIRST_OUT]))
  cohorts <- cohorts[cohorts-K>=FIRST_OUT & cohorts+K<=LAST]                 # 2004..2016
  parts <- list()
  for (g in cohorts){ tr<-d[!is.na(d$year.changed)&d$year.changed==g & d$year>=g-K & d$year<=g+K,]
    if (nrow(tr)==0) next; tr$cohort<-g; tr$treat<-1L
    ct<-nt[nt$year>=g-K & nt$year<=g+K,]; ct$cohort<-g; ct$treat<-0L
    parts[[length(parts)+1L]]<-rbind(tr,ct) }
  rec <- do.call(rbind, parts); rec$year.cohort <- as.numeric(paste0(rec$year, rec$cohort))
  est <- rec[!is.na(rec$any.fatalities), ]; dec <- decomp_stack(est)
  overall <- sum(ifelse(dec$w_s>0 & !is.na(dec$beta_s), dec$w_s*dec$beta_s, 0))
  rows <- lapply(seq_len(nrow(dec)), function(i){ g<-dec$stack[i]
    tids <- unique(d$agency.id[!is.na(d$year.changed) & d$year.changed==g])
    ck <- check_cohort(g, "CS", tids, d, has_controls=(dec$n_c[i]>0))       # NT comparison, like CS
    data.frame(estimator="stacked_cengiz", sample=sample, cohort=g, n_g=dec$n_g[i],
      pre_g=ck$pre_g, post_g=ck$post_g, weight=dec$w_s[i], mde_g=ck$mde_g,
      testable_g=ck$testable, comparison_group="never-treated (Cengiz +/-4)",
      atom_or_stack_estimate=dec$beta_s[i], row.names=NULL) })
  list(tab=do.call(rbind, rows), overall=overall) }

# ============================================================================
# Panel C: BJS/Gardner imputation (droppers only, absorbing)
# ============================================================================
imputation_panel <- function(d, sample){
  d <- d[d$year>=FIRST_OUT & d$year<=LAST & !is.na(d$any.fatalities), ]
  treated <- d$gname %in% 2001:2020 & d$year >= d$gname            # absorbing drop, year>=onset
  untreated <- (d$gname==0) | (d$gname %in% 2001:2020 & d$year < d$gname)  # never + not-yet
  u <- d[untreated, ]; t <- d[treated, ]
  # fixef.rm="none": a lone untreated obs still identifies alpha_i (= Y - gamma_t);
  # the default "perfect_fit" would drop it as a singleton and lose the cohort (e.g. 2001/Denver).
  s1 <- feols(any.fatalities ~ 1 | agency.num + year, data=u, fixef.rm="none")
  t$yhat0 <- predict(s1, newdata=t)                                # NA only if agency truly has no untreated obs
  dropped_na <- sum(is.na(t$yhat0)); t <- t[!is.na(t$yhat0), ]
  t$tau <- t$any.fatalities - t$yhat0
  N1 <- nrow(t); overall <- mean(t$tau)
  cohorts <- sort(unique(t$gname))
  rows <- lapply(cohorts, function(g){ tg <- t[t$gname==g, ]
    tids <- unique(d$agency.id[d$gname==g])
    ck <- check_cohort(g, "imputation", tids, d)
    data.frame(estimator="imputation", sample=sample, cohort=g, n_g=length(unique(tg$agency.id)),
      pre_g=ck$pre_g, post_g=ck$post_g, weight=nrow(tg)/N1, mde_g=ck$mde_g,
      testable_g=ck$testable, comparison_group=ck$comparison,
      atom_or_stack_estimate=mean(tg$tau), row.names=NULL) })
  list(tab=do.call(rbind, rows), overall=overall, N1=N1, dropped_na=dropped_na) }

# ============================================================================
# RUN
# ============================================================================
sf <- read.csv(file.path(DATA,"stacked_fatal.csv"))
csA_deck <- cs_panel(dta,      "deck")
csA_drop <- cs_panel(dta_drop, "droppers")
stB_deck <- stacked_panel(dta, "deck", deck_csv=file.path(OUT,"fwl_decomp_unweighted.csv"), stacked_fatal=sf)
stB_drop <- stacked_panel(dta_drop, "droppers")
cz_deck  <- cengiz_panel(dta,      "deck")
cz_drop  <- cengiz_panel(dta_drop, "droppers")
imC_drop <- imputation_panel(dta_drop, "droppers")

TAB <- rbind(csA_deck$tab, csA_drop$tab, stB_deck$tab, stB_drop$tab,
             cz_deck$tab, cz_drop$tab, imC_drop$tab)
write.csv(TAB, file.path(OUT,"weight_checkability_table.csv"), row.names=FALSE)

cat("=============================================================\n")
cat("overall ATTs: CS deck", round(csA_deck$overall,6), "| CS drop", round(csA_drop$overall,6),
    "| stacked deck", round(stB_deck$overall,6), "| stacked drop", round(stB_drop$overall,6),
    "\n              Cengiz+/-4 NT deck", round(cz_deck$overall,6), "| Cengiz drop", round(cz_drop$overall,6),
    "| imputation drop", round(imC_drop$overall,6), "\n")

# ============================================================================
# VALIDATION
# ============================================================================
cat("\n=== VALIDATION ===\n")
vald <- function(est, samp, target=NA){ d <- TAB[TAB$estimator==est & TAB$sample==samp, ]
  ws <- sum(d$weight, na.rm=TRUE); att <- sum(d$weight*d$atom_or_stack_estimate, na.rm=TRUE)
  cat(sprintf("  %-11s %-8s: sum(w)=%.8f ; sum(w*est)=%.6f%s\n", est, samp, ws, att,
      if(is.finite(target)) sprintf(" (target %.3f, diff %.2e)", target, att-target) else "")) }
vald("CS","deck",0.045); vald("CS","droppers")
vald("stacked","deck",-0.099); vald("stacked","droppers")
vald("stacked_cengiz","deck", cz_deck$overall); vald("stacked_cengiz","droppers", cz_drop$overall)
vald("imputation","droppers", imC_drop$overall)
cat(sprintf("  [1] weights sum to 1 within each estimator x sample (above).\n"))
cat(sprintf("  [4] cohort 2000 present in CS? %s ; in imputation? %s  (dropped: already-treated in first obs year 2000 -> no untreated pre-period)\n",
    2000 %in% TAB$cohort[TAB$estimator=="CS"], 2000 %in% TAB$cohort[TAB$estimator=="imputation"]))
# [3] CS vs imputation simple weights on droppers
csw <- csA_drop$tab[,c("cohort","weight")]; imw <- imC_drop$tab[,c("cohort","weight")]
cmp <- merge(csw, imw, by="cohort", suffixes=c(".cs",".imp"), all=TRUE)
cmp$absdiff <- abs(cmp$weight.cs - cmp$weight.imp)
cat("  [3] CS vs imputation cohort weights (droppers):\n"); print(transform(cmp,
  weight.cs=round(weight.cs,4), weight.imp=round(weight.imp,4), absdiff=round(absdiff,4)), row.names=FALSE)
cat(sprintf("      max abs diff = %.4f\n", max(cmp$absdiff, na.rm=TRUE)))

# ============================================================================
# FIGURES
# ============================================================================
spear <- function(w, pg) if (length(w)>=3) suppressWarnings(cor(w, pg, method="spearman")) else NA_real_

make_fig <- function(sub, file, xvar, xlab, title, width=12){
  base_lab <- c(CS="A: CS (simple)", stacked="B: stacked (R2 FE)",
                stacked_cengiz="B2: stacked Cengiz +/-4 (NT)", imputation="C: imputation (BJS)")
  ord <- intersect(c("CS","stacked","stacked_cengiz","imputation"), unique(sub$estimator))
  # per-panel stats embedded in the strip title (ATT + untestable weight share + rho)
  stat <- sapply(ord, function(e){ p<-sub[sub$estimator==e,]
    sprintf("%s\nATT=%.3f  untestable wt=%.2f  rho(w,pre_g)=%.2f", base_lab[e],
      sum(p$weight*p$atom_or_stack_estimate, na.rm=TRUE), sum(p$weight[!p$testable_g], na.rm=TRUE),
      spear(p$weight, p$pre_g)) })
  sub$estimator <- factor(sub$estimator, levels=ord, labels=stat)
  # x placement: untestable cohorts -> right strip at xmax*1.9 (mde fig only)
  fin <- sub[sub$testable_g & is.finite(sub[[xvar]]), ]
  if (xvar=="mde_g"){
    xmax <- max(fin$mde_g, na.rm=TRUE); x_un <- xmax*1.9
    sub$xpos <- ifelse(sub$testable_g & is.finite(sub$mde_g), sub$mde_g, x_un)
    strip_lo <- xmax*1.35; strip_hi <- xmax*2.6
  } else { sub$xpos <- sub[[xvar]]; x_un <- NA }
  g <- ggplot(sub, aes(xpos, weight)) +
    { if (xvar=="mde_g") annotate("rect", xmin=strip_lo, xmax=strip_hi, ymin=-Inf, ymax=Inf, fill="grey88", alpha=0.6) } +
    geom_point(aes(size=n_g, fill=testable_g), shape=21, color="grey25", alpha=0.9) +
    ggrepel::geom_text_repel(aes(label=cohort), size=2.3, max.overlaps=30, seed=1) +
    scale_fill_manual(values=c(`TRUE`="#0072B2", `FALSE`="#D55E00"), name="testable") +
    scale_size_area(max_size=8, name=expression(n[g])) +
    facet_wrap(~estimator, scales="free_y", nrow=1) +
    labs(x=xlab, y="cohort aggregation weight", title=title) +
    theme_minimal(base_size=11) + theme(legend.position="right", panel.grid.minor=element_blank(),
      strip.text=element_text(size=8, lineheight=1.05))
  if (xvar=="mde_g"){
    g <- g + scale_x_log10() +
      annotate("text", x=sqrt(strip_lo*strip_hi), y=Inf, label="not testable", angle=90, vjust=1.1, hjust=1, size=2.6, color="grey35")
  }
  ggsave(file.path(FIG, paste0(file,".pdf")), g, width=width, height=4.2, device=cairo_pdf)
  ggsave(file.path(FIG, paste0(file,".png")), g, width=width, height=4.2, dpi=200)
  file.copy(file.path(FIG,paste0(file,".png")), file.path(OUT,paste0(file,".png")), overwrite=TRUE)
  invisible(g) }

# version (i): A/B/B2 deck, C droppers ; version (ii): all droppers
v1 <- rbind(csA_deck$tab, stB_deck$tab, cz_deck$tab, imC_drop$tab)
v2 <- rbind(csA_drop$tab, stB_drop$tab, cz_drop$tab, imC_drop$tab)
make_fig(v1, "fig_weight_checkability_d6_v1", "mde_g",
  "minimum detectable pre-trend slope  mde_g  (log; right strip = not testable)",
  "Weight vs. checkability  (i): CS/stacked/Cengiz on deck sample, imputation on droppers", width=15)
make_fig(v2, "fig_weight_checkability_d6_v2", "mde_g",
  "minimum detectable pre-trend slope  mde_g  (log; right strip = not testable)",
  "Weight vs. checkability  (ii): all estimators on droppers sample", width=15)
make_fig(v2, "fig_weight_checkability_d6_preg", "pre_g",
  "observed untreated pre-periods  pre_g = g - 2000",
  "Companion (pre_g axis): all estimators on droppers sample", width=15)

# ---- weight movement between samples (report) --------------------------------
cat("\n=== cohort-weight movement deck -> droppers ===\n")
for (est in c("CS","stacked")){
  a <- TAB[TAB$estimator==est & TAB$sample=="deck", c("cohort","weight")]
  b <- TAB[TAB$estimator==est & TAB$sample=="droppers", c("cohort","weight")]
  m <- merge(a,b,by="cohort",suffixes=c(".deck",".drop"),all=TRUE); m[is.na(m)]<-0
  m$delta <- m$weight.drop - m$weight.deck
  cat(sprintf("  %s: max |delta weight| = %.4f (sum|delta| = %.4f)\n", est, max(abs(m$delta)), sum(abs(m$delta))))
  big <- m[order(-abs(m$delta)),][1:3,]
  cat("    largest shifts:", paste(sprintf("%d: %+.3f", big$cohort, big$delta), collapse=" ; "), "\n") }

# ---- rho per panel per version (report) --------------------------------------
cat("\n=== Spearman rho(weight, pre_g) ===\n")
for (vn in c("v1","v2")){ vv <- get(vn)
  for (est in unique(vv$estimator)){ p <- vv[vv$estimator==est, ]
    cat(sprintf("  %-4s %-11s rho=%.2f ; untestable weight share=%.2f\n", vn, est,
        spear(p$weight,p$pre_g), sum(p$weight[!p$testable_g],na.rm=TRUE))) } }

# ============================================================================
# ITEM 3: bridge the CS-imputation gap on droppers:  CS-NT -> CS-NYT -> imputation
# ============================================================================
cat("\n=== BRIDGE: CS-NT -> CS-NYT -> imputation (droppers) ===\n")
csNYT <- cs_panel(dta_drop, "droppers", cg="notyettreated")
cat(sprintf("overall ATT:  CS-NT %.4f  ->  CS-NYT %.4f  ->  imputation %.4f\n",
    csA_drop$overall, csNYT$overall, imC_drop$overall))
b_nt  <- csA_drop$tab[,c("cohort","atom_or_stack_estimate")]; names(b_nt)[2]<-"CS_NT"
b_nyt <- csNYT$tab[,c("cohort","atom_or_stack_estimate")];    names(b_nyt)[2]<-"CS_NYT"
b_imp <- imC_drop$tab[,c("cohort","atom_or_stack_estimate")]; names(b_imp)[2]<-"imputation"
wts   <- csA_drop$tab[,c("cohort","weight")]; names(wts)[2]<-"weight_CSsimple"
bridge <- Reduce(function(a,b) merge(a,b,by="cohort",all=TRUE), list(wts,b_nt,b_nyt,b_imp))
write.csv(bridge, file.path(OUT,"d6_bridge_cs_imputation.csv"), row.names=FALSE)
cat("per-cohort estimates (tau_g):\n")
print(transform(bridge, weight_CSsimple=round(weight_CSsimple,3), CS_NT=round(CS_NT,3),
  CS_NYT=round(CS_NYT,3), imputation=round(imputation,3)), row.names=FALSE)

# ============================================================================
# ITEM 4: pre-tests for the highest-weight testable cohorts (stacked 2009, 2013)
#   joint pre-trend Wald + linear-trend slope, per estimator's own comparison grp
# ============================================================================
pretest <- function(g, treated_ids, comp_ids, dsrc, label){
  pre <- FIRST_OUT:(g-1)
  sub <- dsrc[dsrc$year %in% pre & !is.na(dsrc$any.fatalities) & dsrc$agency.id %in% c(treated_ids,comp_ids), ]
  sub$D <- as.integer(sub$agency.id %in% treated_ids)
  out <- data.frame(estimator=label, cohort=g, n_pre=length(pre), n_treated=length(treated_ids),
    slope=NA_real_, slope_se=NA_real_, slope_p=NA_real_, mde=NA_real_,
    joint_stat=NA_real_, joint_df=NA_real_, joint_p=NA_real_, row.names=NULL)
  if (length(pre)<3 || !length(comp_ids) || length(unique(sub$D))<2) return(out)
  ml <- tryCatch(feols(any.fatalities ~ year*D, data=sub, cluster=~agency.num), error=function(e) NULL)
  if (!is.null(ml)){ ct<-coeftable(ml); if("year:D" %in% rownames(ct)){
    out$slope<-ct["year:D","Estimate"]; out$slope_se<-ct["year:D","Std. Error"]
    out$slope_p<-ct["year:D","Pr(>|t|)"]; out$mde<-Zc*ct["year:D","Std. Error"] } }
  mj <- tryCatch(feols(any.fatalities ~ i(year, D, ref=g-1) | agency.num + year, data=sub, cluster=~agency.num), error=function(e) NULL)
  if (!is.null(mj)){ w <- tryCatch(fixest::wald(mj, "year::", print=FALSE), error=function(e) NULL)
    if (!is.null(w)){ out$joint_stat<-unname(w$stat); out$joint_df<-unname(w$df1); out$joint_p<-unname(w$p) } }
  out }

cat("\n=== PRE-TESTS: stacked 2009 & 2013, each estimator's own comparison group ===\n")
pt_rows <- list()
for (g in c(2009, 2013)){
  tids <- unique(dta_drop$agency.id[!is.na(dta_drop$year.changed) & dta_drop$year.changed==g])
  tids_cs <- unique(dta_drop$agency.id[dta_drop$gname==g])
  pt_rows[[length(pt_rows)+1L]] <- pretest(g, tids_cs, never_ids(dta_drop), dta_drop, "CS (never-treated)")
  pt_rows[[length(pt_rows)+1L]] <- pretest(g, tids, nochange_ids(dta_drop), dta_drop, "stacked (No-Change controls)")
  pt_rows[[length(pt_rows)+1L]] <- pretest(g, tids_cs, union(never_ids(dta_drop), notyet_ids(dta_drop,g)), dta_drop, "imputation (never+not-yet)")
  pt_rows[[length(pt_rows)+1L]] <- pretest(g, tids, never_ids(dta_drop), dta_drop, "Cengiz (never-treated)")
}
PT <- do.call(rbind, pt_rows)
PT$slope_pass <- ifelse(is.na(PT$slope_p), NA, PT$slope_p>0.05)      # pass = no detectable pre-trend
PT$joint_pass <- ifelse(is.na(PT$joint_p), NA, PT$joint_p>0.05)
write.csv(PT, file.path(OUT,"d6_pretests_2009_2013.csv"), row.names=FALSE)
print(transform(PT, slope=round(slope,4), slope_se=round(slope_se,4), slope_p=round(slope_p,3),
  mde=round(mde,4), joint_stat=round(joint_stat,2), joint_p=round(joint_p,3)), row.names=FALSE)

cat("\nWrote weight_checkability_table.csv, d6_bridge_cs_imputation.csv, d6_pretests_2009_2013.csv,\n")
cat("  fig_weight_checkability_d6_{v1,v2,preg}.*\n")
