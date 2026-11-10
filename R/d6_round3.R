# =============================================================================
# d6_round3.R  --  D6 round-3 follow-ups.
#
# 1. Re-run the 2009/2013 pre-tests (linear slope + joint pre-trend) with
#    (a) restricted WILD CLUSTER BOOTSTRAP (Webb 6-pt weights, B=999 slope /
#        499 joint), and (b) RANDOMIZATION INFERENCE over never-treated placebo
#    cohorts (B=999), reported side by side with the CR1 (cluster-robust t / Wald)
#    p-values from d6_pretests_2009_2013.csv.
# 2. Per-cohort (droppers) PREDICTED CS-imputation gap = slope_g * (pre_g-1)/2
#    (slope_g = estimated differential linear pre-trend vs never-treated) vs the
#    ACTUAL tau_imp - tau_CS; table + predicted-vs-actual plot.
#
# fwildclusterboot is not installable here (CRAN unreachable; R 4.3.3), so WCR is
# implemented manually (restricted bootstrap-t for the slope; restricted Wald
# bootstrap for the joint test).
# =============================================================================
suppressMessages({library(fixest); library(ggplot2); library(ggrepel)})
DATA <- Sys.getenv("DIDREP_DATA", unset="data"); OUT <- Sys.getenv("DIDREP_OUT", unset="output")
FIG  <- Sys.getenv("DIDREP_FIG",  unset="figures")
FIRST_OUT <- 2000; LAST <- 2020
REVERSERS <- c("revere massachusetts","portsmouth new hampshire","memphis tennessee",
               "fall river massachusetts","ramapo new york","springfield massachusetts")
WEBB <- c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))

code <- function(d){ d$agency.num <- as.numeric(factor(d$agency.id))
  d$yc <- ave(d$year.changed, d$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
  d$gname <- ifelse(is.na(d$yc)&d$no.req==0,0, ifelse(is.na(d$yc)&d$no.req==1,1987,d$yc)); d }
dta_drop <- code(read.csv(file.path(DATA,"dta.csv")))
dta_drop <- code(dta_drop[!(dta_drop$agency.id %in% REVERSERS), ])
never_ids <- unique(dta_drop$agency.id[dta_drop$gname==0])
nochange_ids <- unique(dta_drop$agency.id[dta_drop$change.type=="No Change"])
notyet_ids <- function(g) unique(dta_drop$agency.id[dta_drop$gname>g & dta_drop$gname<3000])

pre_sub <- function(g, treated_ids, comp_ids){
  pre <- FIRST_OUT:(g-1)
  s <- dta_drop[dta_drop$year %in% pre & !is.na(dta_drop$any.fatalities) &
                dta_drop$agency.id %in% c(treated_ids, comp_ids), ]
  s$D <- as.integer(s$agency.id %in% treated_ids); s }

# ---- CR1 clustered vcov (manual, for the slope WCR) --------------------------
cr1_vcov <- function(X, resid, cl, XtXinv){
  G <- unique(cl); p <- ncol(X); meat <- matrix(0,p,p)
  for (g in G){ ix<-cl==g; Xg<-X[ix,,drop=FALSE]; sc<-crossprod(Xg, resid[ix]); meat<-meat+tcrossprod(sc) }
  nG<-length(G); n<-nrow(X); adj<-(nG/(nG-1))*((n-1)/(n-p)); XtXinv %*% (adj*meat) %*% XtXinv }

# ---- restricted wild cluster bootstrap-t for the linear slope ----------------
wcr_slope <- function(s, B=999, seed=1){
  X <- model.matrix(~ year*D, s); y <- s$any.fatalities; cl <- s$agency.num
  k <- which(colnames(X)=="year:D"); XtXinv <- solve(crossprod(X)); bhat <- XtXinv %*% crossprod(X,y)
  V <- cr1_vcov(X, y - X%*%bhat, cl, XtXinv); t_obs <- bhat[k]/sqrt(V[k,k])
  Xr <- X[,-k,drop=FALSE]; br <- solve(crossprod(Xr)) %*% crossprod(Xr,y); yhat_r <- Xr%*%br; ur <- y - yhat_r
  G <- unique(cl); set.seed(seed); tb <- numeric(B)
  for (b in 1:B){ wg <- sample(WEBB, length(G), replace=TRUE); names(wg) <- G
    ystar <- yhat_r + wg[as.character(cl)]*ur
    bst <- XtXinv %*% crossprod(X, ystar); Vs <- cr1_vcov(X, ystar - X%*%bst, cl, XtXinv)
    tb[b] <- bst[k] / sqrt(Vs[k,k]) }                       # restricted null -> centered at 0
  (1 + sum(abs(tb) >= abs(t_obs))) / (1 + B) }

# ---- restricted wild cluster bootstrap for the joint pre-trend Wald ----------
wcr_joint <- function(s, g, B=499, seed=1){
  mu <- feols(any.fatalities ~ i(year, D, ref=g-1) | agency.num + year, s, cluster=~agency.num)
  W_obs <- unname(fixest::wald(mu, "year::", print=FALSE)$stat)
  mr <- feols(any.fatalities ~ 1 | agency.num + year, s)     # null: all interactions 0
  yhat_r <- predict(mr); ur <- resid(mr); cl <- s$agency.num; G <- unique(cl)
  set.seed(seed); Wb <- numeric(B)
  for (b in 1:B){ wg <- sample(WEBB, length(G), replace=TRUE); names(wg) <- G
    s$ystar <- yhat_r + wg[as.character(cl)]*ur
    mb <- tryCatch(feols(ystar ~ i(year, D, ref=g-1) | agency.num + year, s, cluster=~agency.num), error=function(e) NULL)
    Wb[b] <- if (is.null(mb)) NA else tryCatch(unname(fixest::wald(mb, "year::", print=FALSE)$stat), error=function(e) NA) }
  (1 + sum(Wb >= W_obs, na.rm=TRUE)) / (1 + sum(!is.na(Wb))) }

# ---- randomization inference null over never-treated placebo cohorts ---------
ri_null <- function(g, n_g, B=999, seed=2){
  pre <- FIRST_OUT:(g-1)
  base <- dta_drop[dta_drop$year %in% pre & !is.na(dta_drop$any.fatalities) & dta_drop$agency.id %in% never_ids, ]
  ids <- unique(base$agency.id); set.seed(seed); sl <- numeric(B); wd <- numeric(B)
  for (b in 1:B){ tr <- sample(ids, n_g); d <- base; d$D <- as.integer(d$agency.id %in% tr)
    ml <- tryCatch(feols(any.fatalities ~ year*D, d), error=function(e) NULL)
    sl[b] <- if (is.null(ml)) NA else unname(coef(ml)["year:D"])
    mj <- tryCatch(feols(any.fatalities ~ i(year, D, ref=g-1) | agency.num + year, d), error=function(e) NULL)
    wd[b] <- if (is.null(mj)) NA else tryCatch(unname(fixest::wald(mj, "year::", print=FALSE)$stat), error=function(e) NA) }
  list(slope=sl, wald=wd) }

# ============================================================================
# ITEM 1: side-by-side CR1 / WCR / RI for 2009 and 2013
# ============================================================================
cr1 <- read.csv(file.path(OUT,"d6_pretests_2009_2013.csv"))
# distinct comparison groups (Cengiz == CS: never-treated)
grp <- list(`CS/Cengiz (never-treated)` = function(g) never_ids,
            `stacked (No-Change)`       = function(g) nochange_ids,
            `imputation (never+not-yet)`= function(g) union(never_ids, notyet_ids(g)))
cr1_key <- c(`CS/Cengiz (never-treated)`="CS (never-treated)",
             `stacked (No-Change)`="stacked (No-Change controls)",
             `imputation (never+not-yet)`="imputation (never+not-yet)")

rows <- list()
for (g in c(2009, 2013)){
  n_g <- length(unique(dta_drop$agency.id[dta_drop$gname==g]))
  tids <- unique(dta_drop$agency.id[dta_drop$gname==g])
  cat(sprintf("\n-- cohort %d (n_g=%d): running WCR + RI --\n", g, n_g))
  ri <- ri_null(g, n_g)                                           # one placebo null per cohort
  for (gl in names(grp)){
    s <- pre_sub(g, tids, grp[[gl]](g))
    ml <- feols(any.fatalities ~ year*D, s, cluster=~agency.num)
    s_obs <- unname(coef(ml)["year:D"])
    W_obs <- unname(fixest::wald(feols(any.fatalities ~ i(year,D,ref=g-1)|agency.num+year, s, cluster=~agency.num), "year::", print=FALSE)$stat)
    cr <- cr1[cr1$cohort==g & cr1$estimator==cr1_key[gl], ]
    rows[[length(rows)+1L]] <- data.frame(cohort=g, comparison=gl,
      slope=s_obs,
      slope_p_CR1 = cr$slope_p, slope_p_WCR = wcr_slope(s), slope_p_RI = (1+sum(abs(ri$slope)>=abs(s_obs),na.rm=TRUE))/(1+sum(!is.na(ri$slope))),
      joint_p_CR1 = cr$joint_p, joint_p_WCR = wcr_joint(s, g), joint_p_RI = (1+sum(ri$wald>=W_obs,na.rm=TRUE))/(1+sum(!is.na(ri$wald))),
      row.names=NULL) }
}
PT3 <- do.call(rbind, rows)
write.csv(PT3, file.path(OUT,"d6_pretests_2009_2013_robust.csv"), row.names=FALSE)
cat("\n=== ITEM 1: pre-test p-values  CR1 / WCR / RI ===\n")
print(transform(PT3, slope=round(slope,4),
  slope_p_CR1=round(slope_p_CR1,3), slope_p_WCR=round(slope_p_WCR,3), slope_p_RI=round(slope_p_RI,3),
  joint_p_CR1=round(joint_p_CR1,3), joint_p_WCR=round(joint_p_WCR,3), joint_p_RI=round(joint_p_RI,3)), row.names=FALSE)

# ============================================================================
# ITEM 2: predicted vs actual CS-imputation gap, per droppers cohort
# ============================================================================
TAB <- read.csv(file.path(OUT,"weight_checkability_table.csv"))
csd <- TAB[TAB$estimator=="CS" & TAB$sample=="droppers", c("cohort","pre_g","atom_or_stack_estimate")]
names(csd)[3] <- "tau_CS"
imp <- TAB[TAB$estimator=="imputation" & TAB$sample=="droppers", c("cohort","atom_or_stack_estimate")]
names(imp)[2] <- "tau_imp"
g2 <- merge(csd, imp, by="cohort")
# slope_g vs never-treated, pre-periods only
g2$slope_g <- sapply(g2$cohort, function(g){ tids <- unique(dta_drop$agency.id[dta_drop$gname==g])
  s <- pre_sub(g, tids, never_ids)
  if (length(unique(s$year))<3 || length(unique(s$D))<2) return(NA_real_)
  m <- tryCatch(feols(any.fatalities ~ year*D, s), error=function(e) NULL)
  if (is.null(m) || !("year:D" %in% names(coef(m)))) NA_real_ else unname(coef(m)["year:D"]) })
g2$predicted_gap <- g2$slope_g * (g2$pre_g - 1)/2
g2$actual_gap <- g2$tau_imp - g2$tau_CS
write.csv(g2, file.path(OUT,"d6_gap_predicted_vs_actual.csv"), row.names=FALSE)
cat("\n=== ITEM 2: predicted vs actual CS-imputation gap (droppers) ===\n")
print(transform(g2, tau_CS=round(tau_CS,3), tau_imp=round(tau_imp,3), slope_g=round(slope_g,4),
  predicted_gap=round(predicted_gap,3), actual_gap=round(actual_gap,3)), row.names=FALSE)
fin <- g2[is.finite(g2$predicted_gap) & is.finite(g2$actual_gap), ]
cat(sprintf("correlation(predicted, actual) = %.3f (n=%d cohorts)\n", cor(fin$predicted_gap, fin$actual_gap), nrow(fin)))

rng <- range(c(fin$predicted_gap, fin$actual_gap))
p <- ggplot(fin, aes(predicted_gap, actual_gap)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey55") +
  geom_hline(yintercept=0, color="grey85") + geom_vline(xintercept=0, color="grey85") +
  geom_point(aes(size=abs(actual_gap)), color="#0072B2", alpha=0.85, show.legend=FALSE) +
  ggrepel::geom_text_repel(aes(label=cohort), size=3, seed=1) +
  coord_equal(xlim=rng, ylim=rng) +
  labs(x=expression("predicted gap  "*slope[g]%*%(pre[g]-1)/2),
       y=expression("actual gap  "*hat(tau)[imp]-hat(tau)[CS]),
       title="Predicted vs. actual CS-imputation per-cohort gap (droppers)",
       subtitle=sprintf("dashed = 45 deg; predicted from pre-trend vs never-treated (r=%.2f, n=%d)",
                        cor(fin$predicted_gap, fin$actual_gap), nrow(fin))) +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"d6_gap_predicted_vs_actual.pdf"), p, width=6.8, height=6, device=cairo_pdf)
ggsave(file.path(FIG,"d6_gap_predicted_vs_actual.png"), p, width=6.8, height=6, dpi=200)
file.copy(file.path(FIG,"d6_gap_predicted_vs_actual.png"), file.path(OUT,"d6_gap_predicted_vs_actual.png"), overwrite=TRUE)

cat("\nWrote d6_pretests_2009_2013_robust.csv, d6_gap_predicted_vs_actual.csv, d6_gap_predicted_vs_actual.*\n")
