# =============================================================================
# task4_pretrend_weight.R  (Task 4.2)
#
# Joint distribution of aggregation weight and pre-trend support, per cohort, per
# estimator. For each estimator we compute a per-cohort joint pre-trend test
# (Wald on that cohort's pre-treatment leads = 0) and pair it with the cohort's
# aggregation weight, so heavily weighted but poorly supported atoms are visible:
#   * CS   : pre-cells ATT(g,t), t<g, from att_gt + mp$V_analytical; weight = CS
#            simple cohort weight.
#   * SA   : pre-cells CATT(g,e), e<0, from sunab + vcov; weight = SA att weight.
#   * stacked : from output/stacked_pretrends.csv (shipped stack); weight = FWL w_s.
# Writes output/task4_pretrend_weight.csv and figures/pretrend_vs_weight.{pdf,png}
# (x = aggregation weight, y = pre-trend p-value, one series per estimator).
# =============================================================================
suppressMessages({library(did); library(fixest); library(ggplot2); library(MASS); library(ggrepel)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")
FIG_DIR  <- Sys.getenv("DIDREP_FIG",  unset = "figures")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987, dta$year.changed))

wald <- function(a, V) {               # joint Wald with generalized inverse
  a <- as.numeric(a); if (!length(a)) return(c(NA, NA))
  Vi <- tryCatch(MASS::ginv(V), error = function(e) NULL); if (is.null(Vi)) return(c(NA, NA))
  stat <- as.numeric(t(a) %*% Vi %*% a); df <- qr(V)$rank
  c(stat, pchisq(stat, df = max(df,1), lower.tail = FALSE))
}

# ---- CS per-cohort pre-trend + simple weight --------------------------------
set.seed(0)
mp <- did::att_gt(yname="any.fatalities", tname="year", idname="agency.num",
                  gname="year.changed", xformla=~1, control_group="nevertreated",
                  clustervars="agency.num", data=dta)
V <- as.matrix(mp$V_analytical); att <- mp$att; G <- mp$group; Tt <- mp$t
dp <- mp$DIDparams; idata <- dp$data; keepu <- !duplicated(idata[[dp$idname]])
ug <- idata[[dp$gname]][keepu]
groups <- sort(unique(G[G>0]))
pg <- sapply(groups, function(g) mean(ug==g))
# simple cohort weight
ispost <- Tt >= G; wcell <- ifelse(ispost & !is.na(att), pg[match(G,groups)], 0)
wcell[is.na(wcell)] <- 0; wCS <- tapply(wcell, G, sum)/sum(wcell)
ntr <- sapply(groups, function(g) length(unique(dta$agency.id[dta$year.changed==g])))
csrows <- lapply(seq_along(groups), function(i){ g <- groups[i]
  idx <- which(G==g & Tt < g & !is.na(att))
  w <- if (length(idx)>=1) wald(att[idx], V[idx,idx,drop=FALSE]) else c(NA,NA)
  data.frame(estimator="CS", cohort=g, weight=as.numeric(wCS[as.character(g)]),
             n_treated=ntr[i], n_leads=length(idx), wald=w[1], p=w[2]) })
cs <- do.call(rbind, csrows)

# ---- SA per-cohort pre-trend + att weight -----------------------------------
d_sa <- dta[!is.na(dta$any.fatalities) & (dta$year.changed==0 | dta$year.changed>=2001),]
d_sa$coh <- ifelse(d_sa$year.changed==0, 10000, d_sa$year.changed)
res <- fixest::feols(any.fatalities ~ sunab(coh, year) | agency.num + year,
                     data=d_sa, cluster=~agency.num)
ct <- summary(res, agg=FALSE)$coeftable; V2 <- vcov(res)[rownames(ct), rownames(ct)]
mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee <- as.integer(vapply(mm, function(z) z[2], character(1)))
gg <- as.integer(vapply(mm, function(z) z[3], character(1)))
keep <- !is.na(ee)&!is.na(gg)
# SA att weights: obs counts per post cell (as in extract_sa_atoms), summed per cohort
mmx <- stats::model.matrix(res); sh <- colSums(abs(sign(mmx)))[rownames(ct)]
wpost <- ifelse(ee>=0 & keep, sh, 0); wpost[is.na(wpost)] <- 0
wSA <- tapply(wpost, gg, sum, na.rm=TRUE); wSA <- wSA/sum(wpost)
sarows <- lapply(sort(unique(gg[keep])), function(g){
  idx <- which(gg==g & ee<0 & keep)
  w <- if (length(idx)>=1) wald(ct[idx,1], V2[idx,idx,drop=FALSE]) else c(NA,NA)
  data.frame(estimator="SA", cohort=g, weight=as.numeric(wSA[as.character(g)]),
             n_treated=ntr[match(g,groups)], n_leads=length(idx), wald=w[1], p=w[2]) })
sa <- do.call(rbind, sarows)

# ---- stacked from committed pretrends + FWL weights -------------------------
pt <- read.csv(file.path(OUT_DIR,"stacked_pretrends.csv"))
pt <- pt[pt$stack=="shipped-11" & pt$reason=="tested",]
fw <- read.csv(file.path(OUT_DIR,"fwl_decomp_unweighted.csv"))
wST <- setNames(fw$w_s/sum(fw$w_s), fw$stack)
st <- data.frame(estimator="stacked", cohort=pt$cohort,
                 weight=as.numeric(wST[as.character(pt$cohort)]),
                 n_treated=pt$n_treat, n_leads=pt$n_leads,
                 wald=pt$wald_stat, p=pt$pretrend_p)

all <- rbind(cs, sa, st); row.names(all) <- NULL
write.csv(all, file.path(OUT_DIR,"task4_pretrend_weight.csv"), row.names=FALSE)
cat("== Task 4.2: pre-trend vs weight (per cohort, per estimator) ==\n")
print(transform(all, weight=round(weight,3), wald=round(wald,2), p=round(p,3)), row.names=FALSE)

# ---- plot -------------------------------------------------------------------
pal <- c(CS="#0072B2", SA="#D55E00", stacked="#009E73")
pd <- all[!is.na(all$p),]
g <- ggplot(pd, aes(weight, p, color=estimator)) +
  geom_hline(yintercept=0.05, linetype="dashed", color="grey55", linewidth=0.4) +
  geom_point(aes(size=n_treated), alpha=0.85) +
  ggrepel::geom_text_repel(aes(label=cohort), size=2.8, show.legend=FALSE, max.overlaps=20, seed=7) +
  scale_color_manual(values=pal, name=NULL) +
  scale_size_continuous(range=c(2,6), name="treated units") +
  labs(x="Aggregation weight on cohort", y="Per-cohort pre-trend p-value") +
  theme_minimal(base_size=12) +
  theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIG_DIR,"pretrend_vs_weight.pdf"), g, width=9, height=6, device=cairo_pdf)
ggsave(file.path(FIG_DIR,"pretrend_vs_weight.png"), g, width=9, height=6, dpi=200)
cat("\nWrote output/task4_pretrend_weight.csv and figures/pretrend_vs_weight.{pdf,png}\n")
