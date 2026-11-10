# =============================================================================
# taskA4_pretrend_magnitude.R
#
# Rebuild the pre-trend figure on MAGNITUDE, not p-values (A1 showed the
# p-values were a V_analytical scaling artifact; A4 asks for coefficient size
# with its CI so that "can't be checked" is visually distinct from "checked and
# flat"). Per cohort, per estimator, summarize the pre-treatment leads by the
# LARGEST ABSOLUTE pre-treatment coefficient and its 95% CI:
#   y = that coefficient +/- 1.96 SE ; x = aggregation weight ;
#   size = number of pre-treatment periods available ; label = cohort ;
#   facet = estimator. Correct scaling throughout: CS uses base_period="universal"
#   (event-study leads vs g-1, excluding the mechanically-zero base) with
#   se=sqrt(diag(V_analytical)/n); SA = sunab leads; stacked = per-cohort event
#   study on the shipped stack (+/-4 window). n_treated annotated.
# Writes output/taskA4_pretrend_magnitude.csv, figures/pretrend_magnitude.{pdf,png}.
# =============================================================================
suppressMessages({library(did); library(fixest); library(ggplot2); library(ggrepel)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
FIG_DIR <- Sys.getenv("DIDREP_FIG", unset="figures")

dta <- read.csv(file.path(DATA_DIR,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed)&dta$no.req==0,0, ifelse(is.na(dta$year.changed)&dta$no.req==1,1987,dta$year.changed))
ci <- read.csv(file.path(OUT_DIR,"task1_cohort_inputs.csv"))
WCS <- setNames(ci$WCS_restr, ci$cohort); wST <- setNames(ci$WStacked_restr, ci$cohort)
ntr <- setNames(sapply(as.integer(names(WCS)), function(g) length(unique(dta$agency.id[dta$year.changed==g]))), names(WCS))

pick <- function(est, se, cohort, e) {            # largest |coef| lead per cohort
  out <- lapply(sort(unique(cohort)), function(g){
    k <- cohort==g; if(!any(k)) return(NULL)
    i <- which(k)[which.max(abs(est[k]))]
    data.frame(cohort=g, n_pre=sum(k), coef=est[i], se=se[i]) })
  do.call(rbind, out)
}

# CS universal (correct /n), leads e<=-2 (exclude base e=-1)
set.seed(0)
mp <- did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="year.changed",
                  xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period="universal",data=dta)
n <- mp$DIDparams$n; seA <- sqrt(diag(as.matrix(mp$V_analytical))/n)
e <- mp$t-mp$group; kk <- e <= -2 & !is.na(mp$att)
cs <- pick(mp$att[kk], seA[kk], mp$group[kk], e[kk]); cs$estimator <- "CS (universal)"

# SA leads
d_sa <- dta[!is.na(dta$any.fatalities)&(dta$year.changed==0|dta$year.changed>=2001),]; d_sa$coh <- ifelse(d_sa$year.changed==0,10000,d_sa$year.changed)
res <- fixest::feols(any.fatalities ~ sunab(coh, year)|agency.num+year, data=d_sa, cluster=~agency.num)
ct <- summary(res, agg=FALSE)$coeftable
mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee <- as.integer(vapply(mm,function(z)z[2],character(1))); gg <- as.integer(vapply(mm,function(z)z[3],character(1)))
k2 <- !is.na(ee)&!is.na(gg)&ee<0; sa <- pick(ct[k2,1], ct[k2,2], gg[k2], ee[k2]); sa$estimator <- "SA"

# stacked per-cohort event study (+/-4), leads
stk <- read.csv(file.path(DATA_DIR,"stacked_fatal.csv")); stk$etime <- stk$year-stk$cohort
strow <- list()
for (g in as.integer(names(WCS))) {
  s <- stk[stk$cohort==g & !is.na(stk$any.fatalities) & abs(stk$etime)<=4,]
  if (length(unique(s$agency.id[s$treat==0]))==0) next           # no controls -> no clean leads
  m <- tryCatch(fixest::feols(any.fatalities ~ i(etime, treat, ref=-1)|agency.id+year, data=s, warn=FALSE, notes=FALSE), error=function(e) NULL)
  if (is.null(m)) next
  c2 <- fixest::coeftable(m); rn <- grep("etime::-", rownames(c2), value=TRUE); if(!length(rn)) next
  i <- rn[which.max(abs(c2[rn,1]))]
  strow[[length(strow)+1]] <- data.frame(cohort=g, n_pre=length(rn), coef=c2[i,1], se=c2[i,2], estimator="stacked")
}
st <- do.call(rbind, strow)

all <- rbind(cs, sa, st)
all$weight <- ifelse(all$estimator=="stacked", wST[as.character(all$cohort)], WCS[as.character(all$cohort)])
all$n_treated <- ntr[as.character(all$cohort)]
all$lo <- all$coef-1.96*all$se; all$hi <- all$coef+1.96*all$se
all <- all[!is.na(all$weight),]
write.csv(all, file.path(OUT_DIR,"taskA4_pretrend_magnitude.csv"), row.names=FALSE)
options(width=160); print(transform(all, coef=round(coef,3), se=round(se,3), weight=round(weight,3)), row.names=FALSE)

pal <- c("CS (universal)"="#0072B2","SA"="#D55E00","stacked"="#009E73")
g <- ggplot(all, aes(weight, coef, color=estimator)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey55", linewidth=0.4) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0, linewidth=0.5, alpha=0.8) +
  geom_point(aes(size=n_pre)) +
  ggrepel::geom_text_repel(aes(label=sprintf("%d (n=%d)", cohort, n_treated)), size=2.5, show.legend=FALSE, max.overlaps=20, seed=7) +
  facet_wrap(~estimator, ncol=3) +
  scale_color_manual(values=pal, guide="none") +
  scale_size_continuous(range=c(1.5,5), name="pre-periods") +
  labs(x="Aggregation weight on cohort", y="Largest |pre-treatment coefficient| with 95% CI") +
  theme_minimal(base_size=12) + theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(FIG_DIR,"pretrend_magnitude.pdf"), g, width=11, height=5, device=cairo_pdf)
ggsave(file.path(FIG_DIR,"pretrend_magnitude.png"), g, width=11, height=5, dpi=200)
cat("\nWrote output/taskA4_pretrend_magnitude.csv and figures/pretrend_magnitude.{pdf,png}\n")
