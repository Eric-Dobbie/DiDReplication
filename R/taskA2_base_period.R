# =============================================================================
# taskA2_base_period.R
#
# Are the CS and SA PRE-treatment estimates the same object? The R5 atom
# equivalence covers POST-treatment atoms; pre-treatment placebos depend on the
# base-period convention. Extract, per cohort and pre-period event time, the
# point estimate from:
#   * CS as currently run  (base_period = "varying")
#   * CS universal base     (base_period = "universal")
#   * SA (sunab) leads vs reference period -1
# Report side by side with differences and the max abs diff per cohort, and the
# two vcovs' SEs + ratio for the matched coefficients.
# =============================================================================
suppressMessages({library(did); library(fixest)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset = "data")
OUT_DIR  <- Sys.getenv("DIDREP_OUT",  unset = "output")

dta <- read.csv(file.path(DATA_DIR, "dta.csv"))
dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id,
  FUN = function(x){v <- x[!is.na(x)]; if (length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed) & dta$no.req == 0, 0,
                    ifelse(is.na(dta$year.changed) & dta$no.req == 1, 1987, dta$year.changed))

cs_pre <- function(bp) {
  set.seed(0)
  mp <- did::att_gt(yname="any.fatalities", tname="year", idname="agency.num",
                    gname="year.changed", xformla=~1, control_group="nevertreated",
                    clustervars="agency.num", base_period=bp, data=dta)
  k <- mp$t < mp$group & !is.na(mp$att)
  n <- mp$DIDparams$n
  data.frame(cohort=mp$group[k], e=mp$t[k]-mp$group[k],
             est=mp$att[k], se=sqrt(pmax(diag(as.matrix(mp$V_analytical))[k],0)/n))  # correct: /n
}
csv_ <- cs_pre("varying"); csu  <- cs_pre("universal")

# SA leads
d_sa <- dta[!is.na(dta$any.fatalities) & (dta$year.changed==0 | dta$year.changed>=2001),]
d_sa$coh <- ifelse(d_sa$year.changed==0, 10000, d_sa$year.changed)
res <- fixest::feols(any.fatalities ~ sunab(coh, year) | agency.num + year, data=d_sa, cluster=~agency.num)
ct <- summary(res, agg=FALSE)$coeftable
mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee <- as.integer(vapply(mm, function(z) z[2], character(1)))
gg <- as.integer(vapply(mm, function(z) z[3], character(1)))
k2 <- !is.na(ee) & !is.na(gg) & ee < 0
sa <- data.frame(cohort=gg[k2], e=ee[k2], est=ct[k2,1], se=ct[k2,2])

m <- merge(merge(csv_, csu, by=c("cohort","e"), suffixes=c(".csV",".csU")),
           sa, by=c("cohort","e")); names(m)[names(m)=="est"]<-"est.SA"; names(m)[names(m)=="se"]<-"se.SA"
m <- m[order(m$cohort, m$e), ]
m$d_csV_SA <- m$est.csV - m$est.SA          # CS-varying vs SA
m$d_csU_SA <- m$est.csU - m$est.SA          # CS-universal vs SA
m$se_ratio_csU_SA <- m$se.csU / m$se.SA      # analytical vs cluster-robust SE

options(width=220)
cat("== Pre-treatment point estimates: CS-varying, CS-universal, SA ==\n")
print(transform(m, est.csV=round(est.csV,4), est.csU=round(est.csU,4), est.SA=round(est.SA,4),
                d_csV_SA=round(d_csV_SA,4), d_csU_SA=round(d_csU_SA,4),
                se.csU=round(se.csU,4), se.SA=round(se.SA,4), se_ratio_csU_SA=round(se_ratio_csU_SA,3))[
      ,c("cohort","e","est.csV","est.csU","est.SA","d_csV_SA","d_csU_SA","se.csU","se.SA","se_ratio_csU_SA")],
      row.names=FALSE)

cat("\n== max |difference| per cohort ==\n")
agg <- aggregate(cbind(maxabs_csV_SA=abs(d_csV_SA), maxabs_csU_SA=abs(d_csU_SA)) ~ cohort, data=m, FUN=max)
print(transform(agg, maxabs_csV_SA=signif(maxabs_csV_SA,3), maxabs_csU_SA=signif(maxabs_csU_SA,3)), row.names=FALSE)
cat(sprintf("\nOVERALL max|CS_universal - SA| = %.2e ; max|CS_varying - SA| = %.3f\n",
            max(abs(m$d_csU_SA)), max(abs(m$d_csV_SA))))
cat(sprintf("median SE ratio (CS analytical / SA cluster-robust) = %.3f\n", median(m$se_ratio_csU_SA)))
write.csv(m, file.path(OUT_DIR, "taskA2_base_period.csv"), row.names=FALSE)
cat("\nWrote output/taskA2_base_period.csv\n")
