# =============================================================================
# taskB4_bootstrap_vcov.R
# CS's apparatus is built on the multiplier bootstrap, not pointwise analytical
# inference. Re-run the per-cohort pre-trend Wald on the multiplier-bootstrap
# vcov and compare to the analytical vcov, side by side, under BOTH base-period
# conventions (varying = adjacent-period changes; universal = long differences
# from g-1, base cell dropped as a structural zero). Report whether they agree
# on which cohorts reject.
# =============================================================================
suppressMessages({library(did); library(Matrix)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
dta <- read.csv(file.path(DATA_DIR,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$yc <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta$yc <- ifelse(is.na(dta$yc)&dta$no.req==0,0, ifelse(is.na(dta$yc)&dta$no.req==1,1987,dta$yc))

waldp <- function(a, V) { a<-as.numeric(a); n<-length(a); if(!n) return(NA_real_)
  V<-as.matrix(V); sv<-svd(V)$d; if (sum(sv>sqrt(.Machine$double.eps)*max(sv))<n) return(NA_real_)
  pchisq(as.numeric(t(a)%*%solve(V)%*%a), n, lower.tail=FALSE) }

run_bp <- function(bp) {
  set.seed(0)
  mp <- did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="yc",
                    xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period=bp,data=dta)
  n<-mp$DIDparams$n; att<-mp$att; G<-mp$group; Tt<-mp$t
  Va <- as.matrix(mp$V_analytical)/n
  IF <- as.matrix(mp$inffunc)
  set.seed(7); Bm<-2000
  boot <- matrix(0, Bm, length(att))
  for (b in 1:Bm) { w<-sample(c(-1,1), nrow(IF), replace=TRUE); boot[b,]<-colSums(IF*w)/n }
  Vb <- cov(boot)
  groups <- sort(unique(G[G>0])); ntr<-sapply(groups,function(g) length(unique(dta$agency.id[dta$yc==g])))
  do.call(rbind, lapply(seq_along(groups), function(i){ g<-groups[i]
    idx <- if (bp=="universal") which(G==g & Tt<g & Tt<(g-1) & !is.na(att)) else which(G==g & Tt<g & !is.na(att))
    data.frame(base=bp, cohort=g, n_pre=length(idx), n_treated=ntr[i],
      p_analytical=waldp(att[idx], Va[idx,idx,drop=FALSE]),
      p_bootstrap =waldp(att[idx], Vb[idx,idx,drop=FALSE])) }))
}
res <- rbind(run_bp("varying"), run_bp("universal"))
res$reject_anal <- ifelse(!is.na(res$p_analytical) & res$p_analytical<0.05,1,0)
res$reject_boot <- ifelse(!is.na(res$p_bootstrap)  & res$p_bootstrap <0.05,1,0)
res$agree <- ifelse(is.na(res$p_analytical)|is.na(res$p_bootstrap), NA, res$reject_anal==res$reject_boot)
options(width=180)
cat("== B4: per-cohort pre-trend p, analytical vs multiplier-bootstrap vcov, both base periods ==\n")
print(transform(res, p_analytical=signif(p_analytical,3), p_bootstrap=signif(p_bootstrap,3)), row.names=FALSE)
for (bp in c("varying","universal")) { s<-res[res$base==bp & !is.na(res$agree),]
  cat(sprintf("\n[%s] analytical vs bootstrap agree on reject/not for %d of %d testable cohorts; both reject: %d\n",
      bp, sum(s$agree), nrow(s), sum(s$reject_anal==1 & s$reject_boot==1))) }
write.csv(res, file.path(OUT_DIR,"taskB4_bootstrap_vcov.csv"), row.names=FALSE)
cat("\nWrote output/taskB4_bootstrap_vcov.csv\n")
