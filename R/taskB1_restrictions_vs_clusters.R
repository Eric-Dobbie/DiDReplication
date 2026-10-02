# =============================================================================
# taskB1_restrictions_vs_clusters.R
#
# Over-rejection diagnostic: joint pre-trend Wald tests impose n_pre restrictions
# estimated off n_treated clusters. With n_pre >> n_treated the Wald test over-
# rejects (Conley-Taber 2011; MacKinnon-Webb 2017). One row per cohort per
# estimator: n_pre, n_treated, n_control, ratio = n_pre/n_treated, and the Wald/
# df/p already computed (corrected V_analytical/n for CS; cluster-robust for SA).
# Reports whether p-values sort with the restriction/cluster ratio.
# =============================================================================
suppressMessages({library(did); library(fixest); library(MASS)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
dta <- read.csv(file.path(DATA_DIR,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
dta$year.changed <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
dta$year.changed <- ifelse(is.na(dta$year.changed)&dta$no.req==0,0, ifelse(is.na(dta$year.changed)&dta$no.req==1,1987,dta$year.changed))

wald_full <- function(a, V) { a<-as.numeric(a); n<-length(a); if(!n) return(c(NA,NA,NA))
  V<-as.matrix(V); sv<-svd(V)$d; rk<-sum(sv>sqrt(.Machine$double.eps)*max(sv))
  if (rk<n) return(c(NA, rk, NA))                       # rank-deficient -> NA by rule
  stat<-as.numeric(t(a)%*%solve(V)%*%a); c(stat, n, pchisq(stat, n, lower.tail=FALSE)) }

n_never <- length(unique(dta$agency.id[dta$year.changed==0]))     # control agencies

# CS (varying, corrected /n)
set.seed(0)
mp <- did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="year.changed",
                  xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period="varying",data=dta)
V <- as.matrix(mp$V_analytical)/mp$DIDparams$n; att<-mp$att; G<-mp$group; Tt<-mp$t
groups <- sort(unique(G[G>0])); ntr <- sapply(groups, function(g) length(unique(dta$agency.id[dta$year.changed==g])))
rows <- list()
for (i in seq_along(groups)) { g<-groups[i]; idx<-which(G==g & Tt<g & !is.na(att))
  w<-wald_full(att[idx], V[idx,idx,drop=FALSE])
  rows[[length(rows)+1]] <- data.frame(estimator="CS", cohort=g, n_pre=length(idx),
    n_treated=ntr[i], n_control=n_never, wald=w[1], df=w[2], p=w[3]) }

# SA
d_sa <- dta[!is.na(dta$any.fatalities)&(dta$year.changed==0|dta$year.changed>=2001),]; d_sa$coh<-ifelse(d_sa$year.changed==0,10000,d_sa$year.changed)
res <- fixest::feols(any.fatalities~sunab(coh,year)|agency.num+year, data=d_sa, cluster=~agency.num)
ct <- summary(res,agg=FALSE)$coeftable; Vsa<-as.matrix(vcov(res))[rownames(ct),rownames(ct)]
mm <- regmatches(rownames(ct), regexec("::(-?[0-9]+):.*::([0-9]+)$", rownames(ct)))
ee<-as.integer(vapply(mm,function(z)z[2],character(1))); gg<-as.integer(vapply(mm,function(z)z[3],character(1))); keep<-!is.na(ee)&!is.na(gg)
for (g in sort(unique(gg[keep]))) { idx<-which(gg==g & ee<0 & keep)
  w<-wald_full(ct[idx,1], Vsa[idx,idx,drop=FALSE])
  rows[[length(rows)+1]] <- data.frame(estimator="SA", cohort=g, n_pre=length(idx),
    n_treated=ntr[match(g,groups)], n_control=n_never, wald=w[1], df=w[2], p=w[3]) }

tab <- do.call(rbind, rows); tab$ratio <- round(tab$n_pre/tab$n_treated,2)
tab$restr_ge_treated <- ifelse(tab$n_pre>=tab$n_treated,"YES","")
tab <- tab[tab$n_pre>0,]
options(width=200)
cat("== B1: restrictions vs treated clusters (sorted by ratio within estimator) ==\n")
for (es in c("CS","SA")) { sub<-tab[tab$estimator==es,]; sub<-sub[order(sub$ratio),]
  cat("\n--",es,"--\n"); print(transform(sub, wald=signif(wald,4), p=signif(p,3))[
    ,c("cohort","n_pre","n_treated","n_control","ratio","restr_ge_treated","wald","df","p")], row.names=FALSE) }
# does p sort with ratio? Spearman on CS testable cohorts
cs <- tab[tab$estimator=="CS" & !is.na(tab$p),]
cat(sprintf("\nCS Spearman(cor) ratio vs p = %.3f  (negative => p smaller where ratio larger = over-rejection signature)\n",
            cor(cs$ratio, cs$p, method="spearman")))
write.csv(tab, file.path(OUT_DIR,"taskB1_restrictions.csv"), row.names=FALSE)
cat("Wrote output/taskB1_restrictions.csv\n")
