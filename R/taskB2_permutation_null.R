# =============================================================================
# taskB2_permutation_null.R  (DECISIVE)
#
# Does the per-cohort pre-trend test hold size under the null? Permutation test
# on the never-treated agencies only, so NO real treatment effect exists by
# construction. Assign placebo adoption years reproducing the observed cohort-
# size structure and adoption-year spacing to a subset of the 41 never-treated
# agencies; the rest are controls; run the IDENTICAL per-cohort pre-trend test
# (att_gt, corrected V_analytical/n Wald). Repeat B times. A well-behaved test
# rejects 5% at alpha=0.05 and its p-values are uniform on [0,1].
#
# Writes output/taskB2_null_pvalues.csv (perm x cohort), a size-by-cohort-size
# summary, and figures/taskB2_null_pvalue_hist.{pdf,png}.
# =============================================================================
suppressMessages({library(did); library(ggplot2)})
DATA_DIR <- Sys.getenv("DIDREP_DATA", unset="data"); OUT_DIR <- Sys.getenv("DIDREP_OUT", unset="output")
FIG_DIR <- Sys.getenv("DIDREP_FIG", unset="figures")
B <- as.integer(Sys.getenv("B2_B", unset="500"))

dta <- read.csv(file.path(DATA_DIR,"dta.csv")); dta$agency.num <- as.numeric(factor(dta$agency.id))
yc <- ave(dta$year.changed, dta$agency.id, FUN=function(x){v<-x[!is.na(x)]; if(length(v)) min(v) else NA_real_})
never_ids <- unique(dta$agency.id[is.na(yc) & dta$no.req==0])     # 41 never-treated
nonchg_ids <- unique(dta$agency.id[dta$change.type=="No Change"]) # 741 non-changers (no real effect)
cat("never-treated:", length(never_ids), " | non-changers:", length(nonchg_ids), "\n")
cat("DESIGN: each permutation draws 41 controls + 28 placebo-treated from the 741 non-changers\n")
cat("        (matches the real control count, keeps all cohorts full-rank; pure null by random assignment).\n")

# observed cohort structure (identified cohorts)
coh_years <- c(2002,2004,2008,2009,2012,2013,2014,2017,2018,2019,2020)
coh_size  <- c(5,2,1,7,1,3,1,3,3,1,1)
N_CTRL <- 41
cat("placebo-treated total:", sum(coh_size), " controls:", N_CTRL, "\n\n")

nd_all <- dta[dta$agency.id %in% nonchg_ids, c("agency.id","agency.num","year","any.fatalities","no.req")]
nd_by <- split(nd_all, nd_all$agency.id)

wald_p <- function(a, V) { a<-as.numeric(a); n<-length(a); if(!n) return(c(NA,n))
  V<-as.matrix(V); sv<-svd(V)$d; if (sum(sv>sqrt(.Machine$double.eps)*max(sv))<n) return(c(NA,n))
  c(pchisq(as.numeric(t(a)%*%solve(V)%*%a), n, lower.tail=FALSE), n) }

run_once <- function() {
  pick <- sample(nonchg_ids, N_CTRL + sum(coh_size))         # distinct agencies
  trt <- pick[1:sum(coh_size)]; ctl <- pick[(sum(coh_size)+1):length(pick)]
  g <- setNames(rep(0, length(pick)), pick); pos <- 1
  for (j in seq_along(coh_years)) { ids <- trt[pos:(pos+coh_size[j]-1)]; g[ids] <- coh_years[j]; pos <- pos+coh_size[j] }
  d <- do.call(rbind, nd_by[pick]); d$gph <- g[as.character(d$agency.id)]
  mp <- tryCatch(did::att_gt(yname="any.fatalities",tname="year",idname="agency.num",gname="gph",
          xformla=~1,control_group="nevertreated",clustervars="agency.num",base_period="varying",data=d),
          error=function(e) NULL)
  if (is.null(mp)) return(setNames(rep(NA_real_,length(coh_years)), coh_years))
  V <- as.matrix(mp$V_analytical)/mp$DIDparams$n; att<-mp$att; G<-mp$group; Tt<-mp$t
  sapply(coh_years, function(cg){ idx<-which(G==cg & Tt<cg & !is.na(att))
    if(!length(idx)) return(NA_real_); wald_p(att[idx], V[idx,idx,drop=FALSE])[1] })
}

set.seed(2024)
P <- matrix(NA_real_, B, length(coh_years), dimnames=list(NULL, as.character(coh_years)))
for (b in 1:B) P[b,] <- run_once()

# ---- size by cohort ---------------------------------------------------------
npre <- coh_years - 2000 - 1 + 1   # approx leads; recompute exactly from one fit below
# exact n_pre per cohort = number of pre periods (t in 2000..cg-1)
npre <- pmin(coh_years - 2000, coh_years - 2000)   # placeholder; set precisely:
npre <- sapply(coh_years, function(cg) length(2000:(cg-1)) - 0)   # pre years count (varying drops none here)
size_by <- data.frame(cohort=coh_years, placebo_size=coh_size, n_pre_approx=npre,
  rej_rate_05 = round(colMeans(P < 0.05, na.rm=TRUE),3),
  n_valid = colSums(!is.na(P)))
options(width=170); cat("== B2: empirical rejection rate at alpha=0.05, per placebo cohort ==\n")
print(size_by, row.names=FALSE)
cat(sprintf("\nOVERALL rejection rate (pooled cells) = %.3f  (nominal 0.05)\n", mean(P < 0.05, na.rm=TRUE)))
cat("rejection rate by placebo cohort size:\n")
agg <- aggregate(rej ~ size, data=data.frame(size=rep(coh_size, each=B), rej=as.vector(P<0.05)), FUN=function(x) mean(x,na.rm=TRUE))
print(transform(agg, rej=round(rej,3)), row.names=FALSE)
# rejection rate vs n_pre
rr_npre <- data.frame(n_pre=npre, rej=round(colMeans(P<0.05,na.rm=TRUE),3))
cat("\nrejection rate vs number of pre-periods:\n"); print(rr_npre[order(rr_npre$n_pre),], row.names=FALSE)

write.csv(as.data.frame(P), file.path(OUT_DIR,"taskB2_null_pvalues.csv"), row.names=FALSE)
write.csv(size_by, file.path(OUT_DIR,"taskB2_size_by_cohort.csv"), row.names=FALSE)

# ---- null distribution plot -------------------------------------------------
pv <- data.frame(p = as.vector(P)); pv <- pv[!is.na(pv$p),,drop=FALSE]
g1 <- ggplot(pv, aes(p)) + geom_histogram(breaks=seq(0,1,0.05), fill="#0072B2", color="white") +
  geom_hline(yintercept=nrow(pv)/20, linetype="dashed", color="grey40") +
  labs(x="pre-trend p-value under the null", y="count") +
  theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG_DIR,"taskB2_null_pvalue_hist.pdf"), g1, width=7, height=4.5, device=cairo_pdf)
ggsave(file.path(FIG_DIR,"taskB2_null_pvalue_hist.png"), g1, width=7, height=4.5, dpi=200)
cat(sprintf("\nUnder a well-behaved test the histogram is flat at %.0f per bin; KS vs uniform p = %.3g\n",
            nrow(pv)/20, suppressWarnings(ks.test(pv$p,"punif")$p.value)))
cat("Wrote taskB2_null_pvalues.csv, taskB2_size_by_cohort.csv, figures/taskB2_null_pvalue_hist.*\n")
